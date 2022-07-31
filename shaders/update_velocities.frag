#version 300 es

precision highp float;
precision highp int;
precision highp isampler2D ;
precision highp usampler2D ;

uniform sampler2D positions_texture, velocities_texture, bests_texture;
uniform usampler2D  itinymtState, itinymtMat;


layout (location = 0) out vec4 new_velocity;
layout (location = 1) out uvec4 otinymtState;

in vec2 cc;

uniform float phi_local, phi_global, omega;
uniform vec4 chi[4], global_best[4];

/*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
 * tinymt.glsl  :   a glsl file to be included the in the shaders
 *
 * PROGRAMMER   :   ABOUZAR KABOUDIAN
 * DATE         :   Tue 31 Mar 2020 14:25:32 (EDT)
 * PLACE        :   Chaos Lab @ GaTech, Atlanta, GA
 *@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
 */

// global variables and the macros for the tinymt algorithm --------------
uvec4 tinymtState, tinymtMat ;

#define TINYMT32_SH0    1
#define TINYMT32_SH1    10
#define TINYMT32_SH8    8
#define TINYMT32_MASK   uint(0x7fffffff)
#define MIN_LOOP        8
#define PRE_LOOP        8

#define MAX_INT         4294967295
#define MAX_INT_F       4294967295.

/*========================================================================
 * initTinyMt : initialize the Tiny Marsenne Twister for the shader
 *========================================================================
 */
void tinymtInit(){
    // Initialize the random number states
    tinymtState = texture( itinymtState, cc ) ;
    tinymtMat   = texture( itinymtMat,   cc ) ;
    return ;
}

/*========================================================================
 * finish the tiny marsenne twister algorithm: should be called before the
 * shader return call ;
 *========================================================================
 */
void tinymtFinish(){
    otinymtState = uvec4(tinymtState) ;
}
#define tinymtReturn  tinymtFinish

/*========================================================================
 * Iterate to the next state
 *========================================================================
 */
void tinymtNextState(){
    uint x,y ;
    y = tinymtState.a ;
    x = (tinymtState.r & TINYMT32_MASK)
        ^ tinymtState.g
        ^ tinymtState.b;
    x ^= (x << TINYMT32_SH0);
    y ^= (y >> TINYMT32_SH0) ^ x;
    tinymtState.r = tinymtState.g;
    tinymtState.g = tinymtState.b;
    tinymtState.b = x ^ (y << TINYMT32_SH1);
    tinymtState.a = y;
    tinymtState.g ^= uint(-int(y & uint(1)) & int(tinymtMat.r));
    tinymtState.b ^= uint(-int(y & uint(1)) & int(tinymtMat.g));

    return ;
}

/*========================================================================
 * temper the results
 *========================================================================
 */
uint  tinymtTemper() {
    uint t0, t1;
    t0 = tinymtState.a;
#if defined(LINEARITY_CHECK)
  t1 = tinymtState.r
      ^ (tinymtState.b >> TINYMT32_SH8);
#else
    t1 = tinymtState.r
        + (tinymtState.b >> TINYMT32_SH8);
#endif
    t0 ^= t1;
    t0 ^= uint(-int(t1 & uint(1)) & int(tinymtMat[2]));
    return t0;
}

/*========================================================================
 * Get a random number using the tiny marsenne twister
 *========================================================================
 */
uint tinymtUrand(){
    tinymtNextState() ;
    return tinymtTemper() ;
}

/*========================================================================
 * tinymtFrand : returns a float number between 0-1 using the tiny mt
 * algorithm.
 *========================================================================
 */
float tinymtFrand(){
    return float(tinymtUrand())/MAX_INT_F ;
}
#define tinymtRand tinymtFrand

/*========================================================================
 * returns an integer with the binomial distribution with a Bernouli trail
 * of "npar" attempts with each attempt a chance of success of "p"
 *========================================================================
 */
uint tinymtBinran(float p, uint npar){
    uint isum = 0u ;
    if ( p > 0.1 ){
        for(uint i=0u ; i<npar ; i++){
            if ( tinymtRand() < p ){
                isum++ ;
            }
        }
    }else{ /* for small probablities of p use a geometric distribution */
        isum = uint( floor ( log( 1. - tinymtRand() ) /
                             log( 1. - p            )   ) ) ;
    }
    return isum ;
}
/*

    END TINYMT COPY/PASTE

*/


void main() {
    tinymtInit();
    vec4 position = texture(positions_texture, cc);
    vec4 velocity = texture(velocities_texture, cc);
    vec4 best = texture(bests_texture, cc);

    int idx = 0;
    if (cc.x > 0.5) idx += 1;
    if (cc.y > 0.5) idx += 2;

    vec4 my_chi = chi[idx];
    vec4 my_global_best = global_best[idx];

    vec4 r_local = vec4(tinymtRand(), tinymtRand(), tinymtRand(), tinymtRand());
    vec4 r_global = vec4(tinymtRand(), tinymtRand(), tinymtRand(), tinymtRand());

    new_velocity = my_chi * (omega * velocity + phi_local * r_local * (best - position) + phi_global * r_global * (my_global_best - position));

    tinymtReturn();
}
