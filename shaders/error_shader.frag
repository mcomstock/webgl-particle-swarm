#version 300 es
precision highp float ;
precision highp int ;

uniform sampler2D error_tex;

out vec4 outcolor ; /*  output of the shader
                        pixel color                 */
in vec2 cc ;        /*  input from vertex shader    */

// Main body of the shader
void main() {
	vec4 error = texture(error_tex,cc);

	outcolor = vec4(error.r, 0., 0., 1.);
    // if ( length(cc-vec2(0.5,0.5)) < 0.2){
    //     outcolor = vec4(cc.x,0.,0.,1.) ;
    // }else{
    //     outcolor = vec4(0.,0.,0.,1.) ;
    // }
    return ;
}