#version 300 es

precision highp float;
precision highp int;

uniform sampler2D velocities_texture;

layout (location = 0) out vec4 new_velocity;

in vec2 cc;

uniform vec4 lower_bounds, upper_bounds;

#define D0_LOWER_BOUND lower_bounds.r
#define D1_LOWER_BOUND lower_bounds.g
#define D2_LOWER_BOUND lower_bounds.b
#define D3_LOWER_BOUND lower_bounds.a

#define D0_UPPER_BOUND upper_bounds.r
#define D1_UPPER_BOUND upper_bounds.g
#define D2_UPPER_BOUND upper_bounds.b
#define D3_UPPER_BOUND upper_bounds.a

#define D0_VEL velocities.r
#define D1_VEL velocities.g
#define D2_VEL velocities.b
#define D3_VEL velocities.a

float bounds_check(float p_val, float p_min, float p_max)
{
    float retVal;
    if(p_val > p_max)
    {
        // retVal = p_min + (0.75 * (p_max-p_min)) * tinymtRand();
        retVal = p_max;
    }
    else if(p_val < p_min)
    {
        // retVal = p_min + (0.25 * (p_max-p_min)) + (0.75 * (p_max-p_min)) * tinymtRand();
        retVal = p_min;
    }
    else
    {
        retVal = p_val;
    }

    return retVal;
}


void main() {

    vec4 velocities = texture(velocities_texture, cc);

    float new_v0 = bounds_check(D0_VEL, D0_LOWER_BOUND, D0_UPPER_BOUND);
    float new_v1 = bounds_check(D1_VEL, D1_LOWER_BOUND, D1_UPPER_BOUND);
    float new_v2 = bounds_check(D2_VEL, D2_LOWER_BOUND, D2_UPPER_BOUND);
    float new_v3 = bounds_check(D3_VEL, D3_LOWER_BOUND, D3_UPPER_BOUND);


    new_velocity = vec4(new_v0, new_v1, new_v2, new_v3);
}
