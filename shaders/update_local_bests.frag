#version 300 es

precision highp float;
precision highp int;

uniform sampler2D local_bests_texture, local_bests_error_texture, cur_vals_texture, cur_error_texture;

layout (location = 0) out vec4 new_local_best;
layout (location = 1) out vec4 new_local_best_error;

in vec2 cc;

void main() {
    ivec2 dims = textureSize(cur_error_texture, 0);
    ivec2 idx = ivec2(floor(cc * vec2(dims)));

    vec4 cur_error = texelFetch(cur_error_texture, idx, 0);
    vec4 local_best_error = texelFetch(local_bests_error_texture, idx, 0);

    vec4 local_bests = texelFetch(local_bests_texture, idx, 0);
    vec4 cur_vals = texelFetch(cur_vals_texture, idx, 0);

    if (cur_error.r < local_best_error.r) {
        new_local_best = cur_vals;
        new_local_best_error = cur_error;
    } else {
        new_local_best = local_bests;
        new_local_best_error = local_best_error;
    }
}
