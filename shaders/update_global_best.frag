#version 300 es

precision highp float;
precision highp int;

uniform sampler2D positions_texture, reduced_error_2_texture, global_best_texture, best_error_value_texture;

in vec2 cc;

layout (location = 0) out vec4 global_best_out;
layout (location = 1) out vec4 best_error_value_out;

void main() {
    vec4 re2 = texelFetch(reduced_error_2_texture, ivec2(0, 0), 0);
    vec4 best_error = texture(best_error_value_texture, cc);

    if (re2.x < best_error.x) {
        ivec2 dims = textureSize(positions_texture, 0);
        ivec2 global_best_idx = ivec2(re2.yz);
        global_best_idx += ivec2(round(cc)) * (dims/2);
        global_best_out = texelFetch(positions_texture, global_best_idx, 0);
        best_error_value_out = vec4(re2.xyz, 0.0);
    } else {
        global_best_out = texture(global_best_texture, cc);
        best_error_value_out = best_error;
    }
}
