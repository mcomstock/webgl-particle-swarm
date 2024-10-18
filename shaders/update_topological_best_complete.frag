#version 300 es

precision highp float;
precision highp int;

uniform sampler2D best_error_value_texture;

layout (location = 0) out uvec4 topological_best_idx;

in vec2 cc;

void main() {
    topological_best_idx = uvec4(texelFetch(best_error_value_texture, ivec2(0, 0), 0).yz, 0u, 0u);
}
