#version 300 es

precision highp float;
precision highp int;

uniform sampler2D best_error_value_texture;

layout (location = 0) out vec4 topological_best_idx;

in vec2 cc;

void main() {
    topological_best_idx = vec4(texelFetch(best_error_value_texture, ivec2(0, 0), 0).yz, 0.0, 0.0);
}
