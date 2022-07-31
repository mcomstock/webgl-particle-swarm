#version 300 es

precision highp float;
precision highp int;

uniform sampler2D error_texture;

layout (location = 0) out vec4 expanded_error_texture;

in vec2 cc;

void main() {
    ivec2 error_size = textureSize(error_texture, 0);
    ivec2 idx = ivec2(floor(mod(cc, vec2(0.5, 0.5)) * 2.0 * vec2(error_size)));

    expanded_error_texture = texelFetch(error_texture, idx, 0);
}
