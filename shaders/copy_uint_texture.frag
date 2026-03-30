#version 300 es

precision highp int;
precision highp float;
precision highp usampler2D;

uniform usampler2D original;

in vec2 cc;

layout (location = 0) out uvec4 copy;

void main() {
    ivec2 tex_size = textureSize(original, 0);
    ivec2 idx = ivec2(floor(cc * vec2(tex_size)));
    copy = texelFetch(original, idx, 0);
}
