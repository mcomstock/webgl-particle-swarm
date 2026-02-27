#version 300 es

precision highp float;
precision highp int;

uniform vec4 ics;

layout (location = 0) out vec4 state_texture;

void main() {
    state_texture = ics;
}
