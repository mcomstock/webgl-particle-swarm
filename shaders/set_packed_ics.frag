#version 300 es

precision highp float;
precision highp int;

uniform vec4 ics_0, ics_1;

layout (location = 0) out vec4 state_texture;

void main() {
    state_texture = vec4(
        uintBitsToFloat(packHalf2x16(ics_0.xy)),
        uintBitsToFloat(packHalf2x16(ics_0.zw)),
        uintBitsToFloat(packHalf2x16(ics_1.xy)),
        uintBitsToFloat(packHalf2x16(ics_1.zw))
    );
}
