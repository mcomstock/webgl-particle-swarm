#version 300 es

precision highp int;
precision highp float;

layout (location = 0) out vec4 state_texture_0;
layout (location = 1) out vec4 state_texture_1;
layout (location = 2) out vec4 state_texture_2;
layout (location = 3) out vec4 state_texture_3;
layout (location = 4) out vec4 state_texture_4;
layout (location = 5) out vec4 state_texture_5;
layout (location = 6) out vec4 state_texture_6;

in vec2 cc;

uniform float V, Na_i, Na_ss, K_i, K_ss, Ca_i, Ca_ss, Ca_nsr, Ca_jsr, m, h, j,
    a, ifast, islow, aCaMK, iCaMKfast, iCaMKslow, d, ffast, fslow, fCafast,
    fCaslow, jCa, n, fCaMKfast, fCaCaMKfast, xrfast, xrslow, xs1, xs2, xK1, JrelNP, JrelCaMK,
    CaMKtrap;

void main() {
    state_texture_0 = vec4(V, Na_i, Na_ss, K_i);

    state_texture_1 = vec4(K_ss, Ca_i, Ca_ss, Ca_nsr);
    state_texture_2 = vec4(Ca_jsr, fslow, xrfast, xrslow);
    state_texture_3 = vec4(xs1, xK1, CaMKtrap, 0.0);

    state_texture_4 = vec4(
        uintBitsToFloat(packHalf2x16(vec2(m, h))),
        uintBitsToFloat(packHalf2x16(vec2(j, a))),
        uintBitsToFloat(packHalf2x16(vec2(ifast, islow))),
        uintBitsToFloat(packHalf2x16(vec2(aCaMK, iCaMKfast)))
    );

    state_texture_5 = vec4(
        uintBitsToFloat(packHalf2x16(vec2(iCaMKslow, d))),
        uintBitsToFloat(packHalf2x16(vec2(ffast, fCafast))),
        uintBitsToFloat(packHalf2x16(vec2(fCaslow, jCa))),
        uintBitsToFloat(packHalf2x16(vec2(n, fCaMKfast)))
    );

    state_texture_6 = vec4(
        uintBitsToFloat(packHalf2x16(vec2(fCaCaMKfast, xs2))),
        uintBitsToFloat(packHalf2x16(vec2(JrelNP, JrelCaMK))),
        0.0,
        0.0
    );
}
