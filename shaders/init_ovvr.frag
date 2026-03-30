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
layout (location = 7) out vec4 state_texture_7;

in vec2 cc;

uniform float V, Na_i, Na_ss, K_i, K_ss, Ca_i, Ca_ss, Ca_nsr, Ca_jsr, m, hfast, hslow, j, hCaMKslow,
    jCaMK, mL, hL, hLCaMK, a, ifast, islow, aCaMK, iCaMKfast, iCaMKslow, d, ffast, fslow, fCafast,
    fCaslow, jCa, n, fCaMKfast, fCaCaMKfast, xrfast, xrslow, xs1, xs2, xK1, JrelNP, JrelCaMK,
    CaMKtrap;

void main() {
    state_texture_0 = vec4(V, Na_i, Na_ss, K_i);
    state_texture_1 = vec4(K_ss, Ca_i, Ca_ss, Ca_nsr);
    state_texture_2 = vec4(Ca_jsr, hLCaMK, fslow, xrfast);
    state_texture_3 = vec4(xrslow, xs1, xK1, CaMKtrap);

    state_texture_4 = vec4(
        uintBitsToFloat(packHalf2x16(vec2(fCafast, j))),
        jCaMK,
        uintBitsToFloat(packHalf2x16(vec2(mL, JrelCaMK))),
        uintBitsToFloat(packHalf2x16(vec2(a, ifast)))
    );

    state_texture_5 = vec4(
        uintBitsToFloat(packHalf2x16(vec2(islow, aCaMK))),
        uintBitsToFloat(packHalf2x16(vec2(iCaMKfast, iCaMKslow))),
        uintBitsToFloat(packHalf2x16(vec2(d, ffast))),
        uintBitsToFloat(packHalf2x16(vec2(jCa, n)))
    );

    state_texture_6 = vec4(
        uintBitsToFloat(packHalf2x16(vec2(fCaMKfast, fCaCaMKfast))),
        uintBitsToFloat(packHalf2x16(vec2(xs2, JrelNP))),
        hslow,
        hCaMKslow
    );

    state_texture_7 = vec4(
        hfast,
        hL,
        fCaslow,
        m
    );
}
