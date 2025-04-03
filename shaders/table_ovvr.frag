#version 300 es

precision highp float;
precision highp int;

layout (location = 0) out vec4 table;

in vec2 cc;

uniform int width, height, npoints;
uniform float vmin, vmax, dt;

const float RR = 8314.0;
const float TT = 310.0;
const float FF = 96485.0;
const float RTF = (RR*TT)/FF;
const float FRT = FF/(RR*TT);
const float FFRT = (FF*FF)/(RR*TT);
const float zCa = 2.0;
const float zNa = 1.0;
const float zK = 1.0;
const float qNa = 0.5224;
const float qCa = 0.1670;
const float Delta = -0.1550;
const float KNaio = 9.073;
const float KNaoo = 27.78;
const float K_o = 5.4;

void main() {
    ivec2 idx_2d = ivec2(floor(cc*vec2(width, height)));
    int idx = idx_2d[1] * width + idx_2d[0];

    int seq = idx / npoints;
    int seq_idx = idx % npoints;

    float V = (float(seq_idx)/float(npoints))*(vmax-vmin) + vmin;

    switch (seq) {
    case 0:
        // minf
        table[0] = 1.0 / (1.0 + exp(-(V + 39.57)/9.871));
        // taum_exp
        float taum = 1.0 / (6.765 * exp((V+11.64)/34.77) + 8.552 * exp(-(V+77.42)/5.955));
        table[1] = exp(-dt/taum);
        // hinf
        table[2] = 1.0 / (1.0 + exp((V + 82.9)/6.086));
        // tauhfast_exp
        float tauhfast = 1.0 / (1.432e-5 * exp(-(V+1.196)/6.285) + 6.149 * exp((V+0.5096)/20.27));
        table[3] = exp(-dt/tauhfast);
        break;
    case 1:
        // tauhslow_exp
        float tauhslow = 1.0 / (0.009764 * exp(-(V+17.95)/28.05) + 0.3343 * exp((V+5.730)/56.66));
        table[0] = exp(-dt/tauhslow);
        // tauj_exp
        float tauj = 2.038 + 1.0 / (0.02136 * exp(-(V+100.6)/8.281) + 0.3052 * exp((V+0.9941)/38.45));
        table[1] = exp(-dt/tauj);
        // hCaMKinf
        table[2] = 1.0 / (1.0 + exp((V+89.1)/6.086));
        // tauhCaMKslow_exp
        float tauhCaMKslow = 3.0 * tauhslow;
        table[3] = exp(-dt/tauhCaMKslow);
        break;
    case 2:
        // mLinf
        table[0] = 1.0 / (1.0 + exp(-(V+42.85)/5.264));
        // hLinf
        table[1] = 1.0 / (1.0 + exp((V+87.61)/7.488));
        // hLCaMKinf
        table[2] = 1.0 / (1.0 + exp((V+93.81)/7.488));
        // ainf
        table[3] = 1.0 / (1.0 + exp(-(V-14.34)/14.82));
        break;
    case 3:
        // taua_exp
        float taua = 1.0515 / ((1.0 / (1.2089 * (1.0 + exp(-(V-18.41)/29.38)))) + (3.5 / (1.0 + exp((V+100.0)/29.38))));
        table[0] = exp(-dt/taua);
        // iinf
        table[1] = 1.0 / (1.0 + exp((V+43.94)/5.711));
        // tauifast_exp
        float tauifast = 4.562 + 1.0 / (0.3933 * exp(-(V+100.0)/100.0) + 0.08004 * exp((V+50.0)/16.59));
        table[2] = exp(-dt/tauifast);
        // tauislow_exp
        float tauislow = 23.62 + 1.0 / (0.001416 * exp(-(V+96.52)/59.05) + 1.7808e-8 * exp((V+114.1)/8.079));
        table[3] = exp(-dt/tauislow);
        break;
    case 4:
        // Aifast
        table[0] = 1.0 / (1.0 + exp((V-213.6)/151.2));
        // aCaMKinf
        table[1] = 1.0 / (1.0 + exp(-(V-24.34)/14.82));
        // tauiCaMKfast_exp
        tauifast = 4.562 + 1.0 / (0.3933 * exp(-(V+100.0)/100.0) + 0.08004 * exp((V+50.0)/16.59));
        float deltaCaMKdevelop = 1.354 + 1e-4 / (exp((V-167.4)/15.89) + exp(-(V-12.23)/0.2154));
        float deltaCaMKrecover = 1.0 - 0.5 / (1.0 + exp((V+70.0)/20.0));
        float tauiCaMKfast = tauifast * deltaCaMKdevelop * deltaCaMKrecover;
        table[2] = exp(-dt/tauiCaMKfast);
        // tauiCaMKslow_exp
        tauislow = 23.62 + 1.0 / (0.001416 * exp(-(V+96.52)/59.05) + 1.7808e-8 * exp((V+114.1)/8.079));
        float tauiCaMKslow = tauislow * deltaCaMKdevelop * deltaCaMKrecover;
        table[3] = exp(-dt/tauiCaMKslow);
        break;
    case 5:
        // dinf
        table[0] = 1.0 / (1.0 + exp(-(V+3.940)/4.230));
        // taud_exp
        float taud = 0.6 + 1.0 / (exp(-0.05 * (V+6.0)) + exp(0.09 * (V+14.0)));
        table[1] = exp(-dt/taud);
        // finf
        table[2] = 1.0 / (1.0 + exp((V+19.58)/3.696));
        // tauffast_exp
        // Note: This equation seems to be implemented correctly, but does not match the plot in Figure
        // 1C of the paper.
        float tauffast = 7.0 + 1.0 / (0.0045 * exp(-(V+20.0)/10.0) + 0.0045 * exp((V+20.0)/10.0));
        table[3] = exp(-dt/tauffast);
        break;
    case 6:
        // taufslow_exp
        float taufslow = 1000.0 + 1.0 / (0.000035 * exp(-(V+5.0)/4.0) + 0.000035 * exp((V+5.0)/6.0));
        table[0] = exp(-dt/taufslow);
        // taufCafast_exp
        float taufCafast = 7.0 + 1.0 / (0.04 * exp(-(V-4.0)/7.0) + 0.04 * exp((V-4.0)/7.0));
        table[1] = exp(-dt/taufCafast);
        // taufCaslow_exp
        float taufCaslow = 100.0 + 1.0 / (0.00012 * exp(-V/3.0) + 0.00012 * exp(V/7.0));
        table[2] = exp(-dt/taufCaslow);
        // AfCafast
        table[3] = 0.3 + 0.6 / (1.0 + exp((V-10.0)/10.0));
        break;
    case 7:
        // taufCaMKfast_exp
        tauffast = 7.0 + 1.0 / (0.0045 * exp(-(V+20.0)/10.0) + 0.0045 * exp((V+20.0)/10.0));
        float taufCaMKfast = 2.5 * tauffast;
        table[0] = exp(-dt/taufCaMKfast);
        // taufCaCaMKfast_exp
        taufCafast = 7.0 + 1.0 / (0.04 * exp(-(V-4.0)/7.0) + 0.04 * exp((V-4.0)/7.0));
        float taufCaCaMKfast = 2.5 * taufCafast;
        table[1] = exp(-dt/taufCaCaMKfast);
        // exp_zcavfrt
        table[2] = exp(zCa*V*FRT);
        // exp_znavfrt
        table[3] = exp(zNa*V*FRT);
        break;
    case 8:
        // exp_zkvfrt
        table[0] = exp(zK*V*FRT);
        // xrinf
        table[1] = 1.0 / (1.0 + exp(-(V+8.337)/6.789));
        // tauxrfast_exp
        float tauxrfast = 12.98 + 1.0 / (0.3652 * exp((V-31.66)/3.869) + 4.123e-5 * exp(-(V-47.78)/20.38));
        table[2] = exp(-dt/tauxrfast);
        // tauxrslow_exp
        float tauxrslow = 1.865 + 1.0 / (0.06629 * exp((V-34.70)/7.355) + 1.128e-5 * exp(-(V-29.74)/25.94));
        table[3] = exp(-dt/tauxrslow);
        break;
    case 9:
        // Axrfast
        table[0] = 1.0 / (1.0 + exp((V+54.81)/38.21));
        // RKr
        table[1] = 1.0 / ((1.0 + exp((V+55.0)/75.0)) * (1.0 + exp((V-10.0)/30.0)));
        // xs1inf
        table[2] = 1.0 / (1.0 + exp(-(V+11.60)/8.932));
        // tauxs1_exp
        float tauxs1 = 817.3 + 1.0 / (2.326e-4 * exp((V+48.28)/17.80) + 0.001292 * exp(-(V+210.0)/230.0));
        table[3] = exp(-dt/tauxs1);
        break;
    case 10:
        // tauxs2_exp
        float tauxs2 = 1.0 / (0.01 * exp((V-50.0)/20.0) + 0.0193 * exp(-(V+66.54)/31.0));
        table[0] = exp(-dt/tauxs2);
        // xK1inf
        table[1] = 1.0 / (1.0 + exp(-(V+2.5538*K_o+144.59)/(1.5692*K_o+3.8115)));
        // tauxK1_exp
        float tauxK1 = 122.2 / (exp(-(V+127.2)/20.36) + exp((V+236.8)/69.33));
        table[2] = exp(-dt/tauxK1);
        // RK1
        table[3] = 1.0 / (1.0 + exp((V+105.8-2.6*K_o)/9.493));
        break;
    case 11:
        // invhCa
        table[0] = 1.0 / exp(qCa*V*FRT);
        // hNa
        table[1] = exp(qNa*V*FRT);
        // invKNai
        table[2] = 1.0 / (KNaio * exp((Delta*V*FRT)/3.0));
        // invKNao
        table[3] = 1.0 / (KNaoo * exp(((1.0-Delta)*V*FRT)/3.0));
        break;
    case 12:
        // INab_coeff
        if (abs(V) < 0.01) {
            // l'hopital
            table[0] = zNa * FF;
        } else {
            table[0] = zNa * zNa * V * FFRT / (exp(zNa*V*FRT) - 1.0);
        }
        // ICab_coeff
        if (V == 0.0) {
            // l'hopital
            table[1] = zCa * FF;
        } else {
            table[1] = zCa * zCa * V * FFRT / (exp(zCa*V*FRT) - 1.0);
        }
        // xKb
        table[2] = 1.0 / (1.0 + exp(-(V-14.48)/18.34));
        // taujCaMK_exp
        tauj = 2.038 + 1.0 / (0.02136 * exp(-(V+100.6)/8.281) + 0.3052 * exp((V+0.9941)/38.45));
        table[3] = exp(-dt/(1.46*tauj));
        break;
    default:
        break;
    }
}
