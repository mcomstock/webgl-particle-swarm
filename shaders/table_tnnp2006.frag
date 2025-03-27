#version 300 es

precision highp float;
precision highp int;

// EPI, MYO, or ENDO
#define ENDO

layout (location = 0) out vec4 table;

in vec2 cc;

uniform int width, height, npoints;
uniform float vmin, vmax, vekmin, vekmax, dt;

const float RR = 8314.3;
const float T = 310.0;
const float FF = 96486.7;
const float KmK = 1.0;
const float KmNa = 40.0;
const float gamma = 0.35;
const float KmCa = 1.38;
const float KmNai = 87.5;
const float Ko = 5.4;
const float Nao = 140.0;
const float Cao = 2.0;
const float ksat = 0.1;
const float alpha = 2.5;

void main() {
    ivec2 idx_2d = ivec2(floor(cc*vec2(width, height)));
    int idx = idx_2d[1] * width + idx_2d[0];

    int seq = idx / npoints;
    int seq_idx = idx % npoints;

    float V = (float(seq_idx)/float(npoints))*(vmax-vmin) + vmin;
    float VEK = (float(seq_idx)/float(npoints))*(vekmax-vekmin) + vekmin;

    switch (seq) {
    case 0:
        // m_inf
        table[0] = pow(1.0 / (1.0 + exp((-56.86 - V)/9.03)), 2.0);
        // tau_m exponential
        float alpha_m = 1.0 / (1.0 + exp((-60.0 - V)/5.0));
        float beta_m = (0.1 / (1.0 + exp((V + 35.0)/5.0))) + (0.1 / (1.0 + exp((V - 50.0)/200.0)));
        float tau_m = alpha_m * beta_m;
        table[1] = exp(-dt/tau_m);
        // h_inf
        table[2] = pow(1.0 / (1.0 + exp((V + 71.55)/7.43)), 2.0);
        // tau_h exponential
        float alpha_h, beta_h;
        if (V >= -40.0) {
            alpha_h = 0.0;
            beta_h = 0.77 / (0.13 * (1.0 + exp(-(V + 10.66)/11.1)));
        } else {
            alpha_h = 0.057 * exp(-(V + 80.0)/6.8);
            beta_h = 2.7 * exp(0.079 * V) + 3.1e5 * exp(0.3485 * V);
        }
        float tau_h = 1.0 / (alpha_h + beta_h);
        table[3] = exp(-dt/tau_h);
        break;
    case 1:
        // j_inf
        table[0] = pow(1.0 / (1.0 + exp((V + 71.55)/7.43)), 2.0);
        // tau_j exponential
        float alpha_j, beta_j;
        if (V >= -40.0) {
            alpha_j = 0.0;
            beta_j = (0.6 * exp(0.057 * V)) / (1.0 + exp(-0.1*(V+32.0)));
        } else {
            alpha_j = ((-2.5428e4*exp(0.2444*V)-6.948e-6*exp(-0.04391*V))*(V+37.78)) / (1.0 + exp(0.311*(V+79.23)));
            beta_j = (0.02424*exp(-0.01052*V)) / (1.0 + exp(-0.1378*(V+40.14)));
        }
        float tau_j = 1.0 / (alpha_j + beta_j);
        table[1] = exp(-dt/tau_j);
        // xs_inf
        table[2] = 1.0 / (1.0 + exp((-5.0 - V)/14.0));
        // tau_xs exponential
        float alpha_xs = 1400.0 / sqrt(1.0 + exp((5.0-V)/6.0));
        float beta_xs = 1.0 / (1.0 + exp((V-35.0)/15.0));
        float tau_xs = alpha_xs * beta_xs + 80.0;
        table[3] = exp(-dt/tau_xs);
        break;
    case 2:
        // d_inf
        table[0] = 1.0 / (1.0 + exp((-8.0 - V)/7.5));
        // tau_d exponential
        float alpha_d = 1.4 / (1.0 + exp((-35.0-V)/13.0)) + 0.25;
        float beta_d = 1.4 / (1.0 + exp((V+5.0)/5.0));
        float gamma_d = 1.0 / (1.0 + exp((50.0-V)/20.0));
        float tau_d = alpha_d * beta_d + gamma_d;
        table[1] = exp(-dt/tau_d);
        // f_inf
        table[2] = 1.0 / (1.0 + exp((V+20.0)/7.0));
        // tau_f exponential
        float alpha_f = 1102.5 * exp(-((V+27.0)/15.0)*((V+27.0)/15.0));
        float beta_f = 200.0 / (1.0 + exp((13.0-V)/10.0));
        float gamma_f = 180.0 / (1.0 + exp((V+30.0)/10.0)) + 20.0;
        float tau_f = alpha_f + beta_f + gamma_f;
        table[3] = exp(-dt/tau_f);
        break;
    case 3:
        // f2_inf
        table[0] = 0.67 / (1.0 + exp((V+35.0)/7.0)) + 0.33;
        // tau_f2 exponential
        // The equation here differs from the paper (original commented out)
        // float alpha_f2 = 600.0 * exp(-((V+25.0)*(V+25.0))/170.0)
        float alpha_f2 = 600.0 * exp(-((V+25.0)*(V+25.0))/49.0);
        float beta_f2 = 31.0 / (1.0 + exp((25.0 - V)/10.0));
        float gamma_f2 = 16.0 / (1.0 + exp((V+30.0)/10.0));
        float tau_f2 = alpha_f2 + beta_f2 + gamma_f2;
        table[1] = exp(-dt/tau_f2);
        // r_inf
        table[2] = 1.0 / (1.0 + exp((20.0 - V)/6.0));
        // tau_r exponential
        float tau_r = 9.5 * exp(-(V+40.0)*(V+40.0)/1800.0) + 0.8;
        table[3] = exp(-dt/tau_r);
        break;
    case 4:
        // s_inf
#if defined ENDO
        table[0] = 1.0 / (1.0 + exp((V+28.0)/5.0));
#else
        table[0] = 1.0 / (1.0 + exp((V+20.0)/5.0));
#endif
        // tau_s exponential
#if defined ENDO
        float tau_s = 1000.0 * exp(-(V+67.0)*(V+67.0)/1000.0) + 8.0;
#else
        float tau_s = 85.0 * exp(-(V+45.0)*(V+45.0)/320.0) + 5.0 / (1.0 + exp((V-20.0)/5.0)) + 3.0;
#endif
        table[1] = exp(-dt/tau_s);
        // xr1_inf
        table[2] = 1.0 / (1.0 + exp((-26.0 - V)/7.0));
        // tau_xr1 exponential
        float alpha_xr1 = 450.0 / (1.0 + exp((-45.0 - V)/10.0));
        float beta_xr1 = 6.0 / (1.0 + exp((V+30.0)/11.5));
        float tau_xr1 = alpha_xr1 * beta_xr1;
        table[3] = exp(-dt/tau_xr1);
        break;
    case 5:
        // xr2_inf
        table[0] = 1.0 / (1.0 + exp((V+88.0)/24.0));
        // tau_xr2 exponential
        float alpha_xr2 = 3.0 / (1.0 + exp((-60.0 - V)/20.0));
        float beta_xr2 = 1.12 / (1.0 + exp((V-60.0)/20.0));
        float tau_xr2 = alpha_xr2 * beta_xr2;
        table[1] = exp(-dt/tau_xr2);
        // ICaL CaSS coefficient
        // L'Hopital
        if (abs(V-15.0) < 0.01) {
            table[2] = 0.5 * FF;
        } else {
            table[2] = (((FF*FF)/(RR*T)) * (V-15.0) * exp(2.0 * (V-15.0) * (FF/(RR*T)))) / (exp(2.0 * (V-15.0) * (FF/(RR*T))) - 1.0);
        }
        // ICaL Cao term
        // L'Hopital
        if (abs(V-15.0) < 0.01) {
            table[3] = 2.0 * FF * Cao;
        } else {
            table[3] = 4.0 * (((FF*FF)/(RR*T)) * (V-15.0) * Cao) / (exp(2.0 * (V-15.0) * (FF/(RR*T))) - 1.0);
        }
        break;
    case 6:
        // IpK coefficient
        table[0] = 1.0 / (1.0 + exp((25.0-V)/5.98));
        // INaCa Nai coeffecient
        float INaCa_denom = (KmNai*KmNai*KmNai + Nao*Nao*Nao) * (KmCa+Cao) * (1.0 + ksat*exp((gamma-1.0)*V*FF/(RR*T)));
        table[1] = Cao * exp(gamma*V*FF/(RR*T)) / INaCa_denom;
        // INaCa Cai coeffecient
        table[2] = alpha * Nao * Nao * Nao * exp((gamma-1.0)*V*FF/(RR*T)) / INaCa_denom;
        // INaK coefficient
        table[3] = Ko / ((Ko+KmK) * (1.0 + 0.1245*exp(-0.1*V*FF/(RR*T)) + 0.0353 * exp(-V*FF/(RR*T))));
        break;
    case 7:
        // xK1_inf (VEK)
        float alpha_K1 = 0.1 / (1.0 + exp(0.06*((VEK)-200.0)));
        float beta_K1 = (3.0 * exp(0.0002*((VEK)+100.0)) + exp(0.1*((VEK)-10.0))) / (1.0 + exp(-0.5*(VEK)));
        table[0] = alpha_K1 / (alpha_K1 + beta_K1);
        break;
    default:
        break;
    }
}
