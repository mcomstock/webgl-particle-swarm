#version 300 es

precision highp float;
precision highp int;

// EPI, MYO, or ENDO
#define ENDO

uniform sampler2D in_particles_1, data_texture;

layout (location = 0) out vec4 error_texture;

in vec2 cc;

uniform float dt, period;
uniform int num_beats, pre_beats, data_type;
uniform float align_thresh;
uniform float sample_interval, apd_thresh, weight;
uniform float stim_dur, stim_mag, stim_offset_1, stim_offset_2, stim_t_scale;
uniform bool stim_biphasic;

float RR = 8314.3;
float T = 310.0;
float FF = 96486.7;
float Cm = 1.0;
float SS = 0.2;
float rho = 162.0;
float Vc = 0.016404;
float Vsr = 0.001094;
float Vss = 0.00005468;
float Ko = 5.4;
float Nao = 140.0;
float Cao = 2.0;
// float GNa = 14.838;
// float GK1 = 5.405;
#if defined ENDO
// float Gto = 0.073;
#else
// float Gto = 0.294;
#endif
// float GKr = 0.153;
#if defined MYO
// float GKs = 0.098;
#else
// float GKs = 0.392;
#endif
float pKNa = 0.03;
// float GCaL = 3.98e-5;
// float kNaCa = 1000.0;
float gamma = 0.35;
float KmCa = 1.38;
float KmNai = 87.5;
float ksat = 0.1;
float alpha = 2.5;
// float PNaK = 2.724;
float KmK = 1.0;
float KmNa = 40.0;
// float GpK = 0.0146;
// float GpCa = 0.1238;
float KpCa = 0.0005;
// float GbNa = 0.00029;
// float GbCa = 0.000592;
float Vmaxup = 0.006375;
float Kup = 0.00025;
float Vrel = 0.102;
float k1p = 0.15;
float k2p = 0.045;
float k3 = 0.06;
float k4 = 0.005;
float EC = 1.5;
float maxsr = 2.5;
float minsr = 1.0;
float Vleak = 0.00036;
float Vxfer = 0.0038;
float Bufc = 0.2;
float Kbufc = 0.001;
float Bufsr = 10.0;
float Kbufsr = 0.3;
float Bufss = 0.4;
float Kbufss = 0.00025;
float capacitance = 0.185;

float biphasic_stim_f(const float t) {
    float a = (t/stim_t_scale - stim_offset_2);

    return -stim_mag * (t / stim_t_scale - stim_offset_1) / (1.0 + a*a*a*a);
}

float square_stim_f(const float t) {
    return stim_mag;
}

void main() {
    int num_period = int(ceil(period/dt));
    int total_beats = pre_beats + num_beats;
    float endtime = ceil(float(total_beats)*period);
    float pre_pace_endtime = ceil(float(pre_beats)*period);
    int pre_pace_steps = int(ceil(pre_pace_endtime/dt));
    int num_steps = int(ceil(endtime/dt));

    ivec2 tex_size = textureSize(in_particles_1, 0);
    ivec2 idx = ivec2(floor(cc * 0.5 * vec2(tex_size)));

    int num_data_points = textureSize(data_texture, 0).x;

    vec4 particles_1 = texelFetch(in_particles_1, idx, 0);
    vec4 particles_2 = texelFetch(in_particles_1, idx + ivec2(tex_size.x/2, 0), 0);
    vec4 particles_3 = texelFetch(in_particles_1, idx + ivec2(0, tex_size.y/2), 0);
    vec4 particles_4 = texelFetch(in_particles_1, idx + ivec2(tex_size.x/2, tex_size.y/2), 0);

    float GNa = particles_1[0];
    float GK1 = particles_1[1];
    float Gto = particles_1[2];
    float GKr = particles_1[3];

    float GKs = particles_2[0];
    float GCaL = particles_2[1];
    float GpK = particles_2[2];
    float GpCa = particles_2[3];

    float GbNa = particles_3[0];
    float GbCa = particles_3[1];
    float PNaK = particles_3[2];
    float kNaCa = particles_3[3];

    float V, Cai, CaSR, CaSS, Nai, Ki, Rhat;
    float m, h, j, r, s, xr1, xr2, xs, d, f, f2, fcass;

    float m_inf, h_inf, j_inf, xs_inf, d_inf, f_inf, f2_inf, fcass_inf, r_inf, s_inf, xr1_inf, xr2_inf, xK1_inf;
    float alpha_m, alpha_h, alpha_j, alpha_xs, alpha_d, alpha_f, alpha_f2, alpha_xr1, alpha_xr2, alpha_K1;
    float beta_m, beta_h, beta_j, beta_xs, beta_d, beta_f, beta_f2, beta_xr1, beta_xr2, beta_K1;
    float tau_m, tau_h, tau_j, tau_xs, tau_d, tau_f, tau_f2, tau_fcass, tau_r, tau_s, tau_xr1, tau_xr2;
    float gamma_d, gamma_f, gamma_f2;

    float Iion, INa, IK1, Ito, IKr, IKs, ICaL, INaCa, INaK, IpCa, IpK, IbCa, IbNa, Istim;
    float Ileak, Iup, Irel, Ixfer;

    float ENa, EK, ECa, EKs;

    float Caibufc, Casrbufsr, Cassbufss, dCaitotal, dCaSRtotal, dCaSStotal, bi, ci, bsr, csr, bss, css, k1, k2, kcasr, O;

    float time;

#if defined EPI
    V = -85.46;
    Rhat = 0.9891;
    Nai = 9.293;
    Ki = 136.2;
    Cai = 0.0001156;
    CaSS = 0.002331;
    CaSR = 3.432;
    m = 0.001633;
    h = 0.7512;
    j = 0.7508;
    xs = 0.003214;
    d = 3.72e-5;
    f = 0.9767;
    f2 = 0.9995;
    fcass = 1.0;
#elif defined MYO
    V = -84.53;
    Rhat = 0.9874;
    Nai = 9.322;
    Ki = 136.0;
    Cai = 0.0001156;
    CaSS = 0.002331;
    CaSR = 4.130;
    m = 0.001694;
    h = 0.7466;
    j = 0.7457;
    xs = 0.003343;
    d = 3.345-5;
    f = 0.9595;
    f2 = 0.9995;
    fcass = 1.0;
#else
    V = -84.70;
    Rhat = 0.9891;
    Nai = 9.413;
    Ki = 136.1;
    Cai = 0.0001021;
    CaSS = 0.002111;
    CaSR = 3.385;
    m = 0.001634;
    h = 0.7512;
    j = 0.7508;
    xs = 0.003213;
    d = 3.27e-5;
    f = 0.9771;
    f2 = 0.9995;
    fcass = 1.0;
#endif
    r = 0.0;
    s = 0.0;
    xr1 = 0.0;
    xr2 = 0.0;

    float compare_stride = round(sample_interval / dt);

    float error = (data_type == 1) ? 0.0 : 10000000000.0;

    int data_index = 0;

    int start_comp = 0;
    int compared_points = 0;
    bool first_align_upstroke = false;
    float saved_value = -1.0;

    float APD_start, APD_end;

    bool activated = false;

    for (int step_count = 1; step_count <= num_steps; ++step_count) {
        float stim = 0.0;
        float stim_t = mod(float(step_count)*dt, period);
        if (stim_t < stim_dur) {
            if (stim_biphasic) {
                stim = biphasic_stim_f(stim_t);
            } else {
                stim = square_stim_f(stim_t);
            }
        }

        Istim = stim;

        m_inf = 1.0 / (1.0 + exp((-56.86 - V)/9.03));
        m_inf = m_inf * m_inf;
        alpha_m = 1.0 / (1.0 + exp((-60.0 - V)/5.0));
        beta_m = (0.1 / (1.0 + exp((V + 35.0)/5.0))) + (0.1 / (1.0 + exp((V - 50.0)/200.0)));
        tau_m = alpha_m * beta_m;
        m = m_inf - (m_inf - m) * exp(-dt/tau_m);

        h_inf = 1.0 / (1.0 + exp((V + 71.55)/7.43));
        h_inf = h_inf * h_inf;
        if (V >= -40.0) {
            alpha_h = 0.0;
            beta_h = 0.77 / (0.13 * (1.0 + exp(-(V + 10.66)/11.1)));
        } else {
            alpha_h = 0.057 * exp(-(V + 80.0)/6.8);
            beta_h = 2.7 * exp(0.079 * V) + 3.1e5 * exp(0.3485 * V);
        }
        tau_h = 1.0 / (alpha_h + beta_h);
        h = h_inf - (h_inf - h) * exp(-dt/tau_h);

        j_inf = 1.0 / (1.0 + exp((V + 71.55)/7.43));
        j_inf = j_inf * j_inf;
        if (V >= -40.0) {
            alpha_j = 0.0;
            beta_j = (0.6 * exp(0.057 * V)) / (1.0 + exp(-0.1*(V+32.0)));
        } else {
            alpha_j = ((-2.5428e4*exp(0.2444*V)-6.948e-6*exp(-0.04391*V))*(V+37.78)) / (1.0 + exp(0.311*(V+79.23)));
            beta_j = (0.02424*exp(-0.01052*V)) / (1.0 + exp(-0.1378*(V+40.14)));
        }
        tau_j = 1.0 / (alpha_j + beta_j);
        j = j_inf - (j_inf - j) * exp(-dt/tau_j);

        xs_inf = 1.0 / (1.0 + exp((-5.0 - V)/14.0));
        alpha_xs = 1400.0 / sqrt(1.0 + exp((5.0-V)/6.0));
        beta_xs = 1.0 / (1.0 + exp((V-35.0)/15.0));
        tau_xs = alpha_xs * beta_xs + 80.0;
        xs = xs_inf - (xs_inf - xs) * exp(-dt/tau_xs);

        d_inf = 1.0 / (1.0 + exp((-8.0 - V)/7.5));
        alpha_d = 1.4 / (1.0 + exp((-35.0-V)/13.0)) + 0.25;
        beta_d = 1.4 / (1.0 + exp((V+5.0)/5.0));
        gamma_d = 1.0 / (1.0 + exp((50.0-V)/20.0));
        tau_d = alpha_d * beta_d + gamma_d;
        d = d_inf - (d_inf - d) * exp(-dt/tau_d);

        f_inf = 1.0 / (1.0 + exp((V+20.0)/7.0));
        alpha_f = 1102.5 * exp(-((V+27.0)/15.0)*((V+27.0)/15.0));
        beta_f = 200.0 / (1.0 + exp((13.0-V)/10.0));
        gamma_f = 180.0 / (1.0 + exp((V+30.0)/10.0)) + 20.0;
        tau_f = alpha_f + beta_f + gamma_f;
        f = f_inf - (f_inf - f) * exp(-dt/tau_f);

        f2_inf = 0.67 / (1.0 + exp((V+35.0)/7.0)) + 0.33;
        // The equation here differs from the paper (original commented out)
        // alpha_f2 = 600.0 * exp(-((V+25.0)*(V+25.0))/170.0)
        alpha_f2 = 600.0 * exp(-((V+25.0)*(V+25.0))/49.0);
        beta_f2 = 31.0 / (1.0 + exp((25.0 - V)/10.0));
        gamma_f2 = 16.0 / (1.0 + exp((V+30.0)/10.0));
        tau_f2 = alpha_f2 + beta_f2 + gamma_f2;
        f2 = f2_inf - (f2_inf - f2) * exp(-dt/tau_f2);

        fcass_inf = 0.6 / (1.0 + (CaSS/0.05)*(CaSS/0.05)) + 0.4;
        tau_fcass = 80.0 / (1.0 + (CaSS/0.05)*(CaSS/0.05)) + 2.0;
        fcass = fcass_inf - (fcass_inf - fcass) * exp(-dt/tau_fcass);

        r_inf = 1.0 / (1.0 + exp((20.0 - V)/6.0));
        tau_r = 9.5 * exp(-(V+40.0)*(V+40.0)/1800.0) + 0.8;
        r = r_inf - (r_inf - r) * exp(-dt/tau_r);

#if defined ENDO
        s_inf = 1.0 / (1.0 + exp((V+28.0)/5.0));
        tau_s = 1000.0 * exp(-(V+67.0)*(V+67.0)/1000.0) + 8.0;
#else
        s_inf = 1.0 / (1.0 + exp((V+20.0)/5.0));
        tau_s = 85.0 * exp(-(V+45.0)*(V+45.0)/320.0) + 5.0 / (1.0 + exp((V-20.0)/5.0)) + 3.0;
#endif
        s = s_inf - (s_inf - s) * exp(-dt/tau_s);

        xr1_inf = 1.0 / (1.0 + exp((-26.0 - V)/7.0));
        alpha_xr1 = 450.0 / (1.0 + exp((-45.0 - V)/10.0));
        beta_xr1 = 6.0 / (1.0 + exp((V+30.0)/11.5));
        tau_xr1 = alpha_xr1 * beta_xr1;
        xr1 = xr1_inf - (xr1_inf - xr1) * exp(-dt/tau_xr1);

        xr2_inf = 1.0 / (1.0 + exp((V+88.0)/24.0));
        alpha_xr2 = 3.0 / (1.0 + exp((-60.0 - V)/20.0));
        beta_xr2 = 1.12 / (1.0 + exp((V-60.0)/20.0));
        tau_xr2 = alpha_xr2 * beta_xr2;
        xr2 = xr2_inf - (xr2_inf - xr2) * exp(-dt/tau_xr2);

        ENa = (RR*T/FF) * log(Nao/Nai);
        EK = (RR*T/FF) * log(Ko/Ki);
        ECa = (RR*T/(2.0*FF)) * log(Cao/Cai);
        EKs = (RR*T/FF) * log((Ko+pKNa*Nao)/(Ki+pKNa*Nai));

        ICaL = GCaL * d * f * f2 * fcass * 4.0 * (((V-15.0)*FF*FF)/(RR*T)) * ((0.25 * CaSS * exp(2.0 * (V-15.0) * FF/(RR*T)) - Cao) / (exp(2.0 * (V-15.0) * FF/(RR*T)) - 1.0));

        IKs = GKs * xs * xs * (V - EKs);
        INa = GNa * m * m * m * h * j * (V - ENa);
        Ito = Gto * r * s * (V - EK);
        IKr = GKr * sqrt(Ko/5.4) * xr1 * xr2 * (V - EK);

        alpha_K1 = 0.1 / (1.0 + exp(0.06*((V-EK)-200.0)));
        beta_K1 = (3.0 * exp(0.0002*((V-EK)+100.0)) + exp(0.1*((V-EK)-10.0))) / (1.0 + exp(-0.5*(V-EK)));
        xK1_inf = alpha_K1 / (alpha_K1 + beta_K1);
        // Note that the sqrt term is 1 with the default parameterization
        IK1 = GK1 * sqrt(Ko/5.4) * xK1_inf * (V-EK);

        INaCa = kNaCa * (exp(gamma*V*FF/(RR*T))*Nai*Nai*Nai*Cao - exp((gamma-1.0)*V*FF/(RR*T))*Nao*Nao*Nao*Cai*alpha) / ((KmNai*KmNai*KmNai + Nao*Nao*Nao) * (KmCa+Cao) * (1.0 + ksat*exp((gamma-1.0)*V*FF/(RR*T))));
        INaK = PNaK*Ko*Nai / ((Ko+KmK) * (Nai+KmNa) * (1.0 + 0.1245*exp(-0.1*V*FF/(RR*T)) + 0.0353 * exp(-V*FF/(RR*T))));
        IpCa = GpCa * Cai / (KpCa + Cai);
        IpK = GpK * (V-EK) / (1.0 + exp((25.0-V)/5.98));
        IbNa = GbNa * (V - ENa);
        IbCa = GbCa * (V - ECa);

        Nai = Nai + dt * (-(INa + IbNa + 3.0*(INaK+INaCa)) / (Vc*FF)) * capacitance;
        Ki = Ki + dt * (-(IK1 + Ito + IKr + IKs + (-2.0*INaK) + IpK + Istim) / (Vc*FF)) * capacitance;

        kcasr = maxsr - ((maxsr - minsr) / (1.0 + (EC/CaSR)*(EC/CaSR)));
        k1 = k1p / kcasr;
        k2 = k2p * kcasr;
        Rhat = Rhat + dt * (k4 * (1.0 - Rhat) - k2 * CaSS * Rhat);
        O = (k1 * CaSS * CaSS * Rhat) / (k3 + k1 * CaSS * CaSS);

        Ileak = Vleak * (CaSR - Cai);
        Iup = Vmaxup / (1.0 + (Kup * Kup)/(Cai * Cai));
        Irel = Vrel * O * (CaSR - CaSS);
        Ixfer = Vxfer * (CaSS - Cai);

        Caibufc = (Cai * Bufc) / (Cai + Kbufc);
        dCaitotal = dt * (-((IbCa + IpCa - 2.0*INaCa) * capacitance)/(2.0*Vc*FF) + (Vsr/Vc)*(Ileak-Iup) + Ixfer);
        bi = Bufc - Caibufc - dCaitotal - Cai + Kbufc;
        ci = Kbufc * (Caibufc + dCaitotal + Cai);
        Cai = (sqrt(bi*bi + 4.0 * ci) - bi) / 2.0;

        Casrbufsr = (CaSR * Bufsr) / (CaSR + Kbufsr);
        dCaSRtotal = dt * (Iup - Ileak - Irel);
        bsr = Bufsr - Casrbufsr - dCaSRtotal - CaSR + Kbufsr;
        csr = Kbufsr * (Casrbufsr + dCaSRtotal + CaSR);
        CaSR = (sqrt(bsr*bsr + 4.0 * csr) - bsr) / 2.0;

        Cassbufss = (CaSS * Bufss) / (CaSS + Kbufss);
        dCaSStotal = dt * (-((ICaL * capacitance) / (2.0 * Vss * FF)) + (Vsr/Vss)*Irel - (Vc/Vss)*Ixfer);
        bss = Bufss - Cassbufss - dCaSStotal - CaSS + Kbufss;
        css = Kbufss * (Cassbufss + dCaSStotal + CaSS);
        CaSS = (sqrt(bss*bss + 4.0 * css) - bss) / 2.0;

        float prev_u = V;
        V = V - dt * (INa + IK1 + Ito + IKr + IKs + ICaL + INaCa + INaK + IpCa + IpK + IbCa + IbNa + Istim) / Cm;
        float u = data_type == 2 ? Cai : V;

        if (step_count > pre_pace_steps) {
            // APD only mode
            if (data_type == 1) {
                if (!activated && u > apd_thresh) {
                    activated = true;
                    float x0 = float((step_count-1))*dt;
                    float x1 = float(step_count)*dt;

                    float y0 = prev_u;
                    float y1 = u;

                    // Linear interpolation of actual crossing of threshold
                    APD_start = (x0*(y1 - apd_thresh) + x1*(apd_thresh - y0)) / (y1-y0);
                } else if (activated && u < apd_thresh) {
                    activated = false;

                    float x0 = float((step_count-1))*dt;
                    float x1 = float(step_count)*dt;

                    float y0 = prev_u;
                    float y1 = u;

                    // Linear interpolation of actual crossing of threshold
                    APD_end = (x0*(y1 - apd_thresh) + x1*(apd_thresh - y0)) / (y1-y0);
                    float sim_APD = APD_end - APD_start;
                    float target_APD = texelFetch(data_texture, ivec2(data_index++, 0), 0).r;
                    error += (target_APD - sim_APD) * (target_APD - sim_APD);
                    compared_points += 1;
                }
            }
            // Curve error only mode
            else {
                if (!first_align_upstroke && u > align_thresh) {
                    first_align_upstroke = true;
                    start_comp = step_count;
                    error = 0.0;
                }
                // Measure curve error
                if (first_align_upstroke && mod(float(step_count - start_comp), compare_stride) == 0.0) {
                    float actual = texelFetch(data_texture, ivec2(data_index++, 0), 0).r;
                    error += (u - actual)*(u - actual);
                    compared_points += 1;
                }
            }
        }

        // Save time series data for plotting
        if (float((step_count - pre_pace_steps) - 1) / float((num_steps - pre_pace_steps) - 1) <= cc.x) {
            saved_value = u;
        }
    }

    // TODO: We're not actually doing RMSE and never have been. This shouldn't
    // affect execution at all because we're only ever using error in
    // relative terms of "better" or "worse" and switching to RMSE shouldn't
    // change that. I'm still not making the change at present, because I
    // want to minimize the amount of things we need to debug at a time.

    // If we had missing activations to 2:1 blocking, I'm just adding the
    // missing APDs as raw error, since this is a very bad solution
    if (data_type == 1) {
        // While there are still leftover target APDs we never matched in the
        // simulation, add them as raw error.
        for (; data_index < num_data_points; data_index++) {
            float missing_APD = texelFetch(data_texture, ivec2(data_index, 0), 0).r;
            error += missing_APD*missing_APD;
            compared_points += 1;
        }
    }

    // Check to make sure that the simulation output was compared with most of the data. This check
    // assumes that no more that half of the first beat should be skipped due to alignment.
    //
    // This check prevents a very specific issue. Certain data traces (mostly downsampled model
    // data) will have such a fast upstroke that alignment occurs with a very large voltage
    // value. The simulation can achieve a very low error by having several beats with a peak just
    // below the maximum and then a later peak that hits it (for example, through physiologically
    // incorrect conductances for later currents or amplitude alternans) which results in very few
    // data points actually being compared.
    int points_to_compare = min(num_data_points, int(floor(((float(num_beats)-0.5)*period)/sample_interval)));
    if (data_type != 1 && data_index < points_to_compare) {
        error += 1e6;
    }

    error_texture = vec4(error, saved_value, 0, compared_points == 0 ? weight : weight / float(compared_points));
}
