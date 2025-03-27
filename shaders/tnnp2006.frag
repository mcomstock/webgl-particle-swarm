#version 300 es

precision highp float;
precision highp int;

// EPI, MYO, or ENDO
#define ENDO

uniform sampler2D in_particles_1, data_texture;

layout (location = 0) out vec4 error_texture;

in vec2 cc;

uniform float dt, period;
uniform int num_beats, pre_beats, data_type, err_type;
uniform float align_thresh;
uniform float sample_interval, apd_thresh, weight;
uniform float stim_dur, stim_mag, stim_offset_1, stim_offset_2, stim_t_scale;
uniform bool stim_biphasic;

uniform sampler2D table;
uniform int table_shift, table_npoints;
uniform float table_vmin, table_vmax, table_vekmin, table_vekmax;

const float RR = 8314.3;
const float T = 310.0;
const float FF = 96486.7;
const float SS = 0.2;
const float rho = 162.0;
const float Vc = 0.016404;
const float Vsr = 0.001094;
const float Vss = 0.00005468;
const float Ko = 5.4;
const float Nao = 140.0;
const float Cao = 2.0;
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
const float pKNa = 0.03;
// float GCaL = 3.98e-5;
// float kNaCa = 1000.0;
const float gamma = 0.35;
const float KmCa = 1.38;
const float KmNai = 87.5;
const float ksat = 0.1;
const float alpha = 2.5;
// float PNaK = 2.724;
const float KmK = 1.0;
const float KmNa = 40.0;
// float GpK = 0.0146;
// float GpCa = 0.1238;
const float KpCa = 0.0005;
// float GbNa = 0.00029;
// float GbCa = 0.000592;
const float Vmaxup = 0.006375;
const float Kup = 0.00025;
const float Vrel = 0.102;
const float k1p = 0.15;
const float k2p = 0.045;
const float k3 = 0.06;
const float k4 = 0.005;
const float EC = 1.5;
const float maxsr = 2.5;
const float minsr = 1.0;
const float Vleak = 0.00036;
const float Vxfer = 0.0038;
const float Bufc = 0.2;
const float Kbufc = 0.001;
const float Bufsr = 10.0;
const float Kbufsr = 0.3;
const float Bufss = 0.4;
const float Kbufss = 0.00025;
const float capacitance = 0.185;

const float Ko54sqrt = sqrt(Ko/5.4);
const float rtof = RR*T/FF;
const float invvcf = 1.0 / (Vc*FF);
const float vsrovc = Vsr/Vc;
const float vsrovss = Vsr/Vss;
const float vcovss = Vc/Vss;
const float invvssf = 1.0 / (Vss*FF);
const float ecsq = EC * EC;

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
    float tau_fcass;

    float Iion, INa, IK1, Ito, IKr, IKs, ICaL, INaCa, INaK, IpCa, IpK, IbCa, IbNa, Istim;
    float Ileak, Iup, Irel, Ixfer;

    float ENa, EK, ECa, EKs;

    float Caibufc, Casrbufsr, Cassbufss, dCaitotal, dCaSRtotal, dCaSStotal, bi, ci, bsr, csr, bss, css, k1, k2, kcasr, O;

    float time;

    float Kibase, Kidiff;
    float Naibase, Naidiff;

    Kibase = 136.0;
    Naibase = 9.413;
#if defined EPI
    V = -85.46;
    Rhat = 0.9891;
    // Nai = 9.293;
    Naidiff = 9.293 - Naibase;
    // Ki = 136.2;
    Kidiff = 0.2;
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
    // Nai = 9.322;
    Naidiff = 9.322 - Naibase;
    // Ki = 136.0;
    Kidiff = 0.0;
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
    // Nai = 9.413;
    Naidiff = 9.413 - Naibase;
    // Ki = 136.1;
    Kidiff = 0.1;
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

    float vidxint, vekidxint, wv1, wv2, wvek1, wvek2;
    int vidx1, vidx2, vekidx1, vekidx2, table_idx1, table_idx2;
    ivec2 table_idx1_2d, table_idx2_2d;
    vec4 table_val1, table_val2;

    float tau_m_exp, tau_h_exp, tau_j_exp, tau_xs_exp, tau_d_exp, tau_f_exp, tau_f2_exp, tau_r_exp, tau_s_exp, tau_xr1_exp, tau_xr2_exp;
    float ICaL_CaSS_coeff, ICaL_Cao, IpK_coeff, INaCa_Nai_coeff, INaCa_Cai_coeff, INaK_coeff;

    float u, prev_u;
    float save_ca = float(data_type == 2);
    float save_v = 1.0 - save_ca;

    float stim, stim_t;

    int table_mask = (1 << table_shift) - 1;
    float invvrange = float(table_npoints) / (table_vmax - table_vmin);
    float invvekrange = float(table_npoints) / (table_vekmax - table_vekmin);

    for (int step_count = 1; step_count <= num_steps; ++step_count) {
        /*
         * Table retrieval
         */

        // Indexing
        vidxint = invvrange * (V - table_vmin);
        vidx1 = int(round(vidxint));

#define SET_TABLE_VALS(SEQ, IDX1) {                             \
            table_idx1 = (SEQ) * table_npoints + (IDX1);        \
            table_idx1_2d[0] = table_idx1 & table_mask;         \
            table_idx1_2d[1] = table_idx1 >> table_shift;       \
            table_val1 = texelFetch(table, table_idx1_2d, 0);   \
        }

        SET_TABLE_VALS(0, vidx1);
        m_inf = table_val1[0];
        tau_m_exp = table_val1[1];
        h_inf = table_val1[2];
        tau_h_exp = table_val1[3];

        SET_TABLE_VALS(1, vidx1);
        j_inf = table_val1[0];
        tau_j_exp = table_val1[1];
        xs_inf = table_val1[2];
        tau_xs_exp = table_val1[3];

        SET_TABLE_VALS(2, vidx1);
        d_inf = table_val1[0];
        tau_d_exp = table_val1[1];
        f_inf = table_val1[2];
        tau_f_exp = table_val1[3];

        SET_TABLE_VALS(3, vidx1);
        f2_inf = table_val1[0];
        tau_f2_exp = table_val1[1];
        r_inf = table_val1[2];
        tau_r_exp = table_val1[3];

        SET_TABLE_VALS(4, vidx1);
        s_inf = table_val1[0];
        tau_s_exp = table_val1[1];
        xr1_inf = table_val1[2];
        tau_xr1_exp = table_val1[3];

        SET_TABLE_VALS(5, vidx1);
        xr2_inf = table_val1[0];
        tau_xr2_exp = table_val1[1];
        ICaL_CaSS_coeff = table_val1[2];
        ICaL_Cao = table_val1[3];

        SET_TABLE_VALS(6, vidx1);
        IpK_coeff = table_val1[0];
        INaCa_Nai_coeff = table_val1[1];
        INaCa_Cai_coeff = table_val1[2];
        INaK_coeff = table_val1[3];

        /*
         * End table retrieval
         */

        stim = 0.0;
        stim_t = mod(float(step_count)*dt, period);
        if (stim_t < stim_dur) {
            if (stim_biphasic) {
                stim = biphasic_stim_f(stim_t);
            } else {
                stim = square_stim_f(stim_t);
            }
        }

        Istim = stim;

        // TODO check performance
        if (abs(Kidiff) > 1.0) {
            Kibase = Kibase + Kidiff;
            Kidiff = 0.0;
        }
        Ki = Kibase + Kidiff;

        if (abs(Naidiff) > 0.1) {
            Naibase = Naibase + Naidiff;
            Naidiff = 0.0;
        }
        Nai = Naibase + Naidiff;

        m = m_inf - (m_inf - m) * tau_m_exp;
        h = h_inf - (h_inf - h) * tau_h_exp;
        j = j_inf - (j_inf - j) * tau_j_exp;
        xs = xs_inf - (xs_inf - xs) * tau_xs_exp;
        d = d_inf - (d_inf - d) * tau_d_exp;
        f = f_inf - (f_inf - f) * tau_f_exp;
        f2 = f2_inf - (f2_inf - f2) * tau_f2_exp;

        fcass_inf = 0.6 / (1.0 + 400.0*CaSS*CaSS) + 0.4;
        tau_fcass = 80.0 / (1.0 + 400.0*CaSS*CaSS) + 2.0;
        fcass = fcass_inf - (fcass_inf - fcass) * exp(-dt/tau_fcass);

        r = r_inf - (r_inf - r) * tau_r_exp;
        s = s_inf - (s_inf - s) * tau_s_exp;
        xr1 = xr1_inf - (xr1_inf - xr1) * tau_xr1_exp;
        xr2 = xr2_inf - (xr2_inf - xr2) * tau_xr2_exp;

        ENa = rtof * log(Nao/Nai);
        EK = rtof * log(Ko/Ki);
        ECa = 0.5*rtof * log(Cao/Cai);
        EKs = rtof * log((Ko+pKNa*Nao)/(Ki+pKNa*Nai));

        /*
         * V-EK table retrieval (must occur after EK is computed)
         */

        // Indexing
        vekidxint = invvekrange * ((V - EK) - table_vekmin);
        vekidx1 = int(round(vekidxint));

        // Lookup
        SET_TABLE_VALS(7, vekidx1);
        xK1_inf = table_val1[0];

        /*
         * End V-EK table retrieval
         */

        ICaL = GCaL * d * f * f2 * fcass * (ICaL_CaSS_coeff * CaSS - ICaL_Cao);

        IKs = GKs * xs * xs * (V - EKs);
        INa = GNa * m * m * m * h * j * (V - ENa);
        Ito = Gto * r * s * (V - EK);
        IKr = GKr * Ko54sqrt * xr1 * xr2 * (V - EK);

        IK1 = GK1 * Ko54sqrt * xK1_inf * (V-EK);

        INaCa = kNaCa * (INaCa_Nai_coeff*Nai*Nai*Nai - INaCa_Cai_coeff*Cai);
        INaK = PNaK*Ko*Nai / ((Ko+KmK) * (Nai+KmNa) * (1.0 + 0.1245*exp(-0.1*V*FF/(RR*T)) + 0.0353 * exp(-V*FF/(RR*T))));
        IpCa = GpCa * Cai / (KpCa + Cai);
        IpK = GpK * (V-EK) * IpK_coeff;
        IbNa = GbNa * (V - ENa);
        IbCa = GbCa * (V - ECa);

        // Track differences to avoid precision issue of adding small differences to large values
        Naidiff = Naidiff + dt * (-(INa + IbNa + 3.0*(INaK+INaCa)) * invvcf) * capacitance;
        Kidiff = Kidiff + dt * (-(IK1 + Ito + IKr + IKs + (-2.0*INaK) + IpK + Istim) * invvcf) * capacitance;

        kcasr = maxsr - ((maxsr - minsr) / (1.0 + ecsq/(CaSR*CaSR)));
        k1 = k1p / kcasr;
        k2 = k2p * kcasr;
        Rhat = Rhat + dt * (k4 * (1.0 - Rhat) - k2 * CaSS * Rhat);
        O = (k1 * CaSS * CaSS * Rhat) / (k3 + k1 * CaSS * CaSS);

        Ileak = Vleak * (CaSR - Cai);
        Iup = Kup / Cai;
        Iup = Vmaxup / (1.0 + Iup*Iup);
        Irel = Vrel * O * (CaSR - CaSS);
        Ixfer = Vxfer * (CaSS - Cai);

        Caibufc = (Cai * Bufc) / (Cai + Kbufc);
        dCaitotal = dt * (-((IbCa + IpCa - 2.0*INaCa) * capacitance) * 0.5 * invvcf + vsrovc*(Ileak-Iup) + Ixfer);
        bi = Bufc - Caibufc - dCaitotal - Cai + Kbufc;
        ci = Kbufc * (Caibufc + dCaitotal + Cai);
        // Use a different form of the equation to avoid cancelation for very small bi/ci values
        Cai = 4.0 * (ci / (bi * bi));
        Cai = (0.5 * bi * Cai) / (sqrt(1.0 + Cai) + 1.0);

        Casrbufsr = (CaSR * Bufsr) / (CaSR + Kbufsr);
        dCaSRtotal = dt * (Iup - Ileak - Irel);
        bsr = Bufsr - Casrbufsr - dCaSRtotal - CaSR + Kbufsr;
        csr = Kbufsr * (Casrbufsr + dCaSRtotal + CaSR);
        CaSR = (sqrt(bsr*bsr + 4.0 * csr) - bsr) * 0.5;

        Cassbufss = (CaSS * Bufss) / (CaSS + Kbufss);
        dCaSStotal = dt * (-((ICaL * capacitance) * 0.5 * invvssf) + vsrovss*Irel - vcovss*Ixfer);
        bss = Bufss - Cassbufss - dCaSStotal - CaSS + Kbufss;
        css = Kbufss * (Cassbufss + dCaSStotal + CaSS);
        CaSS = (sqrt(bss*bss + 4.0 * css) - bss) * 0.5;

        prev_u = V;
        V = V - dt * (INa + IK1 + Ito + IKr + IKs + ICaL + INaCa + INaK + IpCa + IpK + IbCa + IbNa + Istim);
        u = save_ca * Cai + save_v * V;

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
                    error += err_type == 1 ? abs(target_APD - sim_APD) : (target_APD - sim_APD) * (target_APD - sim_APD);
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
                    error += err_type == 1 ? abs(u - actual) : (u - actual) * (u - actual);
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
            error += err_type == 1 ? missing_APD : missing_APD * missing_APD;
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
