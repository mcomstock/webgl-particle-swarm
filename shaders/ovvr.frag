#version 300 es

precision highp float;
precision highp int;

uniform sampler2D in_particles_1, in_particles_2, data_texture;

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

// Model parameters
const float RR = 8314.0;
const float TT = 310.0;
const float FF = 96485.0;
const float RTF = (RR*TT)/FF;
const float FRT = FF/(RR*TT);
const float FFRT = (FF*FF)/(RR*TT);
const float C_M = 1.0;
const float R_CG = 2.0;

const float Na_o = 140.0;
const float Ca_o = 1.8;
const float K_o = 5.4;

const float sqrtko54 = sqrt(K_o/5.4);
const float sqrtko = sqrt(K_o);

const float PRNaK = 0.01833;

const float pi = 4.0 * atan(1.0);
const float len = 0.01;
const float radius = 0.0011;
// The vcell units need to be converted, even though they don't tell you
const float vcell = 1000.0 * pi * radius * radius * len;
const float Ageo = 2.0 * pi * radius * radius + 2.0 * pi * radius * len;
const float Acap = R_CG * Ageo;
const float vmyo = 0.68 * vcell;
const float vnsr = 0.0552 * vcell;
const float vjsr = 0.0048 * vcell;
const float vss = 0.02 * vcell;

const float Ahfast = 0.99;
const float Ahslow = 0.01;
const float AhCaMKfast = Ahfast;
const float AhCaMKslow = Ahslow;
const float KmCaMK = 0.15;
// const float GNafastbar = 75.0;
const float tauhL = 200.0;
const float tauhLCaMK = 3.0 * tauhL;
// const float GNalatebar = 0.0075;
// const float Gtobar = 0.02;
const float Affast = 0.6;
const float Afslow = 1.0 - Affast;
const float taujCa = 75.0;
const float AfCaMKfast = Affast;
const float AfCaMKslow = Afslow;
const float AfCaCaMKfast = AfCaMKfast;
const float AfCaCaMKslow = AfCaMKslow;
const float Kmn = 0.002;
const float kp2n = 1000.0;
// const float PCa = 0.0001;
const float gammaCai = 1.0;
const float gammaCao = 0.341;
const float zCa = 2.0;
// const float PCaNa = 0.00125 * PCa;
const float gammaNai = 0.75;
const float gammaNao = 0.75;
const float zNa = 1.0;
// const float PCaK = 3.574e-4 * PCa;
const float gammaKi = 0.75;
const float gammaKo = 0.75;
const float zK = 1.0;
// const float PCaCaMK = 1.1 * PCa;
// const float PCaNaCaMK = 0.00125 * PCaCaMK;
// const float PCaKCaMK = 3.574e-4 * PCaCaMK;
// const float GKrbar = 0.046;
// const float GKsbar = 0.0034;
// const float GK1bar = 0.1908;
const float kNa1 = 15.0;
const float kNa2 = 5.0;
const float kNa3 = 88.12;
const float kasymm = 12.5;
const float omegaNa = 6e4;
const float omegaCa = 6e4;
const float omegaNaCa = 5e3;
const float kCaon = 1.5e6;
const float kCaoff = 5e3;
const float qNa = 0.5224;
const float qCa = 0.1670;
const float h10 = kasymm + 1.0 + (Na_o/kNa1) * (1.0 + Na_o/kNa2);
const float h11 = (Na_o*Na_o) / (h10 * kNa1 * kNa2);
const float h12 = 1.0 / h10;
const float k2 = kCaoff;
const float k5 = kCaoff;
const float KmCaAct = 150e-6;
// const float GNaCabar = 0.0008;
const float k1pNaK = 949.5;
const float k1mNaK = 182.4;
const float k2pNaK = 687.2;
const float k2mNaK = 39.4;
const float k3pNaK = 1899.0;
const float k3mNaK = 79300.0;
const float k4pNaK = 639.0;
const float k4mNaK = 40.0;
const float KNaio = 9.073;
const float KNaoo = 27.78;
const float Delta = -0.1550;
const float KKi = 0.5;
const float KKo = 0.3582;
const float MgADP = 0.05;
const float MgATP = 9.8;
const float KMgATP = 1.698e-7;
const float Hp = 1e-7;
const float SigmaP = 4.2;
const float KHP = 1.698e-7;
const float KNaP = 224.0;
const float KKP = 292.0;
const float beta1 = k1mNaK * MgADP;
const float alpha2 = k2pNaK;
const float alpha4 = (k4pNaK*(MgATP/KMgATP))/(1.0+MgATP/KMgATP);
const float oneKoKKo2 = pow((1.0+K_o/KKo), 2.0);
const float KoKKo2 = pow((K_o/KKo), 2.0);
// const float PNab = 3.75e-10;
// const float PCab = 2.5e-8;
// const float GKbbar = 0.003;
// const float GpCabar = 0.0005;
const float alphaCaMK = 0.05;
const float betaCaMK = 0.00068;
const float CaMK0 = 0.05;
const float KmCaM = 0.0015;
const float taudiffNa = 2.0;
const float taudiffK = 2.0;
const float taudiffCa = 0.2;
const float betatau = 4.75;
const float alpharel = 0.5 * betatau;
const float betatauCaMK = 1.25 * betatau;
const float alpharelCaMK = 0.5 * betatauCaMK;
const float DeltaKmPLBbar = 0.00017;
const float DeltaJupCaMKbar = 1.75;
const float tautr = 100.0;
const float CMDN = 0.05;
const float TRPN = 0.07;
const float BSR = 0.047;
const float BSL = 1.124;
const float CSQN = 10.0;
const float KmCMDN = 0.00238;
const float KmTRPN = 0.0005;
const float KmBSR = 0.00087;
const float KmBSL = 0.0087;
const float KmCSQN = 0.8;

const float invKKi = 1.0/KKi;
const float b3denom = 1.0 / (1.0 + MgATP/KMgATP);
const float inv15 = 1.0 / 15.0;
const float invtautr = 1.0 / tautr;
const float invtaudiffna = 1.0 / taudiffNa;
const float invtaudiffca = 1.0 / taudiffCa;
const float invtaudiffk = 1.0 / taudiffK;
const float invkNa1 = 1.0 / kNa1;
const float invkNa2 = 1.0 / kNa2;
const float invkNa3 = 1.0 / kNa3;
const float invKHP = 1.0 / KHP;
const float invKNaP = 1.0 / KNaP;
const float invKKP = 1.0 / KKP;

const float zcafgcai = zCa * FF * gammaCai;
const float zcafgcaocao = zCa * FF * gammaCao * Ca_o;
const float znafgnai = zNa * FF * gammaNai;
const float znafgnaonao = zCa * FF * gammaNao * Na_o;
const float zkfgki = zK * FF * gammaKi;
const float zkfgkoko = zK * FF * gammaKo * K_o;
const float zca2frtgcai = zCa * zCa * FFRT * gammaCai;
const float zna2frtgnai = zNa * zNa * FFRT * gammaNai;
const float zk2frtgki = zK * zK * FFRT * gammaKi;
const float zca2frtgcaocao = zCa * zCa * FFRT * gammaCao * Ca_o;
const float zna2frtgnaonao = zNa * zNa * FFRT * gammaNao * Na_o;
const float zk2frtgkoko = zK * zK * FFRT * gammaKo * K_o;

const float acapofvmyo = Acap / (FF*vmyo);
const float acapo2fvmyo = Acap / (2.0*FF*vmyo);
const float vssovmyo = vss / vmyo;
const float acapofvss = Acap / (FF*vss);
const float acapo2fvss = Acap / (2.0*FF*vss);
const float vnsrovmyo = vnsr / vmyo;
const float vjsrovss = vjsr / vss;
const float vjsrovnsr = vjsr / vnsr;
const float cmdnkmcmdn = CMDN * KmCMDN;
const float trpnkmtrpn = TRPN * KmTRPN;
const float bsrkmbsr = BSR * KmBSR;
const float bslkmbsl = BSL * KmBSL;
const float csqnkmcsqn = CSQN * KmCSQN;

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

    vec4 particles_5 = texelFetch(in_particles_2, idx, 0);
    vec4 particles_6 = texelFetch(in_particles_2, idx + ivec2(tex_size.x/2, 0), 0);
    vec4 particles_7 = texelFetch(in_particles_2, idx + ivec2(0, tex_size.y/2), 0);
    vec4 particles_8 = texelFetch(in_particles_2, idx + ivec2(tex_size.x/2, tex_size.y/2), 0);

    float GNafastbar = particles_1[0];
    float GNalatebar = particles_1[1];
    float Gtobar = particles_1[2];
    float PCa = particles_1[3];

    float PCaNa = particles_2[0];
    float PCaK = particles_2[1];
    float PCaCaMK = particles_2[2];
    float PCaNaCaMK = particles_2[3];

    float PCaKCaMK = particles_3[0];
    float GKrbar = particles_3[1];
    float GKsbar = particles_3[2];
    float GK1bar = particles_3[3];

    float GNaCabar = particles_4[0];
    // Note: not a parameter in the model
    float GNaKbar = particles_4[1];
    float PNab = particles_4[2];
    float PCab = particles_4[3];

    float GKbbar = particles_5[0];
    float GpCabar = particles_5[1];

    float INa, Ito, ICaL, ICaNa, ICaK, IKr, IKs, IK1, INaCai, INaCass, INaCa, INaK, INab, ICab, IKb, IpCa;
    float Istim;

    // Reversal potentials
    float ENa, EK, EKs;

    // State variables
    float V;
    float Na_i, Na_ss, K_i, K_ss, Ca_i, Ca_ss, Ca_nsr, Ca_jsr;
    float m, hfast, hslow, j, hCaMKslow, jCaMK, mL, hL, hLCaMK, a, ifast, islow, aCaMK, iCaMKfast, iCaMKslow, d, ffast;
    float fslow, fCafast, fCaslow, jCa, n, fCaMKfast, fCaCaMKfast, xrfast, xrslow, xs1, xs2, xK1, JrelNP, JrelCaMK, CaMKtrap;

    // Intermediate variables
    float time;
    float minf, hinf, h;
    float hCaMKinf, hCaMKfast, hCaMK, OINaCaMK, INafast, mLinf, hLinf;
    float hLCaMKinf, OINaLCaMK, INalate;
    float ainf, iinf, Aifast, Aislow, i, aCaMKinf, deltaCaMKdevelop;
    float deltaCaMKrecover, AiCaMKfast, AiCaMKslow, iCaMK, OItoCaMK;
    float dinf, finf, f, AfCafast, AfCaslow, fCa;
    float fCaMK, fCaCaMK, km2n, alphan, PsiCa;
    float ICaLbar, PsiCaNa, ICaNabar, PsiCaK, ICaKbar, ICaLCaMKbar, ICaNaCaMKbar, ICaKCaMKbar, OICaLCaMK;
    float xrinf, Axrfast, Axrslow, xr, RKr;
    float xs1inf, xs2inf;
    float xK1inf, RK1;
    float hCa, hNa, h1i, h1ss, h2i, h2ss, h3i, h3ss, h4i, h4ss, h5i, h5ss, h6i, h6ss, h7, h8, h9, k1, k3p, k3pp;
    float k3, k4pi, k4pss, k4ppi, k4ppss, k4i, k4ss, k6i, k6ss, k7i, k7ss, k8, x1i, x1ss, x2i, x2ss, x3i, x3ss, x4i, x4ss, E1i;
    float E1ss, E2i, E2ss, E3i, E3ss, E4i, E4ss, alloi, alloss, JNaCaNai, JNaCaNass, JNaCaCai, JNaCaCass;
    float KNai, KNao, P, alpha1, beta2, alpha3, beta3, beta4, x1NaK, x2NaK, x3NaK, x4NaK, E1NaK, E2NaK, E3NaK, E4NaK, JNaKNa;
    float JNaKK, JNaK, NaiKNai3, oneNaiKNai3, oneKiKKi2, NaoKNao3, oneNaoKNao3, oneNoKKo2, KiKKi2;
    float xKb;
    float CaMKbound, CaMKactive;
    float JdiffNa, JdiffCa, JdiffK;
    float JrelNPinf, taurelNP, JrelCaMKinf, taurelCaMK, OrelCaMK, Jrel;
    float JupNP, JupCaMK, OupCaMK, Jleak, Jup;
    float Jtr;
    float betaCai_inv, betaCass_inv, betaCajsr_inv;
    float bc1, bc2;

    float taum_exp, tauhfast_exp, tauhslow_exp, tauj_exp, tauhCaMKslow_exp, taua_exp, tauifast_exp,
        tauislow_exp, tauiCaMKfast_exp, tauiCaMKslow_exp, taud_exp, tauffast_exp, taufslow_exp,
        taufCafast_exp, taufCaslow_exp, taufCaMKfast_exp, taufCaCaMKfast_exp, exp_zcavfrt,
        exp_znavfrt, exp_zkvfrt, tauxrfast_exp, tauxrslow_exp, tauxs1_exp, tauxs2_exp, tauxK1_exp,
        INab_coeff, ICab_coeff, taujCaMK_exp, invKNai, invKNao, invhCa;

    float a1b4denom, b2a3denom, pow_cajsr;
    float invxis, invxsss, invnaks;

    float tauhL_exp = exp(-dt/tauhL);
    float tauhLCaMK_exp = exp(-dt/tauhLCaMK);
    float taujCa_exp = exp(-dt/taujCa);

    float K_i_base, K_i_diff;
    float K_ss_base, K_ss_diff;
    float Na_i_base, Na_i_diff;
    float Na_ss_base, Na_ss_diff;

    Na_i_base = 7.23;
    Na_i_diff = 0.0;
    Na_ss_base = 7.23;
    Na_ss_diff = 0.0;
    K_i_base = 143.79;
    K_i_diff = 0.0;
    K_ss_base = 143.79;
    K_ss_diff = 0.0;

    V = -87.84;
    // Na_i = 7.23;
    // Na_ss = 7.23;
    // K_i = 143.79;
    // K_ss = 143.79;
    Ca_i = 8.54e-5;
    Ca_ss = 8.43e-5;
    Ca_nsr = 1.61;
    Ca_jsr = 1.56;
    m = 0.0074621;
    hfast = 0.692591;
    hslow = 0.692574;
    j = 0.692477;
    hCaMKslow = 0.448501;
    jCaMK = 0.692413;
    mL = 0.000194015;
    hL = 0.496116;
    hLCaMK = 0.265885;
    a = 0.00101185;
    ifast = 0.999542;
    islow = 0.589579;
    aCaMK = 0.000515567;
    iCaMKfast = 0.999542;
    iCaMKslow = 0.641861;
    d = 2.43015e-9;
    ffast = 1.0;
    fslow = 0.910671;
    fCafast = 1.0;
    fCaslow = 0.99982;
    jCa = 0.999977;
    n = 0.00267171;
    fCaMKfast = 1.0;
    fCaCaMKfast = 1.0;
    xrfast = 8.26608e-6;
    xrslow = 0.453268;
    xs1 = 0.270492;
    xs2 = 0.0001963;
    xK1 = 0.996801;
    JrelNP = 2.53943e-5;
    JrelCaMK = 3.17262e-7;
    CaMKtrap = 0.0124065;

    float compare_stride = round(sample_interval / dt);

    float error = (data_type == 1) ? 0.0 : 10000000000.0;

    int data_index = 0;

    int start_comp = 0;
    int compared_points = 0;
    bool first_align_upstroke = false;
    float saved_value = -1.0;

    float APD_start, APD_end;

    bool activated = false;

    float vidxint;
    int vidx1, table_idx1;
    ivec2 table_idx1_2d;
    vec4 table_val1;

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

        vidxint = invvrange * (V - table_vmin);
        vidx1 = int(round(vidxint));

#define SET_TABLE_VALS(SEQ, IDX1) {                             \
            table_idx1 = (SEQ) * table_npoints + (IDX1);        \
            table_idx1_2d[0] = table_idx1 & table_mask;         \
            table_idx1_2d[1] = table_idx1 >> table_shift;       \
            table_val1 = texelFetch(table, table_idx1_2d, 0);   \
        }

        SET_TABLE_VALS(0, vidx1);
        minf = table_val1[0];
        taum_exp = table_val1[1];
        hinf = table_val1[2];
        tauhfast_exp = table_val1[3];

        SET_TABLE_VALS(1, vidx1);
        tauhslow_exp = table_val1[0];
        tauj_exp = table_val1[1];
        hCaMKinf = table_val1[2];
        tauhCaMKslow_exp = table_val1[3];

        SET_TABLE_VALS(2, vidx1);
        mLinf = table_val1[0];
        hLinf = table_val1[1];
        hLCaMKinf = table_val1[2];
        ainf = table_val1[3];

        SET_TABLE_VALS(3, vidx1);
        taua_exp = table_val1[0];
        iinf = table_val1[1];
        tauifast_exp = table_val1[2];
        tauislow_exp = table_val1[3];

        SET_TABLE_VALS(4, vidx1);
        Aifast = table_val1[0];
        aCaMKinf = table_val1[1];
        tauiCaMKfast_exp = table_val1[2];
        tauiCaMKslow_exp = table_val1[3];

        SET_TABLE_VALS(5, vidx1);
        dinf = table_val1[0];
        taud_exp = table_val1[1];
        finf = table_val1[2];
        tauffast_exp = table_val1[3];

        SET_TABLE_VALS(6, vidx1);
        taufslow_exp = table_val1[0];
        taufCafast_exp = table_val1[1];
        taufCaslow_exp = table_val1[2];
        AfCafast = table_val1[3];

        SET_TABLE_VALS(7, vidx1);
        taufCaMKfast_exp = table_val1[0];
        taufCaCaMKfast_exp = table_val1[1];
        exp_zcavfrt = table_val1[2];
        exp_znavfrt = table_val1[3];

        SET_TABLE_VALS(8, vidx1);
        exp_zkvfrt = table_val1[0];
        xrinf = table_val1[1];
        tauxrfast_exp = table_val1[2];
        tauxrslow_exp = table_val1[3];

        SET_TABLE_VALS(9, vidx1);
        Axrfast = table_val1[0];
        RKr = table_val1[1];
        xs1inf = table_val1[2];
        tauxs1_exp = table_val1[3];

        SET_TABLE_VALS(10, vidx1);
        tauxs2_exp = table_val1[0];
        xK1inf = table_val1[1];
        tauxK1_exp = table_val1[2];
        RK1 = table_val1[3];

        SET_TABLE_VALS(11, vidx1);
        invhCa = table_val1[0];
        hNa = table_val1[1];
        invKNai = table_val1[2];
        invKNao = table_val1[3];

        SET_TABLE_VALS(12, vidx1);
        INab_coeff = table_val1[0];
        ICab_coeff = table_val1[1];
        xKb = table_val1[2];
        taujCaMK_exp = table_val1[3];

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

        if (abs(K_i_diff) > 0.1) {
            K_i_base = K_i_base + K_i_diff;
            K_i_diff = 0.0;
        }

        if (abs(K_ss_diff) > 0.1) {
            K_ss_base = K_ss_base + K_ss_diff;
            K_ss_diff = 0.0;
        }

        if (abs(Na_i_diff) > 0.1) {
            Na_i_base = Na_i_base + Na_i_diff;
            Na_i_diff = 0.0;
        }

        if (abs(Na_ss_diff) > 0.1) {
            Na_ss_base = Na_ss_base + Na_ss_diff;
            Na_ss_diff = 0.0;
        }

        K_i = K_i_base + K_i_diff;
        K_ss = K_ss_base + K_ss_diff;
        Na_i = Na_i_base + Na_i_diff;
        Na_ss = Na_ss_base + Na_ss_diff;

        //
        // Calcium/calmodulin-dependent protein kinase (CaMK)
        //

        CaMKbound = CaMK0 * (1.0 - CaMKtrap) / (1.0 + KmCaM/Ca_ss);
        CaMKactive = CaMKbound + CaMKtrap;

        //
        // Reversal potentials
        //

        ENa = RTF * log(Na_o/Na_i);
        EK = RTF * log(K_o/K_i);
        EKs = RTF * log((K_o + PRNaK * Na_o)/(K_i + PRNaK * Na_i));

        //
        // Sodium current (INa)
        //

        // Note: When comparing the INa current from this implementation against Figure 9 of the paper,
        // it appears that the peak is shifted by 1 ms, but is otherwise the same. The difference may be
        // due to how the stimulus is applied.

        m = minf - (minf - m) * taum_exp;

        hfast = hinf - (hinf - hfast) * tauhfast_exp;
        hslow = hinf - (hinf - hslow) * tauhslow_exp;
        h = Ahfast * hfast + Ahslow * hslow;

        j = hinf - (hinf - j) * tauj_exp;

        hCaMKslow = hCaMKinf - (hCaMKinf - hCaMKslow) * tauhCaMKslow_exp;
        hCaMK = AhCaMKfast * hfast + AhCaMKslow * hCaMKslow;

        jCaMK = hinf - (hinf - jCaMK) * taujCaMK_exp;

        OINaCaMK = 1.0 / (1.0 + KmCaMK/CaMKactive);

        INafast = GNafastbar * (V - ENa) * m*m*m * ((1.0 - OINaCaMK) * h * j + OINaCaMK * hCaMK * jCaMK);

        mL = mLinf - (mLinf - mL) * taum_exp;
        hL = hLinf - (hLinf - hL) * tauhL_exp;
        hLCaMK = hLCaMKinf - (hLCaMKinf - hLCaMK) * tauhLCaMK_exp;

        OINaLCaMK = OINaCaMK;

        INalate = GNalatebar * (V - ENa) * mL * ((1.0 - OINaLCaMK) * hL + OINaLCaMK * hLCaMK);

        INa = INafast + INalate;

        //
        // Transient outward potassium current (Ito)
        //

        a = ainf - (ainf - a) * taua_exp;

        ifast = iinf - (iinf - ifast) * tauifast_exp;
        islow = iinf - (iinf - islow) * tauislow_exp;
        Aislow = 1.0 - Aifast;
        i = Aifast * ifast + Aislow * islow;

        aCaMK = aCaMKinf - (aCaMKinf - aCaMK) * taua_exp;

        iCaMKfast = iinf - (iinf - iCaMKfast) * tauiCaMKfast_exp;
        iCaMKslow = iinf - (iinf - iCaMKslow) * tauiCaMKslow_exp;
        iCaMK = Aifast * iCaMKfast + Aislow * iCaMKslow;

        OItoCaMK = OINaCaMK;

        Ito = Gtobar * (V - EK) * ((1.0 - OItoCaMK) * a * i + OItoCaMK * aCaMK * iCaMK);

        //
        // L-type calcium current (ICaL)
        //

        d = dinf - (dinf - d) * taud_exp;

        ffast = finf - (finf - ffast) * tauffast_exp;
        fslow = finf - (finf - fslow) * taufslow_exp;
        f = Affast * ffast + Afslow * fslow;

        fCafast = finf - (finf - fCafast) * taufCafast_exp;
        fCaslow = finf - (finf - fCaslow) * taufCaslow_exp;
        AfCaslow = 1.0 - AfCafast;
        fCa = AfCafast * fCafast + AfCaslow * fCaslow;

        jCa = finf - (finf - jCa) * taujCa_exp;

        fCaMKfast = finf - (finf - fCaMKfast) * taufCaMKfast_exp;
        fCaMK = AfCaMKfast * fCaMKfast + AfCaMKslow * fslow;

        fCaCaMKfast = finf - (finf - fCaCaMKfast) * taufCaCaMKfast_exp;
        fCaCaMK = AfCaCaMKfast * fCaCaMKfast + AfCaCaMKslow * fCaslow;

        km2n = jCa;
        alphan = 1.0 + Kmn/Ca_ss;
        alphan = alphan * alphan;
        alphan = alphan * alphan;
        alphan = 1.0 / (kp2n/km2n + alphan);
        // TODO replace with lookup?
        n = alphan * (kp2n/km2n) - (alphan * (kp2n/km2n) - n) * exp(-km2n*dt);

        if (abs(V) < 0.01) {
            // L'Hopital
            PsiCa = zcafgcai * Ca_ss - zcafgcaocao;
            PsiCaNa = znafgnai * Na_ss - znafgnaonao;
            PsiCaK = zkfgki * K_ss - zkfgkoko;
        } else {
            PsiCa = V * (zca2frtgcai * Ca_ss * exp_zcavfrt - zca2frtgcaocao) / (exp_zcavfrt - 1.0);
            PsiCaNa = V * (zna2frtgnai * Na_ss * exp_znavfrt - zna2frtgnaonao) / (exp_znavfrt - 1.0);
            PsiCaK = V * (zk2frtgki * K_ss * exp_zkvfrt - zk2frtgkoko) / (exp_zkvfrt - 1.0);
        }

        ICaLbar = PCa * PsiCa;
        ICaNabar = PCaNa * PsiCaNa;
        ICaKbar = PCaK * PsiCaK;

        ICaLCaMKbar = PCaCaMK * PsiCa;
        ICaNaCaMKbar = PCaNaCaMK * PsiCaNa;
        ICaKCaMKbar = PCaKCaMK * PsiCaK;

        OICaLCaMK = OINaCaMK;

        ICaL = ICaLbar * d * (1.0 - OICaLCaMK) * (f * (1.0 - n) + fCa * n * jCa) + ICaLCaMKbar * d * OICaLCaMK * (fCaMK * (1.0 - n) + fCaCaMK * n * jCa);
        ICaNa = ICaNabar * d * (1.0 - OICaLCaMK) * (f * (1.0 - n) + fCa * n * jCa) + ICaNaCaMKbar * d * OICaLCaMK * (fCaMK * (1.0 - n) + fCaCaMK * n * jCa);
        ICaK = ICaKbar * d * (1.0 - OICaLCaMK) * (f * (1.0 - n) + fCa * n * jCa) + ICaKCaMKbar * d * OICaLCaMK * (fCaMK * (1.0 - n) + fCaCaMK * n * jCa);

        //
        // Rapid delayed rectifier potassium current (IKr)
        //

        xrfast = xrinf - (xrinf - xrfast) * tauxrfast_exp;
        xrslow = xrinf - (xrinf - xrslow) * tauxrslow_exp;
        Axrslow = 1.0 - Axrfast;
        xr = Axrfast * xrfast + Axrslow * xrslow;

        IKr = GKrbar * sqrtko54 * xr * RKr * (V - EK);

        //
        // Slow delayed rectifier potassium current (IKs)
        //

        xs1 = xs1inf - (xs1inf - xs1) * tauxs1_exp;
        xs2 = xs1inf - (xs1inf - xs2) * tauxs2_exp;

        // TODO: replace with table over Ca vals?
        IKs = GKsbar * (1.0 + 0.6 / (1.0 + pow((3.8e-5/Ca_i), 1.4))) * xs1 * xs2 * (V - EKs);

        //
        // Inward rectifier potassium current (IK1)
        //

        xK1 = xK1inf - (xK1inf - xK1) * tauxK1_exp;

        IK1 = GK1bar * sqrtko * xK1 * RK1 * (V - EK);

        //
        // Sodium/calcium exchange current (INaCa)
        //

        h1i = 1.0 + (Na_i*invkNa3) * (1.0 + hNa);
        h1ss = 1.0 + (Na_ss*invkNa3) * (1.0 + hNa);
        h2i = (Na_i*hNa)/(kNa3*h1i);
        h2ss = (Na_ss*hNa)/(kNa3*h1ss);
        h3i = 1.0/h1i;
        h3ss = 1.0/h1ss;
        h4i = 1.0 + (Na_i*invkNa1) * (1.0 + Na_i*invkNa2);
        h4ss = 1.0 + (Na_ss*invkNa1) * (1.0 + Na_ss*invkNa2);
        h5i = (Na_i*Na_i) / (h4i * kNa1 * kNa2);
        h5ss = (Na_ss*Na_ss) / (h4ss * kNa1 * kNa2);
        h6i = 1.0 / h4i;
        h6ss = 1.0 / h4ss;
        h7 = 1.0 + (Na_o*invkNa3) * (1.0 + 1.0/hNa);
        h8 = Na_o / (kNa3 * hNa * h7);
        h9 = 1.0 / h7;

        k1 = h12 * Ca_o * kCaon;
        k3p = h9 * omegaCa;
        k3pp = h8 * omegaNaCa;
        k3 = k3p + k3pp;
        k4pi = (h3i * omegaCa) * invhCa;
        k4pss = (h3ss * omegaCa) * invhCa;
        k4ppi = h2i * omegaNaCa;
        k4ppss = h2ss * omegaNaCa;
        k4i = k4pi + k4ppi;
        k4ss = k4pss + k4ppss;
        k6i = h6i * Ca_i * kCaon;
        k6ss = h6ss * Ca_ss * kCaon;
        k7i = h5i * h2i * omegaNa;
        k7ss = h5ss * h2ss * omegaNa;
        k8 = h8 * h11 * omegaNa;

        x1i = k2 * k4i * (k7i + k6i) + k5 * k7i * (k2 + k3);
        x1ss = k2 * k4ss * (k7ss + k6ss) + k5 * k7ss * (k2 + k3);
        x2i = k1 * k7i * (k4i + k5) + k4i * k6i * (k1 + k8);
        x2ss = k1 * k7ss * (k4ss + k5) + k4ss * k6ss * (k1 + k8);
        x3i = k1 * k3 * (k7i + k6i) + k8 * k6i * (k2 + k3);
        x3ss = k1 * k3 * (k7ss + k6ss) + k8 * k6ss * (k2 + k3);
        x4i = k2 * k8 * (k4i + k5) + k3 * k5 * (k1 + k8);
        x4ss = k2 * k8 * (k4ss + k5) + k3 * k5 * (k1 + k8);

        invxis = 1.0 / (x1i + x2i + x3i + x4i);
        invxsss = 1.0 / (x1ss + x2ss + x3ss + x4ss);

        E1i = x1i * invxis;
        E1ss = x1ss * invxsss;
        E2i = x2i * invxis;
        E2ss = x2ss * invxsss;
        E3i = x3i * invxis;
        E3ss = x3ss * invxsss;
        E4i = x4i * invxis;
        E4ss = x4ss * invxsss;

        alloi = KmCaAct / Ca_i;
        alloss = KmCaAct / Ca_ss;
        alloi = 1.0 / (1.0 + alloi*alloi);
        alloss = 1.0 / (1.0 + alloss*alloss);

        JNaCaNai = 3.0 * (E4i * k7i - E1i * k8) + E3i * k4ppi - E2i * k3pp;
        JNaCaNass = 3.0 * (E4ss * k7ss - E1ss * k8) + E3ss * k4ppss - E2ss * k3pp;
        JNaCaCai = E2i * k2 - E1i * k1;
        JNaCaCass = E2ss * k2 - E1ss * k1;

        INaCai = GNaCabar * 0.8 * alloi * (zNa * JNaCaNai + zCa * JNaCaCai);
        INaCass = GNaCabar * 0.2 * alloss * (zNa * JNaCaNass + zCa * JNaCaCass);
        INaCa = INaCai + INaCass;

        //
        // Sodium/potassium ATPase current (INaK)
        //

        P = SigmaP / (1.0 + Hp*invKHP + Na_i*invKNaP + K_i*invKKP);

        NaiKNai3 = Na_i * invKNai;
        NaiKNai3 = NaiKNai3 * NaiKNai3 * NaiKNai3;
        oneNaiKNai3 = 1.0 + Na_i * invKNai;
        oneNaiKNai3 = oneNaiKNai3 * oneNaiKNai3 * oneNaiKNai3;
        oneKiKKi2 = 1.0 + K_i * invKKi;
        oneKiKKi2 = oneKiKKi2 * oneKiKKi2;
        NaoKNao3 = Na_o * invKNao;
        NaoKNao3 = NaoKNao3 * NaoKNao3 * NaoKNao3;
        oneNaoKNao3 = 1.0 + Na_o * invKNao;
        oneNaoKNao3 = oneNaoKNao3 * oneNaoKNao3 * oneNaoKNao3;
        KiKKi2 = K_i * invKKi;
        KiKKi2 = KiKKi2 * KiKKi2;

        a1b4denom = 1.0 / (oneNaiKNai3 + oneKiKKi2 - 1.0);
        b2a3denom = 1.0 / (oneNaoKNao3 + oneKoKKo2 - 1.0);

        alpha1 = (k1pNaK * NaiKNai3) * a1b4denom;
        beta2 = (k2mNaK * NaoKNao3) * b2a3denom;
        alpha3 = (k3pNaK * KoKKo2) * b2a3denom;
        beta3 = (k3mNaK * P * Hp) * b3denom;
        beta4 = (k4mNaK * KiKKi2) * a1b4denom;

        // x1NaK = alpha4 * alpha1 * alpha2 + beta2 * beta4 * beta3 + alpha2 * beta4 * beta3 + beta3 * alpha1 * alpha2
        // Correction from
        // https://journals.plos.org/ploscompbiol/article/comment?id=10.1371/annotation/ec788479-502e-4dbb-8985-52b9af657990
        // (beta2 -> beta1)
        x1NaK = alpha4 * alpha1 * alpha2 + beta1 * beta4 * beta3 + alpha2 * beta4 * beta3 + beta3 * alpha1 * alpha2;
        x2NaK = beta2 * beta1 * beta4 + alpha1 * alpha2 * alpha3 + alpha3 * beta1 * beta4 + alpha2 * alpha3 * beta4;
        x3NaK = alpha2 * alpha3 * alpha4 + beta3 * beta2 * beta1 + beta2 * beta1 * alpha4 + alpha3 * alpha4 * beta1;
        x4NaK = beta4 * beta3 * beta2 + alpha3 * alpha4 * alpha1 + beta2 * alpha4 * alpha1 + beta3 * beta2 * alpha1;

        invnaks = 1.0 / (x1NaK + x2NaK + x3NaK + x4NaK);

        E1NaK = x1NaK * invnaks;
        E2NaK = x2NaK * invnaks;
        E3NaK = x3NaK * invnaks;
        E4NaK = x4NaK * invnaks;

        JNaKNa = 3.0 * (E1NaK * alpha3 - E2NaK * beta3);
        JNaKK = 2.0 * (E4NaK * beta1 - E3NaK * alpha1);
        // Note: GNaKbar is not a parameter in the original model; the value of 30 is used in the equation.
        INaK = GNaKbar * (zNa * JNaKNa + zK * JNaKK);

        //
        // Background currents (INab, ICab, IKb) and sarcolemmal calcium pump current (IpCa)
        //

        INab = PNab * INab_coeff * (Na_i * exp_znavfrt - Na_o);
        ICab = PCab * ICab_coeff * (gammaCai * Ca_i * exp_zcavfrt - gammaCao * Ca_o);

        IKb = GKbbar * xKb * (V - EK);
        IpCa = GpCabar * (Ca_i / (0.0005 + Ca_i));

        //
        // Voltage
        //

        prev_u = V;
        V = V - dt * (INa + Ito + ICaL + ICaNa + ICaK + IKr + IKs + IK1 + INaCa + INaK + INab + ICab + IKb + IpCa + Istim);
        float u = save_v * V;

        //
        // Calcium/calmodulin-dependent protein kinase (CaMK) integration
        //

        CaMKtrap = CaMKtrap + dt * (alphaCaMK * CaMKbound * (CaMKbound + CaMKtrap) - betaCaMK * CaMKtrap);

        //
        // Diffusion fluxes (JdiffNa, JdiffCa, JdiffK)
        //

        JdiffNa = (Na_ss - Na_i) * invtaudiffna;
        JdiffCa = (Ca_ss - Ca_i) * invtaudiffca;
        JdiffK = (K_ss - K_i) * invtaudiffk;

        //
        // SR calcium release flux, via ryanodine receptor (Jrel)
        //

        pow_cajsr = 1.5/Ca_jsr;
        pow_cajsr = pow_cajsr * pow_cajsr;
        pow_cajsr = pow_cajsr * pow_cajsr;
        pow_cajsr = pow_cajsr * pow_cajsr;
        pow_cajsr = 1.0 / (1.0 + pow_cajsr);

        JrelNPinf = alpharel * (-ICaL) * pow_cajsr;
        taurelNP = max(betatau / (1.0 + (0.0123 / Ca_jsr)), 0.001);
        JrelNP = JrelNPinf - (JrelNPinf - JrelNP) * exp(-dt/taurelNP);

        JrelCaMKinf = alpharelCaMK * (-ICaL) * pow_cajsr;
        taurelCaMK = max(betatauCaMK / (1.0 + (0.0123 / Ca_jsr)), 0.001);
        JrelCaMK = JrelCaMKinf - (JrelCaMKinf - JrelCaMK) * exp(-dt/taurelCaMK);

        OrelCaMK = OINaCaMK;

        Jrel = (1.0 - OrelCaMK) * JrelNP + OrelCaMK * JrelCaMK;

        //
        // Calcium uptake via SERCA pump (Jup)
        //

        JupNP = (0.004375 * Ca_i) / (0.00092 + Ca_i);
        JupCaMK = (1.0 + DeltaJupCaMKbar) * ((0.004375 * Ca_i) / (0.00092 - DeltaKmPLBbar + Ca_i));

        OupCaMK = OINaCaMK;

        Jleak = (0.0039375 * Ca_nsr) * inv15;

        Jup = (1.0 - OupCaMK) * JupNP + OupCaMK * JupCaMK - Jleak;

        //
        // Calcium translocation from NSR to JSR (Jtr)
        //

        Jtr = (Ca_nsr - Ca_jsr) * invtautr;

        //
        // Concentrations and buffers
        //

        // Na_i = Na_i + dt * (-(INa + INaL + 3.0 * INaCai + 8.0 * INaK + INab) * (Acap / (FF*vmyo)) + JdiffNa * (vss/vmyo))
        // Correction 1 from
        // https://journals.plos.org/ploscompbiol/article/comment?id=10.1371/annotation/0b8121cd-4280-4ff7-91d9-e9887bcce396
        // (Remove INaL)
        // Na_i = Na_i + dt * (-(INa + 3.0 * INaCai + 3.0 * INaK + INab)*(Acap/(FF*vmyo)) + JdiffNa*(vss/vmyo));
        Na_i_diff = Na_i_diff + dt * (-(INa + 3.0 * INaCai + 3.0 * INaK + INab) * acapofvmyo + JdiffNa * vssovmyo);
        // Na_ss = Na_ss + dt * (-(ICaNa + 3.0 * INaCass) * (Acap / (FF*vss)) - JdiffNa);
        Na_ss_diff = Na_ss_diff + dt * (-(ICaNa + 3.0 * INaCass) * acapofvss - JdiffNa);
        // Correction 2 from
        // https://journals.plos.org/ploscompbiol/article/comment?id=10.1371/annotation/0b8121cd-4280-4ff7-91d9-e9887bcce396
        // (IKur -> IKb)
        // K_i = K_i + dt * (-(Ito + IKr + IKs + IK1 + IKb + Istim - 2.0 * INaK)*(Acap/(FF*vmyo)) + JdiffK*(vss/vmyo));
        K_i_diff = K_i_diff + dt * (-(Ito + IKr + IKs + IK1 + IKb + Istim - 2.0 * INaK) * acapofvmyo + JdiffK * vssovmyo);
        // K_ss = K_ss + dt * (-ICaK * (Acap / (FF*vss)) - JdiffK);
        K_ss_diff = K_ss_diff + dt * (-ICaK * acapofvss - JdiffK);

        bc1 = KmCMDN+Ca_i;
        bc1 = bc1*bc1;
        bc2 = KmTRPN+Ca_i;
        bc2 = bc2*bc2;
        betaCai_inv = 1.0 + cmdnkmcmdn/bc1 + trpnkmtrpn/bc2;
        Ca_i = Ca_i + dt * ((-(IpCa + ICab - 2.0 * INaCai) * acapo2fvmyo - Jup * vnsrovmyo + JdiffCa * vssovmyo) / betaCai_inv);
        u += save_ca * Ca_i;

        bc1 = KmBSR+Ca_ss;
        bc1 = bc1*bc1;
        bc2 = KmBSL+Ca_ss;
        bc2 = bc2*bc2;
        betaCass_inv = 1.0 + bsrkmbsr/bc1 + bslkmbsl/bc2;
        Ca_ss = Ca_ss + dt * ((-(ICaL - 2.0 * INaCass) * acapo2fvss + Jrel * vjsrovss - JdiffCa) / betaCass_inv);

        Ca_nsr = Ca_nsr + dt * (Jup - Jtr * vjsrovnsr);

        bc1 = KmCSQN+Ca_jsr;
        bc1 = bc1*bc1;
        betaCajsr_inv = 1.0 + csqnkmcsqn/bc1;
        Ca_jsr = Ca_jsr + dt * ((Jtr - Jrel) / betaCajsr_inv);

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
