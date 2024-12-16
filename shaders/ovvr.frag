#version 300 es

precision highp float;
precision highp int;

uniform sampler2D in_particles_1, in_particles_2, data_texture;

layout (location = 0) out vec4 error_texture;

in vec2 cc;

uniform float dt, period;
uniform int num_beats, pre_beats, data_type;
uniform float align_thresh;
uniform float sample_interval, apd_thresh, weight;
uniform float stim_dur, stim_mag, stim_offset_1, stim_offset_2, stim_t_scale;
uniform bool stim_biphasic;

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
    float minf, taum, hinf, tauhfast, tauhslow, h;
    float jinf, tauj, hCaMKinf, tauhCaMKslow, hCaMKfast, hCaMK, jCaMKinf, taujCaMK, OINaCaMK, INafast, mLinf, taumL, hLinf;
    float hLCaMKinf, OINaLCaMK, INalate;
    float ainf, taua, iinf, tauifast, tauislow, Aifast, Aislow, i, aCaMKinf, tauaCaMK, iCaMKinf, deltaCaMKdevelop;
    float deltaCaMKrecover, tauiCaMKfast, tauiCaMKslow, AiCaMKfast, AiCaMKslow, iCaMK, OItoCaMK;
    float dinf, taud, finf, tauffast, taufslow, f, fCainf, taufCafast, taufCaslow, AfCafast, AfCaslow, fCa, jCainf, fCaMKinf;
    float taufCaMKfast, fCaMKslow, fCaMK, fCaCaMKinf, taufCaCaMKfast, fCaCaMKslow, fCaCaMK, km2n, alphan, PsiCa;
    float ICaLbar, PsiCaNa, ICaNabar, PsiCaK, ICaKbar, ICaLCaMKbar, ICaNaCaMKbar, ICaKCaMKbar, OICaLCaMK;
    float xrinf, tauxrfast, tauxrslow, Axrfast, Axrslow, xr, RKr;
    float xs1inf, tauxs1, xs2inf, tauxs2;
    float xK1inf, tauxK1, RK1;
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
    float betaCai, betaCass, betaCajsr;

    V = -87.84;
    Na_i = 7.23;
    Na_ss = 7.23;
    K_i = 143.79;
    K_ss = 143.79;
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

        minf = 1.0 / (1.0 + exp(-(V + 39.57)/9.871));
        taum = 1.0 / (6.765 * exp((V+11.64)/34.77) + 8.552 * exp(-(V+77.42)/5.955));
        m = minf - (minf - m) * exp(-dt/taum);

        hinf = 1.0 / (1.0 + exp((V + 82.9)/6.086));
        tauhfast = 1.0 / (1.432e-5 * exp(-(V+1.196)/6.285) + 6.149 * exp((V+0.5096)/20.27));
        tauhslow = 1.0 / (0.009764 * exp(-(V+17.95)/28.05) + 0.3343 * exp((V+5.730)/56.66));
        hfast = hinf - (hinf - hfast) * exp(-dt/tauhfast);
        hslow = hinf - (hinf - hslow) * exp(-dt/tauhslow);
        h = Ahfast * hfast + Ahslow * hslow;

        jinf = hinf;
        tauj = 2.038 + 1.0 / (0.02136 * exp(-(V+100.6)/8.281) + 0.3052 * exp((V+0.9941)/38.45));
        j = jinf - (jinf - j) * exp(-dt/tauj);

        hCaMKinf = 1.0 / (1.0 + exp((V+89.1)/6.086));
        tauhCaMKslow = 3.0 * tauhslow;
        hCaMKfast = hfast;
        hCaMKslow = hCaMKinf - (hCaMKinf - hCaMKslow) * exp(-dt/tauhCaMKslow);
        hCaMK = AhCaMKfast * hCaMKfast + AhCaMKslow * hCaMKslow;

        jCaMKinf = jinf;
        taujCaMK = 1.46 * tauj;
        jCaMK = jCaMKinf - (jCaMKinf - jCaMK) * exp(-dt/taujCaMK);

        OINaCaMK = 1.0 / (1.0 + KmCaMK/CaMKactive);

        INafast = GNafastbar * (V - ENa) * m*m*m * ((1.0 - OINaCaMK) * h * j + OINaCaMK * hCaMK * jCaMK);

        mLinf = 1.0 / (1.0 + exp(-(V+42.85)/5.264));
        taumL = taum;
        mL = mLinf - (mLinf - mL) * exp(-dt/taumL);

        hLinf = 1.0 / (1.0 + exp((V+87.61)/7.488));
        hL = hLinf - (hLinf - hL) * exp(-dt/tauhL);

        hLCaMKinf = 1.0 / (1.0 + exp((V+93.81)/7.488));
        hLCaMK = hLCaMKinf - (hLCaMKinf - hLCaMK) * exp(-dt/tauhLCaMK);

        OINaLCaMK = OINaCaMK;

        INalate = GNalatebar * (V - ENa) * mL * ((1.0 - OINaLCaMK) * hL + OINaLCaMK * hLCaMK);

        INa = INafast + INalate;

        //
        // Transient outward potassium current (Ito)
        //

        ainf = 1.0 / (1.0 + exp(-(V-14.34)/14.82));
        taua = 1.0515 / ((1.0 / (1.2089 * (1.0 + exp(-(V-18.41)/29.38)))) + (3.5 / (1.0 + exp((V+100.0)/29.38))));
        a = ainf - (ainf - a) * exp(-dt/taua);

        iinf = 1.0 / (1.0 + exp((V+43.94)/5.711));
        tauifast = 4.562 + 1.0 / (0.3933 * exp(-(V+100.0)/100.0) + 0.08004 * exp((V+50.0)/16.59));
        tauislow = 23.62 + 1.0 / (0.001416 * exp(-(V+96.52)/59.05) + 1.7808e-8 * exp((V+114.1)/8.079));
        ifast = iinf - (iinf - ifast) * exp(-dt/tauifast);
        islow = iinf - (iinf - islow) * exp(-dt/tauislow);
        Aifast = 1.0 / (1.0 + exp((V-213.6)/151.2));
        Aislow = 1.0 - Aifast;
        i = Aifast * ifast + Aislow * islow;

        aCaMKinf = 1.0 / (1.0 + exp(-(V-24.34)/14.82));
        tauaCaMK = taua;
        aCaMK = aCaMKinf - (aCaMKinf - aCaMK) * exp(-dt/tauaCaMK);

        iCaMKinf = iinf;
        deltaCaMKdevelop = 1.354 + 1e-4 / (exp((V-167.4)/15.89) + exp(-(V-12.23)/0.2154));
        deltaCaMKrecover = 1.0 - 0.5 / (1.0 + exp((V+70.0)/20.0));
        tauiCaMKfast = tauifast * deltaCaMKdevelop * deltaCaMKrecover;
        tauiCaMKslow = tauislow * deltaCaMKdevelop * deltaCaMKrecover;
        iCaMKfast = iCaMKinf - (iCaMKinf - iCaMKfast) * exp(-dt/tauiCaMKfast);
        iCaMKslow = iCaMKinf - (iCaMKinf - iCaMKslow) * exp(-dt/tauiCaMKslow);
        AiCaMKfast = Aifast;
        AiCaMKslow = Aislow;
        iCaMK = AiCaMKfast * iCaMKfast + AiCaMKslow * iCaMKslow;

        OItoCaMK = OINaCaMK;

        Ito = Gtobar * (V - EK) * ((1.0 - OItoCaMK) * a * i + OItoCaMK * aCaMK * iCaMK);

        //
        // L-type calcium current (ICaL)
        //

        dinf = 1.0 / (1.0 + exp(-(V+3.940)/4.230));
        taud = 0.6 + 1.0 / (exp(-0.05 * (V+6.0)) + exp(0.09 * (V+14.0)));
        d = dinf - (dinf - d) * exp(-dt/taud);

        finf = 1.0 / (1.0 + exp((V+19.58)/3.696));
        // Note: This equation seems to be implemented correctly, but does not match the plot in Figure
        // 1C of the paper.
        tauffast = 7.0 + 1.0 / (0.0045 * exp(-(V+20.0)/10.0) + 0.0045 * exp((V+20.0)/10.0));
        taufslow = 1000.0 + 1.0 / (0.000035 * exp(-(V+5.0)/4.0) + 0.000035 * exp((V+5.0)/6.0));
        ffast = finf - (finf - ffast) * exp(-dt/tauffast);
        fslow = finf - (finf - fslow) * exp(-dt/taufslow);
        f = Affast * ffast + Afslow * fslow;

        fCainf = finf;
        taufCafast = 7.0 + 1.0 / (0.04 * exp(-(V-4.0)/7.0) + 0.04 * exp((V-4.0)/7.0));
        taufCaslow = 100.0 + 1.0 / (0.00012 * exp(-V/3.0) + 0.00012 * exp(V/7.0));
        fCafast = fCainf - (fCainf - fCafast) * exp(-dt/taufCafast);
        fCaslow = fCainf - (fCainf - fCaslow) * exp(-dt/taufCaslow);
        AfCafast = 0.3 + 0.6 / (1.0 + exp((V-10.0)/10.0));
        AfCaslow = 1.0 - AfCafast;
        fCa = AfCafast * fCafast + AfCaslow * fCaslow;

        jCainf = fCainf;
        jCa = jCainf - (jCainf - jCa) * exp(-dt/taujCa);

        fCaMKinf = finf;
        taufCaMKfast = 2.5 * tauffast;
        fCaMKfast = fCaMKinf - (fCaMKinf - fCaMKfast) * exp(-dt/taufCaMKfast);
        fCaMKslow = fslow;
        fCaMK = AfCaMKfast * fCaMKfast + AfCaMKslow * fCaMKslow;

        fCaCaMKinf = finf;
        taufCaCaMKfast = 2.5 * taufCafast;
        fCaCaMKslow = fCaslow;
        fCaCaMKfast = fCaCaMKinf - (fCaCaMKinf - fCaCaMKfast) * exp(-dt/taufCaCaMKfast);
        fCaCaMK = AfCaCaMKfast * fCaCaMKfast + AfCaCaMKslow * fCaCaMKslow;

        km2n = jCa;
        alphan = 1.0 / (kp2n/km2n + pow((1.0 + Kmn/Ca_ss), 4.0));
        n = alphan * (kp2n/km2n) - (alphan * (kp2n/km2n) - n) * exp(-km2n*dt);

        PsiCa = zCa*zCa * V*FFRT * (gammaCai * Ca_ss * exp(zCa*V*FRT) - gammaCao * Ca_o) / (exp(zCa*V*FRT) - 1.0);
        ICaLbar = PCa * PsiCa;

        PsiCaNa = zNa*zNa * V*FFRT * (gammaNai * Na_ss * exp(zNa*V*FRT) - gammaNao * Na_o) / (exp(zNa*V*FRT) - 1.0);
        ICaNabar = PCaNa * PsiCaNa;

        PsiCaK = zK*zK * V*FFRT * (gammaKi * K_ss * exp(zK*V*FRT) - gammaKo * K_o) / (exp(zK*V*FRT) - 1.0);
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

        xrinf = 1.0 / (1.0 + exp(-(V+8.337)/6.789));
        tauxrfast = 12.98 + 1.0 / (0.3652 * exp((V-31.66)/3.869) + 4.123e-5 * exp(-(V-47.78)/20.38));
        tauxrslow = 1.865 + 1.0 / (0.06629 * exp((V-34.70)/7.355) + 1.128e-5 * exp(-(V-29.74)/25.94));
        xrfast = xrinf - (xrinf - xrfast) * exp(-dt/tauxrfast);
        xrslow = xrinf - (xrinf - xrslow) * exp(-dt/tauxrslow);
        Axrfast = 1.0 / (1.0 + exp((V+54.81)/38.21));
        Axrslow = 1.0 - Axrfast;
        xr = Axrfast * xrfast + Axrslow * xrslow;

        RKr = 1.0 / ((1.0 + exp((V+55.0)/75.0)) * (1.0 + exp((V-10.0)/30.0)));

        IKr = GKrbar * sqrt(K_o/5.4) * xr * RKr * (V - EK);

        //
        // Slow delayed rectifier potassium current (IKs)
        //

        // Note: The tail below -15mV looks slightly different than Figure 3D from the paper.
        xs1inf = 1.0 / (1.0 + exp(-(V+11.60)/8.932));
        tauxs1 = 817.3 + 1.0 / (2.326e-4 * exp((V+48.28)/17.80) + 0.001292 * exp(-(V+210.0)/230.0));
        xs1 = xs1inf - (xs1inf - xs1) * exp(-dt/tauxs1);

        xs2inf = xs1inf;
        tauxs2 = 1.0 / (0.01 * exp((V-50.0)/20.0) + 0.0193 * exp(-(V+66.54)/31.0));
        xs2 = xs2inf - (xs2inf - xs2) * exp(-dt/tauxs2);

        IKs = GKsbar * (1.0 + 0.6 / (1.0 + pow((3.8e-5/Ca_i), 1.4))) * xs1 * xs2 * (V - EKs);

        //
        // Inward rectifier potassium current (IK1)
        //

        xK1inf = 1.0 / (1.0 + exp(-(V+2.5538*K_o+144.59)/(1.5692*K_o+3.8115)));
        // Note: The plot in Figure 2C of the paper has larger values, but the same shape, than found
        // here.
        tauxK1 = 122.2 / (exp(-(V+127.2)/20.36) + exp((V+236.8)/69.33));
        xK1 = xK1inf - (xK1inf - xK1) * exp(-dt/tauxK1);

        RK1 = 1.0 / (1.0 + exp((V+105.8-2.6*K_o)/9.493));

        IK1 = GK1bar * sqrt(K_o) * xK1 * RK1 * (V - EK);

        //
        // Sodium/calcium exchange current (INaCa)
        //

        hCa = exp(qCa*V*FRT);
        hNa = exp(qNa*V*FRT);

        h1i = 1.0 + (Na_i/kNa3) * (1.0 + hNa);
        h1ss = 1.0 + (Na_ss/kNa3) * (1.0 + hNa);
        h2i = (Na_i*hNa)/(kNa3*h1i);
        h2ss = (Na_ss*hNa)/(kNa3*h1ss);
        h3i = 1.0/h1i;
        h3ss = 1.0/h1ss;
        h4i = 1.0 + (Na_i/kNa1) * (1.0 + Na_i/kNa2);
        h4ss = 1.0 + (Na_ss/kNa1) * (1.0 + Na_ss/kNa2);
        h5i = (Na_i*Na_i) / (h4i * kNa1 * kNa2);
        h5ss = (Na_ss*Na_ss) / (h4ss * kNa1 * kNa2);
        h6i = 1.0 / h4i;
        h6ss = 1.0 / h4ss;
        h7 = 1.0 + (Na_o/kNa3) * (1.0 + 1.0/hNa);
        h8 = Na_o / (kNa3 * hNa * h7);
        h9 = 1.0 / h7;

        k1 = h12 * Ca_o * kCaon;
        k3p = h9 * omegaCa;
        k3pp = h8 * omegaNaCa;
        k3 = k3p + k3pp;
        k4pi = (h3i * omegaCa) / hCa;
        k4pss = (h3ss * omegaCa) / hCa;
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

        E1i = x1i / (x1i + x2i + x3i + x4i);
        E1ss = x1ss / (x1ss + x2ss + x3ss + x4ss);
        E2i = x2i / (x1i + x2i + x3i + x4i);
        E2ss = x2ss / (x1ss + x2ss + x3ss + x4ss);
        E3i = x3i / (x1i + x2i + x3i + x4i);
        E3ss = x3ss / (x1ss + x2ss + x3ss + x4ss);
        E4i = x4i / (x1i + x2i + x3i + x4i);
        E4ss = x4ss / (x1ss + x2ss + x3ss + x4ss);

        alloi = 1.0 / (1.0 + pow((KmCaAct/Ca_i), 2.0));
        alloss = 1.0 / (1.0 + pow((KmCaAct/Ca_ss), 2.0));

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

        KNai = KNaio * exp((Delta*V*FRT)/3.0);
        KNao = KNaoo * exp(((1.0-Delta)*V*FRT)/3.0);

        P = SigmaP / (1.0 + Hp/KHP + Na_i/KNaP + K_i/KKP);

        NaiKNai3 = pow((Na_i/KNai), 3.0);
        oneNaiKNai3 = pow((1.0 + Na_i/KNai), 3.0);
        oneKiKKi2 = pow((1.0 + K_i/KKi), 2.0);
        NaoKNao3 = pow((Na_o/KNao), 3.0);
        oneNaoKNao3 = pow((1.0 + Na_o/KNao), 3.0);
        KiKKi2 = pow((K_i/KKi), 2.0);

        alpha1 = (k1pNaK * NaiKNai3) / (oneNaiKNai3 + oneKiKKi2 - 1.0);
        beta2 = (k2mNaK * NaoKNao3) / (oneNaoKNao3 + oneKoKKo2 - 1.0);
        alpha3 = (k3pNaK * KoKKo2) / (oneNaoKNao3 + oneKoKKo2 - 1.0);
        beta3 = (k3mNaK * P * Hp) / (1.0 + MgATP/KMgATP);
        beta4 = (k4mNaK * KiKKi2) / (oneNaiKNai3 + oneKiKKi2 - 1.0);

        // x1NaK = alpha4 * alpha1 * alpha2 + beta2 * beta4 * beta3 + alpha2 * beta4 * beta3 + beta3 * alpha1 * alpha2
        // Correction from
        // https://journals.plos.org/ploscompbiol/article/comment?id=10.1371/annotation/ec788479-502e-4dbb-8985-52b9af657990
        // (beta2 -> beta1)
        x1NaK = alpha4 * alpha1 * alpha2 + beta1 * beta4 * beta3 + alpha2 * beta4 * beta3 + beta3 * alpha1 * alpha2;
        x2NaK = beta2 * beta1 * beta4 + alpha1 * alpha2 * alpha3 + alpha3 * beta1 * beta4 + alpha2 * alpha3 * beta4;
        x3NaK = alpha2 * alpha3 * alpha4 + beta3 * beta2 * beta1 + beta2 * beta1 * alpha4 + alpha3 * alpha4 * beta1;
        x4NaK = beta4 * beta3 * beta2 + alpha3 * alpha4 * alpha1 + beta2 * alpha4 * alpha1 + beta3 * beta2 * alpha1;

        E1NaK = x1NaK / (x1NaK + x2NaK + x3NaK + x4NaK);
        E2NaK = x2NaK / (x1NaK + x2NaK + x3NaK + x4NaK);
        E3NaK = x3NaK / (x1NaK + x2NaK + x3NaK + x4NaK);
        E4NaK = x4NaK / (x1NaK + x2NaK + x3NaK + x4NaK);

        JNaKNa = 3.0 * (E1NaK * alpha3 - E2NaK * beta3);
        JNaKK = 2.0 * (E4NaK * beta1 - E3NaK * alpha1);
        // Note: GNaKbar is not a parameter in the original model; the value of 30 is used in the equation.
        INaK = GNaKbar * (zNa * JNaKNa + zK * JNaKK);

        //
        // Background currents (INab, ICab, IKb) and sarcolemmal calcium pump current (IpCa)
        //

        INab = PNab * zNa * zNa * V * FFRT * (Na_i * exp(zNa*V*FRT) - Na_o) / (exp(zNa*V*FRT) - 1.0);
        ICab = PCab * zCa * zCa * V * FFRT * (gammaCai * Ca_i * exp(zCa*V*FRT) - gammaCao * Ca_o) / (exp(zCa*V*FRT) - 1.0);

        xKb = 1.0 / (1.0 + exp(-(V-14.48)/18.34));

        IKb = GKbbar * xKb * (V - EK);
        IpCa = GpCabar * (Ca_i / (0.0005 + Ca_i));

        //
        // Voltage
        //

        float prev_u = V;
        V = V - dt * (INa + Ito + ICaL + ICaNa + ICaK + IKr + IKs + IK1 + INaCa + INaK + INab + ICab + IKb + IpCa + Istim) / C_M;
        float u = V;

        //
        // Calcium/calmodulin-dependent protein kinase (CaMK) integration
        //

        CaMKtrap = CaMKtrap + dt * (alphaCaMK * CaMKbound * (CaMKbound + CaMKtrap) - betaCaMK * CaMKtrap);

        //
        // Diffusion fluxes (JdiffNa, JdiffCa, JdiffK)
        //

        JdiffNa = (Na_ss - Na_i) / taudiffNa;
        JdiffCa = (Ca_ss - Ca_i) / taudiffCa;
        JdiffK = (K_ss - K_i) / taudiffK;

        //
        // SR calcium release flux, via ryanodine receptor (Jrel)
        //

        JrelNPinf = (alpharel * (-ICaL)) / (1.0 + pow((1.5/Ca_jsr), 8.0));
        taurelNP = max(betatau / (1.0 + (0.0123 / Ca_jsr)), 0.001);
        JrelNP = JrelNPinf - (JrelNPinf - JrelNP) * exp(-dt/taurelNP);

        JrelCaMKinf = (alpharelCaMK * (-ICaL)) / (1.0 + pow((1.5/Ca_jsr), 8.0));
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

        Jleak = (0.0039375 * Ca_nsr) / 15.0;

        Jup = (1.0 - OupCaMK) * JupNP + OupCaMK * JupCaMK - Jleak;

        //
        // Calcium translocation from NSR to JSR (Jtr)
        //

        Jtr = (Ca_nsr - Ca_jsr) / tautr;

        //
        // Concentrations and buffers
        //

        // Na_i = Na_i + dt * (-(INa + INaL + 3.0 * INaCai + 8.0 * INaK + INab) * (Acap / (FF*vmyo)) + JdiffNa * (vss/vmyo))
        // Correction 1 from
        // https://journals.plos.org/ploscompbiol/article/comment?id=10.1371/annotation/0b8121cd-4280-4ff7-91d9-e9887bcce396
        // (Remove INaL)
        Na_i = Na_i + dt * (-(INa + 3.0 * INaCai + 3.0 * INaK + INab)*(Acap/(FF*vmyo)) + JdiffNa*(vss/vmyo));
        Na_ss = Na_ss + dt * (-(ICaNa + 3.0 * INaCass) * (Acap / (FF*vss)) - JdiffNa);
        // Correction 2 from
        // https://journals.plos.org/ploscompbiol/article/comment?id=10.1371/annotation/0b8121cd-4280-4ff7-91d9-e9887bcce396
        // (IKur -> IKb)
        K_i = K_i + dt * (-(Ito + IKr + IKs + IK1 + IKb + Istim - 2.0 * INaK)*(Acap/(FF*vmyo)) + JdiffK*(vss/vmyo));
        K_ss = K_ss + dt * (-ICaK * (Acap / (FF*vss)) - JdiffK);

        betaCai = 1.0 / (1.0 + (CMDN*KmCMDN)/(pow((KmCMDN+Ca_i), 2.0)) + (TRPN*KmTRPN)/(pow((KmTRPN+Ca_i), 2.0)));
        Ca_i = Ca_i + dt * (betaCai * (-(IpCa + ICab - 2.0 * INaCai) * (Acap/(2.0*FF*vmyo)) - Jup * (vnsr/vmyo) + JdiffCa * (vss/vmyo)));

        betaCass = 1.0 / (1.0 + (BSR*KmBSR)/(pow((KmBSR+Ca_ss), 2.0)) + (BSL*KmBSL)/(pow((KmBSL+Ca_ss), 2.0)));
        Ca_ss = Ca_ss + dt * (betaCass * (-(ICaL - 2.0 * INaCass) * (Acap/(2.0*FF*vss)) + Jrel * (vjsr/vss) - JdiffCa));

        Ca_nsr = Ca_nsr + dt * (Jup - Jtr * (vjsr/vnsr));

        betaCajsr = 1.0 / (1.0 + (CSQN*KmCSQN)/(pow((KmCSQN+Ca_jsr), 2.0)));
        Ca_jsr = Ca_jsr + dt * (betaCajsr * (Jtr - Jrel));

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
