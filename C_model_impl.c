// Raw C impl

#include <stdlib.h>
#include <stdio.h>
#include <math.h>


// #define TR_POS 33.33
// #define TSI_POS 29.0
// #define TWP_POS 870.0
// #define TD_POS 0.25
// #define TVP_POS 3.33
// #define TV1M_POS 19.6
// #define TV2M_POS 1250.0
// #define TWM_POS 41.0
// #define TO_POS 12.5
// #define XK_POS 10.0
// #define UCSI_POS 0.85
// #define UC_POS 0.13
// #define UV_POS 0.04

// #define TR_POS 100.01567
// #define TSI_POS 75.586
// #define TWP_POS 304.35117
// #define TD_POS 0.16317
// #define TVP_POS 4.8449
// #define TV1M_POS 11.08933
// #define TV2M_POS 609.8374633
// #define TWM_POS 16.6054
// #define TO_POS 7.5254
// #define XK_POS 1.65253
// #define UCSI_POS 0.55996
// #define UC_POS 0.101079
// #define UV_POS 0.019021

#define TR_POS 101.4594955444336
#define TSI_POS 98.33734893798828
#define TWP_POS 268.5346984863281
#define TD_POS 0.22014516592025757
#define TVP_POS 6.206625938415527
#define TV1M_POS 47.53730773925781
#define TV2M_POS 1000.5355834960938
#define TWM_POS 91.08232116699219
#define TO_POS 48.045433044433594
#define XK_POS 12.532750129699707
#define UCSI_POS 0.5306475162506104
#define UC_POS 0.13606366515159607
#define UV_POS 0.049632322043180466


#define TANH(x) ((exp(2.0*(x)) - 1.0) / (exp(2.0*(x)) + 1.0))
// #define TANH(x) tanh(x)

float actual_data[400] = {
6.67E-02,
8.86E-02,
1.54E-01,
3.75E-01,
8.26E-01,
9.88E-01,
1.01E+00,
1.00E+00,
1.00E+00,
1.01E+00,
1.01E+00,
1.01E+00,
1.02E+00,
1.02E+00,
1.02E+00,
1.02E+00,
1.03E+00,
1.03E+00,
1.04E+00,
1.03E+00,
1.03E+00,
1.03E+00,
1.03E+00,
1.04E+00,
1.04E+00,
1.04E+00,
1.04E+00,
1.04E+00,
1.04E+00,
1.03E+00,
1.04E+00,
1.03E+00,
1.03E+00,
1.04E+00,
1.04E+00,
1.03E+00,
1.03E+00,
1.03E+00,
1.03E+00,
1.03E+00,
1.02E+00,
1.02E+00,
1.02E+00,
1.02E+00,
1.02E+00,
1.02E+00,
1.01E+00,
1.01E+00,
1.01E+00,
1.01E+00,
1.01E+00,
1.01E+00,
1.00E+00,
1.00E+00,
9.99E-01,
9.99E-01,
9.97E-01,
9.98E-01,
9.89E-01,
9.91E-01,
9.88E-01,
9.88E-01,
9.78E-01,
9.83E-01,
9.81E-01,
9.79E-01,
9.79E-01,
9.74E-01,
9.75E-01,
9.71E-01,
9.61E-01,
9.67E-01,
9.64E-01,
9.53E-01,
9.61E-01,
9.59E-01,
9.57E-01,
9.50E-01,
9.53E-01,
9.51E-01,
9.49E-01,
9.44E-01,
9.31E-01,
9.39E-01,
9.35E-01,
9.33E-01,
9.29E-01,
9.34E-01,
9.27E-01,
9.20E-01,
9.19E-01,
9.16E-01,
9.16E-01,
9.14E-01,
9.13E-01,
9.11E-01,
9.04E-01,
9.03E-01,
9.02E-01,
9.01E-01,
8.98E-01,
8.89E-01,
8.92E-01,
8.89E-01,
8.86E-01,
8.85E-01,
8.81E-01,
8.74E-01,
8.68E-01,
8.68E-01,
8.70E-01,
8.66E-01,
8.62E-01,
8.63E-01,
8.56E-01,
8.56E-01,
8.48E-01,
8.47E-01,
8.48E-01,
8.37E-01,
8.32E-01,
8.34E-01,
8.33E-01,
8.32E-01,
8.26E-01,
8.24E-01,
8.20E-01,
8.12E-01,
8.16E-01,
8.12E-01,
8.05E-01,
8.03E-01,
7.98E-01,
7.97E-01,
7.94E-01,
7.89E-01,
7.94E-01,
7.78E-01,
7.75E-01,
7.72E-01,
7.72E-01,
7.66E-01,
7.61E-01,
7.59E-01,
7.54E-01,
7.51E-01,
7.49E-01,
7.48E-01,
7.38E-01,
7.41E-01,
7.30E-01,
7.24E-01,
7.24E-01,
7.11E-01,
7.09E-01,
7.05E-01,
6.96E-01,
6.95E-01,
6.90E-01,
6.83E-01,
6.84E-01,
6.72E-01,
6.68E-01,
6.65E-01,
6.62E-01,
6.60E-01,
6.52E-01,
6.50E-01,
6.39E-01,
6.30E-01,
6.27E-01,
6.23E-01,
6.16E-01,
6.10E-01,
6.02E-01,
5.99E-01,
5.89E-01,
5.81E-01,
5.79E-01,
5.69E-01,
5.68E-01,
5.60E-01,
5.60E-01,
5.48E-01,
5.40E-01,
5.28E-01,
5.20E-01,
5.13E-01,
5.10E-01,
4.99E-01,
4.86E-01,
4.84E-01,
4.73E-01,
4.68E-01,
4.63E-01,
4.52E-01,
4.49E-01,
4.37E-01,
4.25E-01,
4.07E-01,
4.09E-01,
3.96E-01,
3.90E-01,
3.76E-01,
3.66E-01,
3.55E-01,
3.46E-01,
3.34E-01,
3.22E-01,
3.08E-01,
2.99E-01,
2.91E-01,
2.79E-01,
2.69E-01,
2.56E-01,
2.43E-01,
2.37E-01,
2.26E-01,
2.12E-01,
2.03E-01,
1.91E-01,
1.88E-01,
1.76E-01,
1.65E-01,
1.62E-01,
1.51E-01,
1.44E-01,
1.36E-01,
1.30E-01,
1.27E-01,
1.21E-01,
1.17E-01,
1.10E-01,
1.03E-01,
9.89E-02,
8.97E-02,
9.40E-02,
8.49E-02,
8.17E-02,
8.00E-02,
7.88E-02,
7.40E-02,
7.55E-02,
6.09E-02,
7.33E-02,
6.69E-02,
6.46E-02,
6.83E-02,
6.23E-02,
6.00E-02,
5.79E-02,
5.42E-02,
5.52E-02,
6.04E-02,
4.77E-02,
5.50E-02,
4.71E-02,
4.51E-02,
4.34E-02,
4.16E-02,
4.33E-02,
4.30E-02,
5.31E-02,
3.82E-02,
4.06E-02,
3.89E-02,
3.94E-02,
4.29E-02,
3.86E-02,
3.62E-02,
3.59E-02,
3.19E-02,
3.40E-02,
3.44E-02,
3.22E-02,
3.02E-02,
3.13E-02,
2.91E-02,
2.83E-02,
3.02E-02,
3.09E-02,
2.07E-02,
2.88E-02,
2.58E-02,
2.47E-02,
2.50E-02,
2.57E-02,
2.20E-02,
2.05E-02,
2.05E-02,
2.46E-02,
2.35E-02,
2.43E-02,
2.22E-02,
2.00E-02,
2.00E-02,
2.15E-02,
2.33E-02,
2.89E-02,
2.50E-02,
2.27E-02,
2.24E-02,
2.47E-02,
2.33E-02,
2.49E-02,
2.12E-02,
2.41E-02,
1.80E-02,
2.24E-02,
2.24E-02,
2.15E-02,
2.20E-02,
1.91E-02,
2.04E-02,
2.07E-02,
1.96E-02,
3.49E-02,
2.01E-02,
1.71E-02,
2.09E-02,
2.07E-02,
1.57E-02,
1.71E-02,
1.72E-02,
2.15E-02,
2.06E-02,
1.91E-02,
1.77E-02,
1.67E-02,
2.09E-02,
1.85E-02,
1.44E-02,
2.11E-02,
1.86E-02,
1.47E-02,
1.51E-02,
1.58E-02,
1.81E-02,
1.63E-02,
2.07E-02,
1.75E-02,
1.22E-02,
1.38E-02,
2.04E-02,
1.79E-02,
1.43E-02,
1.90E-02,
2.00E-02,
1.78E-02,
1.72E-02,
1.80E-02,
1.53E-02,
1.45E-02,
1.27E-02,
1.38E-02,
1.57E-02,
1.58E-02,
1.46E-02,
1.59E-02,
1.19E-02,
1.88E-02,
1.61E-02,
2.04E-02,
1.92E-02,
1.56E-02,
8.85E-03,
1.44E-02,
1.41E-02,
1.54E-02,
1.39E-02,
1.19E-02,
1.07E-02,
6.20E-03,
1.58E-02,
1.10E-02,
1.50E-02,
1.18E-02,
1.35E-02,
1.04E-02,
1.33E-02,
1.36E-02,
1.24E-02,
1.51E-02,
1.38E-02,
1.07E-02,
0,
0,
0,
0,
0,
0,
3.05E-02,
4.00E-02,
4.27E-02,
3.94E-02,
3.68E-02,
3.89E-02,
3.47E-02,
4.59E-02,
5.57E-02
};


int main(int argc, char const *argv[])
{
    // if(argc != 2)
    // {
        // printf("Give me an output file name.\n");
    // }

    double act_max = -10;
    for (int i = 0; i < 400; ++i)
    {
        if(actual_data[i] > act_max)
        {
            act_max = actual_data[i];
        }
    }

    for (int i = 0; i < 400; ++i)
    {
        actual_data[i] /= act_max;
    }

    // FILE *u_out = fopen(argv[1],"w");
    FILE *u_out = fopen("U_out.data","w");
    FILE *real_out = fopen("zebrafish.data","w");

    float dt, period, stim_start, stim_end, stim_mag;
    int num_beats;
    float v_init, w_init;

    dt = 0.02;
    period = 400.0;
    stim_start = 2.0;
    stim_end = 7.0;
    stim_mag = 0.1;
    // stim_mag = 1.0;
    num_beats = 1;
    v_init = 1.0;
    w_init = 1.0;



      // PSO derived parameters
    int num_period = (int)(ceil(period/dt));
    float endtime = ceil((float)(num_beats)*period);
    int num_steps = (int)(ceil(endtime/dt));

    // // Get the relevant color from each texture
    // vec4 particles_1 = texture(in_particles_1, cc);
    // vec4 particles_2 = texture(in_particles_2, cc);
    // vec4 particles_3 = texture(in_particles_3, cc);
    // vec4 particles_4 = texture(in_particles_4, cc);

    // Initialize values for the simulation
    float u = 0.0;
    float v = v_init;
    float w = w_init;

    float error = 0.0;

    int data_index = 0;

    // Run the simulation with the current swarm parameters
    for (int step_count = 1; step_count <= num_steps; ++step_count) {
        float p = u >= UC_POS ? 1.0 : 0.0;
        float q = u >= UV_POS ? 1.0 : 0.0;

        float dv =
            // p false
            (1.0 - p)*(1.0 - v) / (
                // q false
                (1.0 - q)*TV1M_POS
                // q true
                + q*TV2M_POS
                )
            // p true
            - p*(v/TVP_POS);

        float dw =
            // p false
            (1.0-p)*(1.0-w) / TWM_POS
            // p true
            - p*(w/TWP_POS);

        v += dt*dv;
        w += dt*dw;

        float jfi =
            // p true
            p * -v * (u - UC_POS)*(1.0-u)/TD_POS;

        float jso =
            // p false
            (1.0 - p)*u / TO_POS
            // p true
            + p/TR_POS;

        float jsi = -w * (1.0 + TANH(XK_POS * (u-UCSI_POS))) / (2.0 * TSI_POS);

        // Apply stimulus
        // 
        float stim = 0.0;
        int stimCheck = step_count % (int)(period/dt);
        
        // if (mod(float(step_count), period/dt) > stim_start/dt && mod(float(step_count), period/dt) < stim_end/dt) {
        if(stimCheck > stim_start/dt && stimCheck < stim_end/dt)
        {
            // printf("Stimulus!\n");
            stim = stim_mag;
        }

        u -= (jfi+jso+jsi-stim)*dt;

        // Measure error
        // if (mod(float(step_count), round(1.0/dt)) == 0.0)
        if(step_count % (int)(1/dt) == 0)
        {
            // float actual = texelFetch(data_texture, ivec2(data_index++, 0), 0).r;
            // error += (u - actual)*(u - actual);
            fprintf(u_out, "%f\n", u);

            // printf("%f\n",error);
            float actual = actual_data[data_index++];
            fprintf(real_out, "%f\n", actual);
            error+= (u-actual) * (u-actual);
        }
    }

    // error_texture = vec4(error, 0, 0, 0);
    fclose(u_out);
    fclose(real_out);
    printf("%f\n",error);
    return 0;
}