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

#define TR_POS 100.01567
#define TSI_POS 75.586
#define TWP_POS 304.35117
#define TD_POS 0.16317
#define TVP_POS 4.8449
#define TV1M_POS 11.08933
#define TV2M_POS 609.8374633
#define TWM_POS 16.6054
#define TO_POS 7.5254
#define XK_POS 1.65253
#define UCSI_POS 0.55996
#define UC_POS 0.101079
#define UV_POS 0.019021



#define TANH(x) ((exp(2.0*(x)) - 1.0) / (exp(2.0*(x)) + 1.0))
// #define TANH(x) tanh(x)

int main(int argc, char const *argv[])
{
    // if(argc != 2)
    // {
        // printf("Give me an output file name.\n");
    // }

    // FILE *u_out = fopen(argv[1],"w");
    FILE *u_out = fopen("U_out.data","w");

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
    for (int step_count = 0; step_count < num_steps; ++step_count) {
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
            printf("Stimulus!\n");
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
            error+=1;
        }
    }

    // error_texture = vec4(error, 0, 0, 0);
    fclose(u_out);
    printf("%f\n",error);
    return 0;
}