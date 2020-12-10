#version 300 es

precision highp float;
precision highp int;

uniform sampler2D in_particles_1, in_particles_2, in_particles_3, in_particles_4;
uniform sampler2D data_texture;

layout (location = 0) out vec4 error_texture;

in vec2 cc;

// Simulation parameters
uniform float dt, period, stim_start, stim_end, stim_mag;
uniform int num_beats;
uniform float v_init, w_init;

// Macros to get the particles from the textures
#define TR_POS particles_1.r
#define TSI_POS particles_1.g
#define TWP_POS particles_1.b
#define TD_POS particles_1.a

#define TVP_POS particles_2.r
#define TV1M_POS particles_2.g
#define TV2M_POS particles_2.b
#define TWM_POS particles_2.a

#define TO_POS particles_3.r
#define XK_POS particles_3.g
#define UCSI_POS particles_3.b
#define UC_POS particles_3.a

#define UV_POS particles_4.r

// Useful functions
#define TANH(x) ((exp(2.0*(x)) - 1.0) / (exp(2.0*(x)) + 1.0))

void main() {
    // PSO derived parameters
    int num_period = int(ceil(period/dt));
    float endtime = ceil(float(num_beats)*period);
    int num_steps = int(ceil(endtime/dt));

    // Get the relevant color from each texture
    vec4 particles_1 = texture(in_particles_1, cc);
    vec4 particles_2 = texture(in_particles_2, cc);
    vec4 particles_3 = texture(in_particles_3, cc);
    vec4 particles_4 = texture(in_particles_4, cc);

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
        float stim = 0.0;
        if (mod(float(step_count), period/dt) > stim_start/dt && mod(float(step_count), period/dt) < stim_end/dt) {
            stim = stim_mag;
        }

        u -= (jfi+jso+jsi-stim)*dt;

        // Measure error
        if (mod(float(step_count), round(1.0/dt)) == 0.0) {
            float actual = texelFetch(data_texture, ivec2(data_index++, 0), 0).r;
            error += (u - actual)*(u - actual);
        }
    }

    error_texture = vec4(error, 0, 0, 0);
}
