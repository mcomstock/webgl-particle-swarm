#version 300 es

precision highp float;
precision highp int;

uniform sampler2D in_particles_1;
uniform sampler2D data_texture;

layout (location = 0) out vec4 error_texture;

in vec2 cc;

// Simulation parameters
uniform float dt, period, stim_start, stim_end, stim_mag;
uniform int num_beats, pre_beats;
uniform float v_init, w_init;
uniform float align_thresh;
uniform float sample_interval;

float stim_f(const float t) {
    const float stim_scale = 0.4;
    const float stim_dur = 10.0;
    const float offset_1 = 7.0;
    const float offset_2 = offset_1 * 0.96;
    const float t_scale = 0.725;

    return ( -stim_scale * ( t / t_scale - offset_1) / (1.0 + pow((t/t_scale - offset_2) , 4.0) ) );
}

void main() {
    // PSO derived parameters
    int num_period = int(ceil(period/dt));
    int total_beats = pre_beats + num_beats;
    float endtime = ceil(float(total_beats)*period);
    float pre_pace_endtime = ceil(float(pre_beats)*period);
    int pre_pace_steps = int(ceil(pre_pace_endtime/dt));
    int num_steps = int(ceil(endtime/dt));

    const float stim_dur = 10.0;

    ivec2 tex_size = textureSize(in_particles_1, 0);
    ivec2 idx = ivec2(floor(cc * 0.5 * vec2(tex_size)));

    vec4 particles_1 = texelFetch(in_particles_1, idx, 0);
    vec4 particles_2 = texelFetch(in_particles_1, idx + ivec2(tex_size.x/2, 0), 0);
    vec4 particles_3 = texelFetch(in_particles_1, idx + ivec2(0, tex_size.y/2), 0);
    vec4 particles_4 = texelFetch(in_particles_1, idx + ivec2(tex_size.x/2, tex_size.y/2), 0);

    float tr = particles_1.r;
    float tsi = particles_1.g;
    float twp = particles_1.b;
    float td = particles_1.a;

    float tvp = particles_2.r;
    float tv1m = particles_2.g;
    float tv2m = particles_2.b;
    float twm = particles_2.a;

    float to = particles_3.r;
    float xk = particles_3.g;
    float ucsi = particles_3.b;
    float uc = particles_3.a;

    float uv = particles_4.r;

    // Initialize values for the simulation
    float u = 0.0;
    float v = v_init;
    float w = w_init;

    float compare_stride = round(sample_interval / dt);

    float error = 10000000000.0;

    int data_index = 0;

    int start_comp = 0;
    bool first_upstroke = false;
    float saved_value = -1.0;

    // Run the simulation with the current swarm parameters
    for (int step_count = 1; step_count <= num_steps; ++step_count) {
        float p = u >= uc ? 1.0 : 0.0;
        float q = u >= uv ? 1.0 : 0.0;

        float dv =
            // p false
            (1.0 - p)*(1.0 - v) / (
                // q false
                (1.0 - q)*tv1m
                // q true
                + q*tv2m
                )
            // p true
            - p*(v/tvp);

        float dw =
            // p false
            (1.0-p)*(1.0-w) / twm
            // p true
            - p*(w/twp);

        v += dt*dv;
        w += dt*dw;

        float jfi =
            // p true
            p * -v * (u - uc)*(1.0-u)/td;

        float jso =
            // p false
            (1.0 - p)*u / to
            // p true
            + p/tr;

        float jsi = -w * (1.0 + tanh(xk * (u-ucsi))) / (2.0 * tsi);

        // Apply stimulus
        float stim = 0.0;
        float stim_step = mod(float(step_count), period/dt);
        if (stim_step < stim_dur/dt) {
            stim = stim_f(stim_step * dt);
        }

        u -= (jfi+jso+jsi-stim)*dt;

        if (step_count > pre_pace_steps && !first_upstroke && u > align_thresh) {
            first_upstroke = true;
            start_comp = step_count;
            error = 0.0;
        }



        // Measure error
        if (first_upstroke && mod(float(step_count - start_comp), compare_stride) == 0.0) {
            float actual = texelFetch(data_texture, ivec2(data_index++, 0), 0).r;
            error += (u - actual)*(u - actual);
        }

        if (float((step_count - pre_pace_steps) - 1) / float((num_steps - pre_pace_steps) - 1) <= cc.x) {
            saved_value = u;
        }
    }

    error_texture = vec4(error, saved_value, 0, 1.0/period);
}
