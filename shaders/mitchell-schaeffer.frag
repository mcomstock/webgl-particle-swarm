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
uniform float h_init;
uniform float align_thresh;
uniform float sample_interval;

float stim_f(const float t) {
    const float stim_scale = 0.4;
    const float stim_dur = 10.0;
    const float offset_1 = 7.0;
    const float offset_2 = offset_1 * 0.96;
    const float t_scale = 0.725;

    return (-stim_scale * (t/t_scale - offset_1) / (1.0 + pow((t/t_scale - offset_2), 4.0)));
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

    float gna = particles_1.r;
    float gk = particles_1.g;
    float tclose = particles_1.b;
    float topen = particles_1.a;
    float vgate = particles_2.r;

    float v = 0.0;
    float h = h_init;
    float f, m, thf, ah, bh, stim, stim_step, dv, dh;
    float kappa = 100.0;

    float compare_stride = round(sample_interval / dt);

    float error = 10000000000.0;

    int data_index = 0;

    int start_comp = 0;
    bool first_upstroke = false;
    float saved_value = -1.0;

    // Run the simulation with the current swarm parameters
    for (int step_count = 1; step_count <= num_steps; ++step_count) {
        m = max(0.0, v);
        m = min(1.0, m);

        f = 0.5 * (1.0 + tanh(kappa * (v - vgate)));
        thf = topen + (tclose - topen) * f;
        ah = (1.0 - f) / thf;
        bh = f / thf;

        stim = 0.0;
        stim_step = mod(float(step_count), period/dt);
        if (stim_step < stim_dur/dt) {
            stim = stim_f(stim_step * dt);
        }

        dv = stim + gna * m * m * h * (1.0 - v) - gk * v;
        dh = ah * (1.0 - h) - bh * h;

        v += dv * dt;
        h += dh * dt;

        if (step_count > pre_pace_steps && !first_upstroke && v > align_thresh) {
            first_upstroke = true;
            start_comp = step_count;
            error = 0.0;
        }

        // Measure error
        if (first_upstroke && mod(float(step_count - start_comp), compare_stride) == 0.0) {
            float actual = texelFetch(data_texture, ivec2(data_index++, 0), 0).r;
            error += (v - actual) * (v - actual);
        }

        if (float((step_count - pre_pace_steps) - 1) / float((num_steps - pre_pace_steps) - 1) <= cc.x) {
            saved_value = v;
        }
    }

    error_texture = vec4(error, saved_value, 0, 1.0/period);
}
