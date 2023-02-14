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

// % this function should give stimulus value from time in ms
// %stimscale = 0.35*180;
// stimscale = 0.1;
// stimdur = 10;
// offset1 = 7;
// offset2 = offset1*.96;
// tscale = 0.725;
// f = @(t) -stimscale*(t/tscale-offset1)./(1+(t/tscale-offset2).^4);


float stim_f(const float t) {
    const float stim_scale = 0.1 * 200.0;
    const float stim_dur = 2.0;
    const float offset_1 = 7.0;
    const float offset_2 = offset_1 * 0.96;
    const float t_scale = 0.725;

    return 2.0 * (-stim_scale * (t/t_scale - offset_1) / (1.0 + pow((t/t_scale - offset_2), 4.0)));
}

void main() {
    // PSO derived parameters
    int num_period = int(ceil(period/dt));
    int total_beats = pre_beats + num_beats;
    float endtime = ceil(float(total_beats)*period);
    float pre_pace_endtime = ceil(float(pre_beats)*period);
    int pre_pace_steps = int(ceil(pre_pace_endtime/dt));
    int num_steps = int(ceil(endtime/dt));

    // const float stim_dur = 10.0;
    const float stim_dur = 2.0;

    ivec2 tex_size = textureSize(in_particles_1, 0);
    ivec2 idx = ivec2(floor(cc * 0.5 * vec2(tex_size)));

    vec4 particles_1 = texelFetch(in_particles_1, idx, 0);
    vec4 particles_2 = texelFetch(in_particles_1, idx + ivec2(tex_size.x/2, 0), 0);

    float alpha = particles_1.r;
    float beta = particles_1.g;
    float eps = particles_1.b;
    float mu = particles_1.a;

    float gamma = particles_2.r;
    float theta = particles_2.g;
    float delta = particles_2.b;

    // float v = 0.0;
    // float h = h_init;

    // float u = 0.0;
    // float u = 0.25;
    // float v = 0.05;
    float u = 0.0;
    float v = 0.0;

    // float f, m, thf, ah, bh, stim, stim_step, dv, dh;

    // float kappa = 100.0;

    float du, dv, stim, stim_step;

    float compare_stride = round(sample_interval / dt);

    float error = 0.0;
    error = 10000000000.0;
    // error = 100000.0;

    int data_index = 0;

    int start_comp = 0;
    bool first_upstroke = false;
    float saved_value = -1.0;

    // Run the simulation with the current swarm parameters
    for (int step_count = 1; step_count <= num_steps; ++step_count) {
        // m = max(0.0, v);
        // m = min(1.0, m);

        // f = 0.5 * (1.0 + tanh(kappa * (v - vgate)));
        // thf = topen + (tclose - topen) * f;
        // ah = (1.0 - f) / thf;
        // bh = f / thf;

        stim = 0.0;
        stim_step = mod(float(step_count), period/dt);
        if (stim_step < stim_dur/dt) {
            stim = stim_f(stim_step * dt);
        }


        // if(mod(float(step_count),float(num_period)) == 0.0)
        // {
            // u = 0.5;
        // }

        // eps      [.001, 1]
        // alpha    [0.05, 0.6]
        // beta     [0.2, 2]
        // gamma    [0.01, 1]
        // mu       [0.2, 2]
        // theta    [-0.1, 0.1]
        // delta    [0.5, 1.5]


        du = stim + mu*u*(1.0-u)*(u-alpha)-u*v;
        dv = eps * ( (beta - u) * (u - gamma) - delta*v - theta );


        // du = mu*u*(1.0-u)*(u-alpha)-u*v;
        // dv = eps * ( (beta - u) * (u - gamma) - delta*v - theta );


        // u += dt * du;
        // v += dt * dv;

        u = u + dt*du;
        v = v + dt*dv;

        // dv = stim + gna * m * m * h * (1.0 - v) - gk * v;
        // dh = ah * (1.0 - h) - bh * h;

        // v += dv * dt;
        // h += dh * dt;

        if (step_count > pre_pace_steps && !first_upstroke && u > align_thresh) {
            first_upstroke = true;
            start_comp = step_count;
            error = 0.0;
        }

        // Measure error
        if (first_upstroke && mod(float(step_count - start_comp), compare_stride) == 0.0) {
            float actual = texelFetch(data_texture, ivec2(data_index++, 0), 0).r;
            error += (u - actual) * (u - actual);
        }

        if (float((step_count - pre_pace_steps) - 1) / float((num_steps - pre_pace_steps) - 1) <= cc.x) {
            saved_value = u;
        }
    }

    error_texture = vec4(error, saved_value, 0, 1.0/period);
}
