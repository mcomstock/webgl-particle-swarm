#version 300 es

precision highp float;
precision highp int;

uniform sampler2D in_particles_1, in_particles_2;
uniform sampler2D data_texture;
uniform sampler2D var_init_texture;

layout (location = 0) out vec4 error_texture;

in vec2 cc;

// Simulation parameters
uniform float dt, period;
uniform int num_beats, pre_beats, data_type;
uniform float v_init, w_init;
uniform float align_thresh;
uniform float sample_interval, apd_thresh, weight;
uniform float stim_dur, stim_mag, stim_offset_1, stim_offset_2, stim_t_scale;
uniform bool stim_biphasic;

#define thwinf thw
#define thvinf thvm
#define thso thw
#define thsi thso

float biphasic_stim_f(const float t) {
    float a = (t/stim_t_scale - stim_offset_2);

    return -stim_mag * (t / stim_t_scale - stim_offset_1) / (1.0 + a*a*a*a);
}

float square_stim_f(const float t) {
    return stim_mag;
}

void main() {
    // PSO derived parameters
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

    // 4v params
    float thv   = particles_1.r;
    float tv1m  = particles_1.g;
    float tv2m  = particles_1.b;
    float tvp   = particles_1.a;

    float uwm   = particles_2.r;
    float tso1  = particles_2.g;
    float kso   = particles_2.b;
    float ts1   = particles_2.a;

    float ts2   = particles_3.r;
    float ks    = particles_3.g;
    float tw1m  = particles_3.b;
    float tw2m  = particles_3.a;

    float tw1p  = particles_4.r;
    float tw2p  = tw1p;
    float tfi   = particles_4.g;
    float to1   = particles_4.b;
    float to2   = particles_4.a;

    float tso2  = particles_5.r;
    float uso   = particles_5.g;
    float us   = particles_5.b;
    float tsi1  = particles_5.a;
    float tsi2  = tsi1;

    float thw   = particles_6.r;
    float thvm  = particles_6.g;
    float tho   = particles_6.b;
    float kwm   = particles_6.a;

    float twinf     = particles_7.r;
    float winfstar  = particles_7.g;
    float uu        = particles_7.b;

    float kwp   = 5.7;
    float wcp   = 0.15;
    float ksi   = 97.8;
    float sc    = 0.007;
    float uo = 0.0;

    vec4 var_init = texture(var_init_texture, cc);

    float u = var_init[0];
    float v = var_init[1];
    float w = var_init[2];
    float s = var_init[3];

    float uNew;

    float dv;
    float dw;
    float ds;

    float tvm;
    float ts;
    float to;
    float twp;
    float twm;
    float tso;
    float tsi;
    float vinf;
    float winf;

    float xfi;
    float xso;
    float xsi;

    float compare_stride = round(sample_interval / dt);

    // TODO: should really be doing error such that we can always start at
    // 0.0. I'll update the curve error part to allow this later.
    float error = (data_type == 1) ? 0.0 : 10000000000.0;

    int data_index = 0;

    int start_comp = 0;
    int compared_points = 0;
    bool first_align_upstroke = false;
    float saved_value = -1.0;

    float APD_start, APD_end;

    bool activated = false;

    // Run the simulation with the current swarm parameters
    for (int step_count = 1; step_count <= num_steps; ++step_count) {
        // Begin 4v update
        tvm = (u >= thvm) ? tv2m : tv1m;
        ts = (u >= thw) ? ts2 : ts1;
        to = (u >= tho) ? to2 : to1;

        twp = tw1p + (tw2p - tw1p) * (1.0+tanh((w-wcp)*kwp))/2.0;
        twm = tw1m + (tw2m - tw1m) * (1.0+tanh((u-uwm)*kwm))/2.0;
        tso = tso1 + (tso2 - tso1) * (1.0+tanh((u-uso)*kso))/2.0;
        tsi = tsi1 + (tsi2 - tsi1) * (1.0+tanh((s-sc)*ksi))/2.0;

        vinf = (u >= thvinf) ? 0.0 : 1.0;
        winf = (u >= thwinf) ? winfstar : (1.0-u/twinf);

        // Gate Evolution
        dv = (u >= thv) ? (-(v/tvp)) : ((vinf-v)/tvm);
        dw = (u >= thw) ? (-(w/twp)) : ((winf-w)/twm);
        ds = ((1.0+tanh((u-us)*ks))/2.0-s)/ts;

        v += dt * dv;
        w += dt * dw;
        s += dt * ds;

        // Currents
        xfi = (u >= thv) ? -v*(u-thv)*(uu-u)/tfi : 0.0;
        xso = (u >= thso) ? 1.0/tso : (u-uo)/to;
        xsi = (u >= thsi) ? -w*s/tsi : 0.0;

        float stim = 0.0;
        float stim_t = mod(float(step_count)*dt, period);
        if (stim_t < stim_dur) {
            if (stim_biphasic) {
                stim = biphasic_stim_f(stim_t);
            } else {
                stim = square_stim_f(stim_t);
            }
        }

        float prev_u = u;
        u = u - dt*(xfi + xso + xsi - stim);

        // End 4v update

        if (step_count > pre_pace_steps) {
            //APD only mode
            if (data_type == 1) {
                if (!activated && u > apd_thresh) {
                    activated = true;
                    float x0 = float(step_count-1)*dt;
                    float x1 = float(step_count)*dt;

                    float y0 = prev_u;
                    float y1 = u;

                    // Linear interpolation of actual crossing of threshold
                    APD_start = (x0*(y1 - apd_thresh) + x1*(apd_thresh - y0)) / (y1-y0);
                } else if (activated && u < apd_thresh) {
                    activated = false;

                    float x0 = float(step_count-1)*dt;
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

    error_texture = vec4(error, saved_value, 0, compared_points == 0 ? weight : weight / float(compared_points));
}
