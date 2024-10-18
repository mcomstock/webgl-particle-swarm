#version 300 es

precision highp float;
precision highp int;

uniform sampler2D local_bests_error_texture;

layout (location = 0) out uvec4 topological_best_idx;

in vec2 cc;

void main() {
    ivec2 texture_dims = textureSize(local_bests_error_texture, 0);
    ivec2 particle_dims = texture_dims / 2;
    ivec2 my_particle_idx = ivec2(floor(cc*vec2(particle_dims)));
    int my_particle_num = my_particle_idx.y * particle_dims.x + my_particle_idx.x;
    int tex_size = particle_dims.x * particle_dims.y;

    int pnums[2];
    pnums[0] = my_particle_num == 0 ? tex_size - 1 : my_particle_num - 1;
    pnums[1] = my_particle_num == tex_size-1 ? 0 : my_particle_num + 1;

    float cur_error = texelFetch(local_bests_error_texture, my_particle_idx, 0).x;
    ivec2 best_idx = my_particle_idx;

    for (int i = 0; i < 2; ++i) {
        ivec2 location = ivec2(pnums[i] % particle_dims.x, pnums[i] / particle_dims.x);
        float new_error = texelFetch(local_bests_error_texture, location, 0).x;
        if (new_error < cur_error) {
            cur_error = new_error;
            best_idx = location;
        }
    }

    topological_best_idx = uvec4(best_idx, 0u, 0u);
}
