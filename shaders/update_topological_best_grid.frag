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
    int tex_size = particle_dims.x * particle_dims.y;

    ivec2 locations[4];
    locations[0] = my_particle_idx.x == 0 ? ivec2(tex_size-1, my_particle_idx.y) : ivec2(my_particle_idx.x-1, my_particle_idx.y);
    locations[1] = my_particle_idx.x == tex_size-1 ? ivec2(0, my_particle_idx.y) : ivec2(my_particle_idx.x+1, my_particle_idx.y);
    locations[2] = my_particle_idx.y == 0 ? ivec2(my_particle_idx.x, tex_size-1) : ivec2(my_particle_idx.x, my_particle_idx.y-1);
    locations[3] = my_particle_idx.y == tex_size-1 ? ivec2(my_particle_idx.x, 0) : ivec2(my_particle_idx.x, my_particle_idx.y+1);

    float cur_error = texelFetch(local_bests_error_texture, my_particle_idx, 0).x;
    ivec2 best_idx = my_particle_idx;

    for (int i = 0; i < 4; ++i) {
        float new_error = texelFetch(local_bests_error_texture, locations[i], 0).x;
        if (new_error < cur_error) {
            cur_error = new_error;
            best_idx = locations[i];
        }
    }

    topological_best_idx = uvec4(best_idx, 0u, 0u);
}
