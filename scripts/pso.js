define('scripts/pso', [
  'scripts/gl_helper',
  'text!shaders/copy.frag',
  'text!shaders/copy_uint_texture.frag',
  'text!shaders/default.vert',
  'text!shaders/mitchell-schaeffer.frag',
  'text!shaders/expand_error.frag',
  'text!shaders/run_simulation_0d.frag',
  'text!shaders/reduce_error_s1.frag',
  'text!shaders/reduce_error_s2.frag',
  'text!shaders/update_velocities.frag',
  'text!shaders/update_particles.frag',
  'text!shaders/update_local_bests.frag',
  'text!shaders/round_float.frag',
  'text!shaders/hector_fhn.frag',
], function(
  GlHelper,
  CopyShader,
  CopyUintShader,
  DefaultVertexShader,
  MitchellSchaefferShader,
  ExpandErrorShader,
  RunSimulationShader,
  ReduceErrorS1Shader,
  ReduceErrorS2Shader,
  UpdateVelocitiesShader,
  UpdateParticlesShader,
  UpdateLocalBestsShader,
  RoundFloatShader,
  HectorFHNShader,
) {
  'use strict';

  return class Pso {
    constructor(particles_width, particles_height) {
      this.particles_width = particles_width;
      this.particles_height = particles_height;

      this.tex_width = 2 * particles_width;
      this.tex_height = 2 * particles_height;

      const canvas = document.createElement('canvas');
      canvas.width = this.particles_width;
      canvas.height = this.particles_height;

      this.gl_helper = new GlHelper(canvas);
    }

    static getEnv() {
      const env = {
        simulation: {
          model: 'fk',
          dt: 0.02,
          period: [],
          stim_start: 2.0,
          stim_end: 7.0,
          stim_mag: 0.1,
          num_beats: 1,
          pre_beats: 4,
          v_init: 1.0,
          w_init: 1.0,
          h_init: 1.0,
          align_thresh: [],
          trimmed_data: [],
          data_arrays: [],
          sample_interval: 1.0,
        },
        particles: {
          phi_local: 2.05,
          phi_global: 2.05,
          global_bests: [],
          best_error_value: 1e10,
          chi: [],
          lower_bounds:[],
          upper_bounds:[],
          learning_rate: 0.1,
          omega: 0.05,
          parameter_textures: 1,
        },
        fk_bounds: [
          [25, 10, 50, 0.15, 1, 10, 500, 5, 5, 1, 0.2, 0.1, 0.005],
          [200, 300, 500, 0.4, 20, 50, 1500, 100, 50, 15, 0.9, 0.1, 0.05],
        ],
        ms_bounds: [
          [5.0/3.0, 1.0/12.0, 75, 60, 0.065],
          [20.0/3.0, 1.0/3.0, 300, 240, 0.26],
        ],
        fhn_bounds: [
          [0.0142, 0.0759, 0.0072, 0.5176, 0.1906, -0.06, 1.0],
          [0.0142, 0.0759, 0.0072, 0.5176, 0.1906, -0.06, 1.0],
        ],
        velocity_update: {},
      };

      return env;
    }

    setupEnv(model, bounds, input_cls, pre_beats, num_beats, sample_interval) {
      this.env = Pso.getEnv();
      const env = this.env;

      env.simulation.model = model;
      env.simulation.period = input_cls;

      if (Number(num_beats)) {
        env.simulation.num_beats = Number(num_beats);
      }

      if (Number(pre_beats)) {
        env.simulation.pre_beats = Number(pre_beats);
      }

      if (Number(sample_interval)) {
        env.simulation.sample_interval = Number(sample_interval);
      }

      // 16 parameters per texture for now
      env.particles.parameter_textures = Math.ceil(bounds[0].length/16);

      const phi = env.particles.phi_global + env.particles.phi_local;
      const chi = 0.05 * 2 / (phi - 2 + Math.sqrt(phi * (phi - 4)));

      env.particles.chi = [];
      // The bounds arrays had better be the same length
      for (let i = 0; i < bounds[0].length; ++i) {
        env.particles.chi.push(chi * (bounds[1][i] - bounds[0][i])/2);
      }

      env.particles.lower_bounds = bounds[0];
      env.particles.upper_bounds = bounds[1];

      // Pad out these arrays so chunks of 16 can always be used as uniforms
      while (env.particles.chi.length % 16 !== 0) {
        env.particles.chi.push(0);
        env.particles.lower_bounds.push(0);
        env.particles.upper_bounds.push(0);
      }

      // Initialize all global best values to 0
      env.particles.global_bests = env.particles.lower_bounds.map(() => 0);
    }

    readData(raw_input_data, normalize) {
      const trimmed_data = [];
      const data_arrays = [];
      const align_thresh = [];

      const normalization = Number(normalize) || 1;
      const delta = 0.001;

      for (let i = 0; i < raw_input_data.length; ++i) {
        const raw_text = raw_input_data[i];

        const split_data = raw_text.split('\n');
        const actual_data = split_data.filter(x => !(x.trim() === ""));

        const full_parsed_data = actual_data.map(x => parseFloat(x.trim()));
        const maxVal = Math.max(...full_parsed_data);
        const full_normalized_data = full_parsed_data.map(x => (x * (normalization / maxVal)));

        const first_compare_index = full_normalized_data.findIndex(number => number > 0.15);

        const left_trimmed_data = full_normalized_data.slice(first_compare_index);

        // Pad out the extra pixel values. The data could be stored more densely by using the full pixel
        // value and by using a two-dimensional texture, but for now there is not enough to require that.
        const data_length = left_trimmed_data.length;
        const data_array = new Float32Array(4 * data_length);
        for (let j = 0; j < data_length; ++j) {
          data_array[4*j] = left_trimmed_data[j];
        }

        trimmed_data.push(left_trimmed_data);
        data_arrays.push(data_array);
        align_thresh.push(left_trimmed_data[0] - delta);
      }

      this.env.simulation.data_arrays = data_arrays;
      this.env.simulation.trimmed_data = trimmed_data;
      this.env.simulation.align_thresh = align_thresh;
    }

    initializeParticles() {
      const tex_width = this.tex_width;
      const tex_height = this.tex_height;
      const asize = 4 * tex_width * tex_height;
      const init_arrays = [];
      const { lower_bounds, upper_bounds } = this.env.particles;

      const get_random = (idx) => Math.random() * (upper_bounds[idx] - lower_bounds[idx]) + lower_bounds[idx];

      for (let tex = 0; tex < lower_bounds.length / 16; ++tex) {
        const init_array = new Float32Array(asize);

        for (let i = 0; i < tex_width; ++i) {
          for (let j = 0; j < tex_height; ++j) {
            const idx = tex*16 + Math.floor(j/this.particles_height)*8 + Math.floor(i/this.particles_width)*4;

            for (let p = 0; p < 4; ++p) {
              init_array[4*(tex_width*j+i)+p] = get_random(idx+p);
            }
          }
        }

        init_arrays.push(init_array);
      }

      return init_arrays;
    }

    initializeTextures() {
      const gl_helper = this.gl_helper;
      const particles_width = this.particles_width;
      const particles_height = this.particles_height;
      const tex_width = this.tex_width;
      const tex_height = this.tex_height;
      const { num_beats, period, sample_interval } = this.env.simulation;

      const data_arrays = this.env.simulation.data_arrays;
      const init_arrays = this.initializeParticles();
      const zero_array = new Float32Array(tex_width*tex_height*4);

      this.simulation_lengths = [];
      this.data_textures = [];

      console.log("Period length is "+period.length);
      for (let i = 0; i < period.length; i++) {
        this.simulation_lengths.push(Math.ceil(Math.ceil(num_beats * period[i]) / sample_interval));
        this.data_textures.push(gl_helper.loadFloatTexture(data_arrays[i].length/4, 1, data_arrays[i]));
      }

      this.particles_textures = [];
      this.velocities_textures = [];
      this.bests_textures = [];
      this.particles_out_textures = [];
      this.velocities_out_textures = [];
      this.bests_out_textures = [];

      for (const init_array of init_arrays) {
        this.particles_textures.push(gl_helper.loadFloatTexture(tex_width, tex_height, init_array));
        this.velocities_textures.push(gl_helper.loadFloatTexture(tex_width, tex_height, zero_array));
        this.bests_textures.push(gl_helper.loadFloatTexture(tex_width, tex_height, zero_array));
        this.particles_out_textures.push(gl_helper.loadFloatTexture(tex_width, tex_height, null));
        this.velocities_out_textures.push(gl_helper.loadFloatTexture(tex_width, tex_height, null));
        this.bests_out_textures.push(gl_helper.loadFloatTexture(tex_width, tex_height, null));
      }

      // The error textures are used to reduce the error quantities of each particles from each
      // simulation run down to a global best.
      const local_error_init = new Float32Array(tex_width * tex_height * 4);
      for (let i = 0; i < tex_width * tex_height * 4; i += 4) {
        local_error_init[i] = 100000.0;
      }

      this.local_bests_error_texture = gl_helper.loadFloatTexture(tex_width, tex_height, local_error_init);
      this.local_bests_error_texture_out = gl_helper.loadFloatTexture(tex_width, tex_height, null);

      this.error_texture = gl_helper.loadFloatTexture(particles_width, particles_height, null);
      this.expanded_error_texture = gl_helper.loadFloatTexture(tex_width, tex_height, null);
      this.simulation_texture = gl_helper.loadFloatTexture(Math.max(...this.simulation_lengths), 1, null);

      this.error_sum_texture_0 = gl_helper.loadFloatTexture(particles_width, particles_height, null);
      this.error_sum_texture_1 = gl_helper.loadFloatTexture(particles_width, particles_height, null);

      this.reduced_error_1_texture = gl_helper.loadFloatTexture(particles_width, particles_height, null);
      this.reduced_error_2_texture = gl_helper.loadFloatTexture(particles_width, particles_height, null);

      const env = this.env;

      env.velocity_update.istate  = new Uint32Array(tex_width*tex_height*4);
      env.velocity_update.imat    = new Uint32Array(tex_width*tex_height*4);

      // Modifies the entries of state
      const next_tmt_state = (mat, state) => {
        let x = (state[0] & 0x7fffffff) ^ state[1] ^ state[2];
        let y = state[3];

        x ^= (x << 1);
        y ^= (y >>> 1) ^ x;
        state[0] = state[1];
        state[1] = state[2];
        state[2] = x ^ (y << 10);
        state[3] = y;

        state[1] ^= (-(y & 1) >>> 0) & mat[0];
        state[2] ^= (-(y & 1) >>> 0) & mat[1];
      };

      let p = 0;
      const seed = Date.now();
      const mat = [0, 0, 0, seed];
      const state = [0, 0, 0, 0];
      for (let j = 0; j < tex_height; ++j) {
        for (let i = 0; i < tex_width; ++i) {
          mat[0] = i;
          mat[1] = j;

          state[0] = mat[3];
          state[1] = mat[0];
          state[2] = mat[1];
          state[3] = mat[2];

          for (let k = 1; k < 8; ++k) {
            const a = k & 3;
            const b = (k-1) & 3;
            state[a] ^= k + Math.imul(1812433253, (state[b] ^ (state[b] >>> 30)));
          }

          for (let k = 0; k < 8; ++k) {
            next_tmt_state(mat, state);
          }

          for (let k = 0; k < 4; ++k) {
            env.velocity_update.istate[p] = state[k];
            env.velocity_update.imat[p] = mat[k];
            p++;
          }
        }
      }

      env.velocity_update.ftinymtState = gl_helper.loadUintTexture(tex_width, tex_height, env.velocity_update.istate);
      env.velocity_update.stinymtState = gl_helper.loadUintTexture(tex_width, tex_height, env.velocity_update.istate);
      // mat state for each point of the generator .............................
      env.velocity_update.tinymtMat = gl_helper.loadUintTexture(tex_width, tex_height, env.velocity_update.imat);
    }

    getDefaultShaderMap() {
      const makeUpdateLocalBestsSolver = (num) => {
        return {
          vert: DefaultVertexShader,
          frag: UpdateLocalBestsShader,
          uniforms: [
            ['local_bests_texture', 'tex', () => this.bests_textures[num]],
            ['local_bests_error_texture', 'tex', () => this.local_bests_error_texture],
            ['cur_vals_texture', 'tex', () => this.particles_textures[num]],
            ['cur_error_texture', 'tex', () => this.expanded_error_texture],
          ],
          out: [this.bests_out_textures[num], this.local_bests_error_texture_out],
          run: this.gl_helper.runBigProgram,
        };
      };

      const makeVelocityUpdateSolver = (num) => {
        return {
          vert: DefaultVertexShader,
          frag: UpdateVelocitiesShader,
          uniforms: [
            ['positions_texture', 'tex', () => this.particles_textures[num]],
            ['velocities_texture', 'tex', () => this.velocities_textures[num]],
            ['bests_texture', 'tex', () => this.bests_textures[num]],
            ['itinymtState', 'tex', () => this.env.velocity_update.ftinymtState],
            ['itinymtMat', 'tex', () => this.env.velocity_update.itinymtMat],
            ['chi', '4fv_a', () => [this.env.particles.chi, num*16, 16]],
            ['phi_local', '1f', () => this.env.particles.phi_local],
            ['phi_global', '1f', () => this.env.particles.phi_global],
            ['global_best', '4fv_a', () => [this.env.particles.global_bests, num*16, 16]],
            ['omega', '1f', () => this.env.particles.omega],
          ],
          out: [this.velocities_out_textures[num], this.env.velocity_update.stinymtState],
          run: this.gl_helper.runBigProgram,
        };
      };

      const makeParticleUpdateSolver = (num) => {
        return {
          vert: DefaultVertexShader,
          frag: UpdateParticlesShader,
          uniforms: [
            ['positions_texture', 'tex', () => this.particles_textures[num]],
            ['velocities_texture', 'tex', () => this.velocities_textures[num]],
            ['itinymtState', 'tex', () => this.env.velocity_update.ftinymtState],
            ['itinymtMat', 'tex', () => this.env.velocity_update.itinymtMat],
            ['lower_bounds', '4fv_a', () => [this.env.particles.lower_bounds, num*16, 16]],
            ['upper_bounds', '4fv_a', () => [this.env.particles.upper_bounds, num*16, 16]],
            ['learning_rate', '1f', () => this.env.particles.learning_rate],
          ],
          out: [this.particles_out_textures[num], this.env.velocity_update.stinymtState],
          run: this.gl_helper.runBigProgram,
        };
      };

      const makeCopySolver = (original, copy, idx) => {
        return {
          vert: DefaultVertexShader,
          frag: CopyShader,
          uniforms: [
            ['original', 'tex', idx === undefined ? () => this[original] : () => this[original][idx]],
          ],
          out: [idx === undefined ? this[copy] : this[copy][idx]],
          run: this.gl_helper.runBigProgram,
        };
      };

      const makeRunSimulationSolver = (final) => {

        let model_frag;
        switch (String(this.env.simulation.model)) {
          case 'fk':
            model_frag = RunSimulationShader;
            break;
          case 'ms':
            model_frag = MitchellSchaefferShader;
            break;
          case 'fhn':
            model_frag = HectorFHNShader;
            break;
          default:
            console.log("How could no model be selected oh no!");
        }

        const solver = {
          vert: DefaultVertexShader,
          frag: model_frag,
          uniforms: [
            ['data_texture', 'tex', (cl_idx) => this.data_textures[cl_idx]],
            ['dt', '1f', () => this.env.simulation.dt],
            ['period', '1f', (cl_idx) => this.env.simulation.period[cl_idx]],
            ['stim_start', '1f', () => this.env.simulation.stim_start],
            ['stim_end', '1f', () => this.env.simulation.stim_end],
            ['stim_mag', '1f', () => this.env.simulation.stim_mag],
            ['num_beats', '1i', () => this.env.simulation.num_beats],
            ['pre_beats', '1i', () => this.env.simulation.pre_beats],
            ['align_thresh', '1f', (cl_idx) => this.env.simulation.align_thresh[cl_idx]],
            ['sample_interval', '1f', () => this.env.simulation.sample_interval],
          ],
          out: [final ? this.simulation_texture : this.error_texture],
          run: final ? this.gl_helper.runFinal : this.gl_helper.runSimulation,
        };

        for (let i = 0; i < this.particles_textures.length; ++i) {
          if (final) {
            solver.uniforms.push(['in_particles_' + (i+1), 'tex', () => this.final_particles_textures[i]]);
          } else {
            solver.uniforms.push(['in_particles_' + (i+1), 'tex', () => this.particles_textures[i]]);
          }
        }

        if (this.env.simulation.model === 'ms') {
          solver.uniforms.push(['h_init', '1f', () => this.env.simulation.h_init]);
        } else {
          solver.uniforms.push(['v_init', '1f', () => this.env.simulation.v_init]);
          solver.uniforms.push(['w_init', '1f', () => this.env.simulation.w_init]);
        }

        return solver;
      };

      const shader_map = {
        run_simulation: makeRunSimulationSolver(false),

        run_final_simulation: makeRunSimulationSolver(true),

        reduce_error_1: {
          vert: DefaultVertexShader,
          frag: ReduceErrorS1Shader,
          uniforms: [
            ['error_texture', 'tex', () => this.error_texture],
          ],
          out: [this.reduced_error_1_texture],
          run: this.gl_helper.runProgram,
        },

        reduce_error_2: {
          vert: DefaultVertexShader,
          frag: ReduceErrorS2Shader,
          uniforms: [
            ['reduced_error_1', 'tex', () => this.reduced_error_1_texture],
          ],
          out: [this.reduced_error_2_texture],
          run: this.gl_helper.runProgram,
        },

        expand_error: {
          vert: DefaultVertexShader,
          frag: ExpandErrorShader,
          uniforms: [
            ['error_texture', 'tex', () => this.error_texture],
          ],
          out: [this.expanded_error_texture],
          run: this.gl_helper.runBigProgram,
        },

        tinymt_copy: {
          vert: DefaultVertexShader,
          frag: CopyUintShader,
          uniforms: [
            ['original', 'tex', () => this.env.velocity_update.stinymtState],
          ],
          out: [this.env.velocity_update.ftinymtState],
          run: this.gl_helper.runBigProgram,
        },

        local_error_copy: makeCopySolver('local_bests_error_texture_out', 'local_bests_error_texture'),
      };

      for (let i = 0; i < this.env.particles.parameter_textures; ++i) {
        shader_map['local_best_update_' + i] = makeUpdateLocalBestsSolver(i);
        shader_map['velocity_' + i] = makeVelocityUpdateSolver(i);
        shader_map['position_' + i] = makeParticleUpdateSolver(i);
        shader_map['local_bests_copy_' + i] = makeCopySolver('bests_out_textures', 'bests_textures', i);
        shader_map['positions_copy_' + i] = makeCopySolver('particles_out_textures', 'particles_textures', i);
        shader_map['velocities_copy_' + i] = makeCopySolver('velocities_out_textures', 'particles_textures', i);
      }

      return shader_map;
    }

    setupAllSolvers() {
      const shader_map = this.getDefaultShaderMap();
      const program_map = {};

      this.gl_helper.initDefaultVertexBuffer();

      for (const key in shader_map) {
        program_map[key] = this.gl_helper.setupDefault(shader_map[key], this);
      }

      this.program_map = program_map;
    }

    updateGlobalBest() {
      const env = this.env;

      const error_array = new Float32Array(this.particles_width*this.particles_height*4);
      this.gl_helper.getFloatTextureArray(this.reduced_error_2_texture, this.particles_width, this.particles_height, error_array);

      const [best_error, best_x_index, best_y_index] = error_array.slice(-4, -1);

      const particles_array = new Float32Array(this.tex_width*this.tex_height*4);
      if (best_error < env.particles.best_error_value) {
        const best_particle_index = [
          4 * (best_y_index * this.tex_width + best_x_index),
          4 * (best_y_index * this.tex_width + best_x_index + this.particles_width),
          4 * ((best_y_index + this.particles_width) * this.tex_width + best_x_index),
          4 * ((best_y_index + this.particles_width) * this.tex_width + best_x_index + this.particles_width),
        ];

        for (let i = 0; i < this.particles_textures.length; ++i) {
          this.gl_helper.getFloatTextureArray(this.particles_textures[i], this.tex_width, this.tex_height, particles_array);
          for (let j = 0; j < 4; ++j) {
            for (let k = 0; k < 4; ++k) {
              env.particles.global_bests[16*i+4*j+k] = particles_array[best_particle_index[j]+k];
            }
          }
        }

        env.particles.best_error_value = best_error;
      }
    }

    runOneIteration() {
      const program_map = this.program_map;

      for (let i = 0; i < this.env.simulation.period.length; ++i) {
        program_map.run_simulation(i, i === 0);
      }

      program_map.reduce_error_1();
      program_map.reduce_error_2();

      program_map.expand_error();

      this.updateGlobalBest();

      for (let i = 0; i < this.env.particles.parameter_textures; ++i) {
        program_map['local_best_update_' + i]();
      }

      for (let i = 0; i < this.env.particles.parameter_textures; ++i) {
        program_map['local_bests_copy_' + i]();
      }

      program_map.local_error_copy();

      for (let i = 0; i < this.env.particles.parameter_textures; ++i) {
        program_map['velocity_' + i]();
        program_map.tinymt_copy();
      }

      for (let i = 0; i < this.env.particles.parameter_textures; ++i) {
        program_map['velocities_copy_' + i]();
      }

      for (let i = 0; i < this.env.particles.parameter_textures; ++i) {
        program_map['position_' + i]();
        program_map.tinymt_copy();
      }

      for (let i = 0; i < this.env.particles.parameter_textures; ++i) {
        program_map['positions_copy_' + i]();
      }
    }

    setFinalPosition(sim_length, values) {
      const final_width = 2*sim_length;
      const final_height = 2;
      const final_size = 4 * final_width * final_height;

      const num_textures = Math.ceil(values.length / 16);

      this.final_particles_textures = [];

      for (let t = 0; t < num_textures; ++t) {
        const particles_array = new Float32Array(final_size);

        for (let i = 0; i < sim_length; ++i) {
          for (let j = 0; j < 4; ++j) {
            for (let k = 0; k < 4; ++k) {
              particles_array[4*(i+j*sim_length)+k] = values[16*t + 4*j + k];
            }
          }
        }

        this.final_particles_textures.push(this.gl_helper.loadFloatTexture(final_width, final_height, particles_array));
      }
    }

    runFinalSimulationSolver(cl_idx, values) {
      const simsize = this.simulation_lengths[cl_idx];
      const texsize = simsize*4;

      const parameter_values = values || this.env.particles.global_bests;

      while (parameter_values.length % 16 !== 0) {
        parameter_values.push(0);
      }

      this.setFinalPosition(simsize, parameter_values);
      this.program_map.run_final_simulation(cl_idx, simsize);

      const texture_array = new Float32Array(texsize);
      this.gl_helper.getFloatTextureArray(this.simulation_texture, simsize, 1, texture_array);

      const simulation_data = new Float32Array(simsize);
      for (let i = 0; i < simsize; ++i) {
        simulation_data[i] = texture_array[4*i+1];
      }

      return simulation_data;
    }
  };
});
