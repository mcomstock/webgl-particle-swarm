/* global require */
require([
  'libs/Abubu.js',
  'text!data/zebrafish_onecl.txt',
  'text!shaders/run_simulation.frag',
  'text!shaders/reduce_error_s1.frag',
  'text!shaders/reduce_error_s2.frag',
  'text!shaders/update_velocities.frag',
  'text!shaders/update_particles.frag',
  'text!shaders/update_local_bests.frag',
  'text!shaders/copy_uint_texture.frag'
], function(
  Abubu,
  ActualData,
  RunSimulationShader,
  ReduceErrorS1Shader,
  ReduceErrorS2Shader,
  UpdateVelocitiesShader,
  UpdateParticlesShader,
  UpdateLocalBestsShader,
  CopyUIntTextureShader,
) {
  'use strict';

  //
  // Get canvas information
  //

  var canvas_1 = document.getElementById('canvas_1');
  var particles_width = parseInt(canvas_1.getAttribute('width'));
  var particles_height = parseInt(canvas_1.getAttribute('height'));

  var num_particles = particles_width * particles_height;

  //
  // Environment for solver with default values
  //

  var env = {
    simulation: {
      dt: 0.02,
      period: 400.0,
      stim_start: 2.0,
      stim_end: 7.0,
      stim_mag: 0.1,
      num_beats: 1,
      v_init: 1.0,
      w_init: 1.0,
    },
    particles: {
      phi_local: 2.05,
      phi_global: 2.05,
      global_bests: [
        [0.0, 0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0, 0.0],
      ],
      best_error_value: 1e10,
      chi: [],
      lower_bounds:[],
      upper_bounds:[],
      learning_rate: 0.1,
      // omega: 0.05,
      omega: 1.0,
    },
    bounds: [
      [[25, 200], [10, 300], [50, 500], [0.15, 0.4]],
      [[1, 20], [10, 50], [500, 1500], [5, 100]],
      [[5, 50], [1, 15], [0.2, 0.9], [0.1, 0.25]],
      [[0.005, 0.05], [0, 0], [0, 0], [0, 0]],
    ],
    velocity_update: {},
  };

  var phi = env.particles.phi_global + env.particles.phi_local;
  var chi = 0.05 * 2 / (phi - 2 + Math.sqrt(phi * (phi - 4)));
  // console.log("Chi is:\t"+chi);
  for (var i = 0; i < 4; ++i) {
    env.particles.chi.push(env.bounds[i].map(([min, max]) => chi * (max + min)/2));
  }
  // console.log(env.particles.chi);

  for(var i  =0; i < 4; i++)
  {
    var addArr_low = [];
    var addArr_up = [];
    for(var j = 0; j < 4; j++)
    {
      addArr_low.push(env.bounds[i][j][0]);
      addArr_up.push(env.bounds[i][j][1]);
    }
    env.particles.lower_bounds.push(addArr_low);
    env.particles.upper_bounds.push(addArr_up);
  }

  //
  // Set up data used to evaluate the each simulation run
  //

  var actual_data = ActualData.split('\n');
  var data_length = actual_data.length;
  var data_array = new Float32Array(data_length*4);

  // Pad out the extra pixel values. The data could be stored more densely by using the full pixel
  // value and by using a two-dimensional texture, but for now there is not enough to require that.
  var p = 0;
  for (var i = 0; i < data_length; ++i) {
    data_array[p++] = parseFloat(actual_data[i].trim());
    data_array[p++] = 0.0;
    data_array[p++] = 0.0;
    data_array[p++] = 0.0;
  }

  var init_array_1 = new Float32Array(num_particles * 4);
  var init_array_2 = new Float32Array(num_particles * 4);
  var init_array_3 = new Float32Array(num_particles * 4);
  var init_array_4 = new Float32Array(num_particles * 4);

  var v_init_arr_1 = new Float32Array(num_particles * 4);
  var v_init_arr_2 = new Float32Array(num_particles * 4);
  var v_init_arr_3 = new Float32Array(num_particles * 4);
  var v_init_arr_4 = new Float32Array(num_particles * 4);

  function random_init(i, j) {
    var [min, max] = env.bounds[i][j];
    return Math.random() * (max - min) + min;
  }

  var v_init_scale = 0.00005;
  // var v_init_scale = 0.5;
  function random_v_init(i,j)
  {
    // var [min, max] = env.bounds[i][j];
    // var realMin = -1 * (Math.abs(max - min));
    // var realMax = Math.abs(max-min);
    // return (Math.random() * (realMax - realMin) + realMin) * v_init_scale;
    // return (Math.random() * (realMax - realMin) + realMin);
    return Math.random() * v_init_scale;
    // return 0;
  }

  var p = 0;
  var p_v = 0;
  for (var i = 0; i < num_particles; ++i) {
    init_array_1[p] = random_init(0, 0);
    init_array_2[p] = random_init(1, 0);
    init_array_3[p] = random_init(2, 0);
    init_array_4[p++] = random_init(3, 0);

    init_array_1[p] = random_init(0, 1);
    init_array_2[p] = random_init(1, 1);
    init_array_3[p] = random_init(2, 1);
    init_array_4[p++] = 0;

    init_array_1[p] = random_init(0, 2);
    init_array_2[p] = random_init(1, 2);
    init_array_3[p] = random_init(2, 2);
    init_array_4[p++] = 0;

    init_array_1[p] = random_init(0, 3);
    init_array_2[p] = random_init(1, 3);
    init_array_3[p] = random_init(2, 3);
    init_array_4[p++] = 0;


    v_init_arr_1[p_v] = random_v_init(0,0);
    v_init_arr_2[p_v] = random_v_init(1,0);
    v_init_arr_3[p_v] = random_v_init(2,0);
    v_init_arr_4[p_v++] = random_v_init(3,0);
  

    v_init_arr_1[p_v] = random_v_init(0,1);
    v_init_arr_2[p_v] = random_v_init(1,1);
    v_init_arr_3[p_v] = random_v_init(2,1);
    v_init_arr_4[p_v++] = 0;

    v_init_arr_1[p_v] = random_v_init(0,2);
    v_init_arr_2[p_v] = random_v_init(1,2);
    v_init_arr_3[p_v] = random_v_init(2,2);
    v_init_arr_4[p_v++] = 0;

    v_init_arr_1[p_v] = random_v_init(0,3);
    v_init_arr_2[p_v] = random_v_init(1,3);
    v_init_arr_3[p_v] = random_v_init(2,3);
    v_init_arr_4[p_v++] = 0;


  }

  //
  // Textures
  //

  // The recorded data used to evaluate the correctness of simulation runs

  var data_texture = new Abubu.Float32Texture(data_length, 1, {
    pairable: true,
    data: data_array,
  });

  // These textures record the position (or particle values), velocity, and global best of each
  // particle



  var particles_texture_1 = new Abubu.Float32Texture(particles_width, particles_height, {
    pariable: true,
    data: init_array_1,
  });
  var particles_texture_2 = new Abubu.Float32Texture(particles_width, particles_height, {
    pariable: true,
    data: init_array_2,
  });
  var particles_texture_3 = new Abubu.Float32Texture(particles_width, particles_height, {
    pariable: true,
    data: init_array_3,
  });
  var particles_texture_4 = new Abubu.Float32Texture(particles_width, particles_height, {
    pariable: true,
    data: init_array_4,
  });

  var velocities_texture_1 = new Abubu.Float32Texture(particles_width, particles_height, {
    pariable: true,
    data: v_init_arr_1,
  });
  var velocities_texture_2 = new Abubu.Float32Texture(particles_width, particles_height, {
    pariable: true,
    data: v_init_arr_2,
  });
  var velocities_texture_3 = new Abubu.Float32Texture(particles_width, particles_height, {
    pariable: true,
    data: v_init_arr_3,
  });
  var velocities_texture_4 = new Abubu.Float32Texture(particles_width, particles_height, {
    pariable: true,
    data: v_init_arr_4,
  });

  var bests_texture_1 = new Abubu.Float32Texture(particles_width, particles_height, {
    pariable: true,
  });
  var bests_texture_2 = new Abubu.Float32Texture(particles_width, particles_height, {
    pariable: true,
  });
  var bests_texture_3 = new Abubu.Float32Texture(particles_width, particles_height, {
    pariable: true,
  });
  var bests_texture_4 = new Abubu.Float32Texture(particles_width, particles_height, {
    pariable: true,
  });

  // The out textures are used to get the values from updates

  var particles_out_texture_1 = new Abubu.Float32Texture(particles_width, particles_height, {
    pariable: true,
  });
  var particles_out_texture_2 = new Abubu.Float32Texture(particles_width, particles_height, {
    pariable: true,
  });
  var particles_out_texture_3 = new Abubu.Float32Texture(particles_width, particles_height, {
    pariable: true,
  });
  var particles_out_texture_4 = new Abubu.Float32Texture(particles_width, particles_height, {
    pariable: true,
  });

  var velocities_out_texture_1 = new Abubu.Float32Texture(particles_width, particles_height, {
    pariable: true,
  });
  var velocities_out_texture_2 = new Abubu.Float32Texture(particles_width, particles_height, {
    pariable: true,
  });
  var velocities_out_texture_3 = new Abubu.Float32Texture(particles_width, particles_height, {
    pariable: true,
  });
  var velocities_out_texture_4 = new Abubu.Float32Texture(particles_width, particles_height, {
    pariable: true,
  });

  var bests_out_texture_1 = new Abubu.Float32Texture(particles_width, particles_height, {
    pariable: true,
  });
  var bests_out_texture_2 = new Abubu.Float32Texture(particles_width, particles_height, {
    pariable: true,
  });
  var bests_out_texture_3 = new Abubu.Float32Texture(particles_width, particles_height, {
    pariable: true,
  });
  var bests_out_texture_4 = new Abubu.Float32Texture(particles_width, particles_height, {
    pariable: true,
  });


  // The error textures are used to reduce the error quantities of each particles from each
  // simulation run down to a global best.

  var local_best_error_init = new Float32Array(particles_width * particles_height * 4);
  // local_best_error_init.fill(100000, 0, particles_width * particles_height * 4);
  for(var i = 0; i < particles_width * particles_height * 4; i+=4)
  {
    local_best_error_init[i] = 100000;
  }

  var local_bests_error_texture_in = new Abubu.Float32Texture(particles_width, particles_height, {
    pariable: true,
    data: local_best_error_init,
  });
  // console.log(local_bests_error_texture_in.value);

  var local_bests_error_texture_out = new Abubu.Float32Texture(particles_width, particles_height, {
    pariable: true,
  });


  var error_texture = new Abubu.Float32Texture(particles_width, particles_height, {
    pariable: true,
  });


  var reduced_error_1_texture = new Abubu.Float32Texture(particles_width,1,{
    pairable: true,
  });
  var reduced_error_2_texture = new Abubu.Float32Texture(1,1,{
    pairable: true,
  });



  // var reduced_error_1_texture = new Abubu.Float32Texture(particles_width, particles_height, {
  //   pairable: true,
  // });

  // var reduced_error_2_texture = new Abubu.Float32Texture(particles_width, particles_height, {
  //   pairable: true,
  // });

  env.velocity_update.istate  = new Uint32Array(particles_width*particles_height*4);
  env.velocity_update.imat    = new Uint32Array(particles_width*particles_height*4);


  // Initialize the local bests with the inital positions.
  var init_bests_1 = new Abubu.Copy(particles_texture_1, bests_texture_1);
  var init_bests_2 = new Abubu.Copy(particles_texture_2, bests_texture_2);
  var init_bests_3 = new Abubu.Copy(particles_texture_3, bests_texture_3);
  var init_bests_4 = new Abubu.Copy(particles_texture_4, bests_texture_4);

  init_bests_1.render();
  init_bests_2.render();
  init_bests_3.render();
  init_bests_4.render();
  // console.log(bests_texture_1.value);


  var p=0;
  var seed = Date.now();
  // var seed = 0;
  var tm = new Abubu.TinyMT({vmat:0});

  for(var j=0 ; j<particles_height ; j++){
    for(var i=0 ; i<particles_width ; i++){
      //  mat1            mat2            seed
      tm.mat[0] = i ;     tm.mat[1] = j ; tm.mat[3] = seed ;
      tm.init() ;

      for(var k=0 ; k<4 ; k++){
        env.velocity_update.istate[p] = tm.state[k] ;
        env.velocity_update.imat[p] = tm.mat[k] ;
        p++ ;
      }
    }
  }

  env.velocity_update.ftinymtState = new Abubu.Uint32Texture(particles_width, particles_height, {
    data: env.velocity_update.istate,
    pairable : true,
  });

  env.velocity_update.stinymtState = new Abubu.Uint32Texture(particles_width, particles_height, {
    data: env.velocity_update.istate,
    pairable: true,
  });

  // mat state for each point of the generator .............................
  env.velocity_update.tinymtMat = new Abubu.Uint32Texture(particles_width, particles_height, {
    data: env.velocity_update.imat,
    pairable: true,
  });


  //
  // Simulation running solver
  //

  // Solver for running a simulation once for each particle
  var run_simulations_solver = new Abubu.Solver({
    fragmentShader: RunSimulationShader,
    uniforms: {
      in_particles_1: {
        type: 't',
        value: particles_texture_1,
      },
      in_particles_2: {
        type: 't',
        value: particles_texture_2,
      },
      in_particles_3: {
        type: 't',
        value: particles_texture_3,
      },
      in_particles_4: {
        type: 't',
        value: particles_texture_4,
      },
      data_texture: {
        type: 't',
        value: data_texture,
      },
      dt: {
        type: 'f',
        value: env.simulation.dt,
      },
      period: {
        type: 'f',
        value: env.simulation.period,
      },
      stim_start: {
        type: 'f',
        value: env.simulation.stim_start,
      },
      stim_end: {
        type: 'f',
        value: env.simulation.stim_end,
      },
      stim_mag: {
        type: 'f',
        value: env.simulation.stim_mag,
      },
      num_beats: {
        type: 'i',
        value: env.simulation.num_beats,
      },
      v_init: {
        type: 'f',
        value: env.simulation.v_init,
      },
      w_init: {
        type: 'f',
        value: env.simulation.w_init,
      },
    },
    targets: {
      error_texture: {
        location: 0,
        target: error_texture,
      },
    },
  });


  // var copy_tinymt_solver = new Abubu.Solver({
  //   fragmentShader: CopyUIntTextureShader,
  //   uniforms: {
  //     src_texture: {
  //       type: 't',
  //       value: env.velocity_update.stinymtState,
  //     },
  //   },
  //   targets: {
  //     dest_texture: {
  //       location: 0,
  //       target: env.velocity_update.ftinymtState,
  //     },
  //   },
  // });



  //
  // Solvers to find the global best simulation results
  //

  var reduce_error_1_solver = new Abubu.Solver({
    fragmentShader: ReduceErrorS1Shader,
    uniforms: {
      error_texture: {
        type: 't',
        value: error_texture,
      },
    },
    targets: {
      reduced_error_1: {
        location: 0,
        target: reduced_error_1_texture,
      },
    },
  });

  var reduce_error_2_solver = new Abubu.Solver({
    fragmentShader: ReduceErrorS2Shader,
    uniforms: {
      reduced_error_1: {
        type: 't',
        value: reduced_error_1_texture,
      },
    },
    targets: {
      reduced_error_2: {
        location: 0,
        target: reduced_error_2_texture,
      },
    },
  });

//
//  Solvers to update the local best positions
//

  /*

  So, this is for updating the local best errors. As it stands, this will copy
the new error value multiple times in the event that the new error is lower.
However, I tend to think that this redundancy will still perform better than
making a separate solver just to update the error?

*/
  function makeUpdateLocalBestsSolver(lbt, cvt, new_lbt)
  {
    return new Abubu.Solver({
      fragmentShader: UpdateLocalBestsShader,
      uniforms:{
        local_bests_texture: {
          type: 't',
          value: lbt,
        },
        local_bests_error_texture: {
          type: 't',
          value: local_bests_error_texture_in,
        },
        cur_vals_texture: {
          type: 't',
          value: cvt
        },
        cur_error_texture: {
          type: 't',
          value: error_texture,
        },
      },
      targets: {
        new_local_best: {
          location: 0,
          target: new_lbt,
        },
        new_local_best_error: {
          location: 1,
          target: local_bests_error_texture_out,
        },
      },
    });
  }

  var local_best_update_1_solver = makeUpdateLocalBestsSolver(bests_texture_1, particles_texture_1, bests_out_texture_1);
  var local_best_update_2_solver = makeUpdateLocalBestsSolver(bests_texture_2, particles_texture_2, bests_out_texture_2);
  var local_best_update_3_solver = makeUpdateLocalBestsSolver(bests_texture_3, particles_texture_3, bests_out_texture_3);
  var local_best_update_4_solver = makeUpdateLocalBestsSolver(bests_texture_4, particles_texture_4, bests_out_texture_4);



  //
  // Solvers to update the velocities
  //
  function makeVelocityUpdateSolver(pt, vt, bt, vto, num) {
    return new Abubu.Solver({
      fragmentShader: UpdateVelocitiesShader,
      uniforms: {
        positions_texture: {
          type: 't',
          value: pt,
        },
        velocities_texture: {
          type: 't',
          value: vt,
        },
        bests_texture: {
          type: 't',
          value: bt,
        },
        itinymtState: {
          type: 't',
          value: env.velocity_update.ftinymtState,
        },
        itinymtMat: {
          type: 't',
          value: env.velocity_update.tinymtMat,
        },
        chi: {
          type: 'v4',
          value: env.particles.chi[num],
        },
        phi_local: {
          type: 'f',
          value: env.particles.phi_local,
        },
        phi_global: {
          type: 'f',
          value: env.particles.phi_global,
        },
        global_best: {
          type: 'v4',
          value: env.particles.global_bests[num],
        },
        omega: {
          type: 'f',
          value: env.particles.omega,
        },
      },
      targets: {
        new_velocity: {
          location: 0,
          target: vto,
        },
        otinymtState: {
          location: 1,
          target: env.velocity_update.stinymtState,
        },
      },
    });
  }

  var velocity_1_solver = makeVelocityUpdateSolver(particles_texture_1, velocities_texture_1, bests_texture_1, velocities_out_texture_1, 0);
  var velocity_2_solver = makeVelocityUpdateSolver(particles_texture_2, velocities_texture_2, bests_texture_2, velocities_out_texture_2, 1);
  var velocity_3_solver = makeVelocityUpdateSolver(particles_texture_3, velocities_texture_3, bests_texture_3, velocities_out_texture_3, 2);
  var velocity_4_solver = makeVelocityUpdateSolver(particles_texture_4, velocities_texture_4, bests_texture_4, velocities_out_texture_4, 3);

  //
  // Solvers to update the positions
  //

  function makeParticleUpdateSolver(pt, vt, pto, num) {
    return new Abubu.Solver({
      fragmentShader: UpdateParticlesShader,
      uniforms: {
        positions_texture: {
          type: 't',
          value: pt,
        },
        velocities_texture: {
          type: 't',
          value: vt,
        },
        itinymtState: {
          type: 't',
          value: env.velocity_update.ftinymtState,
        },
        itinymtMat: {
          type: 't',
          value: env.velocity_update.tinymtMat,
        },
        lower_bounds: {
          type: 'v4',
          value: env.particles.lower_bounds[num],
        },
        upper_bounds: {
          type: 'v4',
          value: env.particles.upper_bounds[num],
        },
        learning_rate: {
          type: 'f',
          value: env.particles.learning_rate,
        },
      },
      targets: {
        new_position: {
          location: 0,
          target: pto,
        },
        otinymtState: {
          location: 1,
          target: env.velocity_update.stinymtState,
        },
      },
    });
  }

  var position_1_solver = makeParticleUpdateSolver(particles_texture_1, velocities_texture_1, particles_out_texture_1, 0);
  var position_2_solver = makeParticleUpdateSolver(particles_texture_2, velocities_texture_2, particles_out_texture_2, 1);
  var position_3_solver = makeParticleUpdateSolver(particles_texture_3, velocities_texture_3, particles_out_texture_3, 2);
  var position_4_solver = makeParticleUpdateSolver(particles_texture_4, velocities_texture_4, particles_out_texture_4, 3);

  var tinymt_copy = new Abubu.Copy(env.velocity_update.stinymtState, env.velocity_update.ftinymtState);

  var local_bests_1_copy = new Abubu.Copy(bests_out_texture_1, bests_texture_1);
  var local_bests_2_copy = new Abubu.Copy(bests_out_texture_2, bests_texture_2);
  var local_bests_3_copy = new Abubu.Copy(bests_out_texture_3, bests_texture_3);
  var local_bests_4_copy = new Abubu.Copy(bests_out_texture_4, bests_texture_4);

  var local_bests_error_copy = new Abubu.Copy(local_bests_error_texture_out, local_bests_error_texture_in);

  var positions_1_copy = new Abubu.Copy(particles_out_texture_1, particles_texture_1);
  var positions_2_copy = new Abubu.Copy(particles_out_texture_2, particles_texture_2);
  var positions_3_copy = new Abubu.Copy(particles_out_texture_3, particles_texture_3);
  var positions_4_copy = new Abubu.Copy(particles_out_texture_4, particles_texture_4);

  var velocities_1_copy = new Abubu.Copy(velocities_out_texture_1, velocities_texture_1);
  var velocities_2_copy = new Abubu.Copy(velocities_out_texture_2, velocities_texture_2);
  var velocities_3_copy = new Abubu.Copy(velocities_out_texture_3, velocities_texture_3);
  var velocities_4_copy = new Abubu.Copy(velocities_out_texture_4, velocities_texture_4);

  function update_global_best() {
    var [best_error, best_x_index, best_y_index] = reduced_error_2_texture.value.slice(-4, -1);
    // console.log(reduced_error_2_texture.value);
    // console.log(best_error);
    // console.log(best_x_index);
    // console.log(best_y_index);

    if (best_error < env.particles.best_error_value) {
      var best_particle_index = 4 * (best_y_index * particles_width + best_x_index);

      env.particles.global_bests[0] = particles_texture_1.value.slice(best_particle_index, best_particle_index + 4);
      env.particles.global_bests[1] = particles_texture_2.value.slice(best_particle_index, best_particle_index + 4);
      env.particles.global_bests[2] = particles_texture_3.value.slice(best_particle_index, best_particle_index + 4);
      env.particles.global_bests[3] = particles_texture_4.value.slice(best_particle_index, best_particle_index + 4);

      env.particles.best_error_value = best_error;
    }
  }

  function run() {
    run_simulations_solver.render();

    reduce_error_1_solver.render();
    reduce_error_2_solver.render();

    update_global_best();

    local_best_update_1_solver.render();
    local_best_update_2_solver.render();
    local_best_update_3_solver.render();
    local_best_update_4_solver.render();

    local_bests_1_copy.render();
    local_bests_2_copy.render();
    local_bests_3_copy.render();
    local_bests_4_copy.render();

    local_bests_error_copy.render();

    velocity_1_solver.render();
    env.velocity_update.ftinymtState.data = env.velocity_update.stinymtState.value;
    velocity_2_solver.render();
    env.velocity_update.ftinymtState.data = env.velocity_update.stinymtState.value;
    velocity_3_solver.render();
    env.velocity_update.ftinymtState.data = env.velocity_update.stinymtState.value;
    velocity_4_solver.render();
    env.velocity_update.ftinymtState.data = env.velocity_update.stinymtState.value;

    velocities_1_copy.render();
    velocities_2_copy.render();
    velocities_3_copy.render();
    velocities_4_copy.render();

    position_1_solver.render();
    env.velocity_update.ftinymtState.data = env.velocity_update.stinymtState.value;
    position_2_solver.render();
    env.velocity_update.ftinymtState.data = env.velocity_update.stinymtState.value;
    position_3_solver.render();
    env.velocity_update.ftinymtState.data = env.velocity_update.stinymtState.value;
    position_4_solver.render();
    env.velocity_update.ftinymtState.data = env.velocity_update.stinymtState.value;

    positions_1_copy.render();
    positions_2_copy.render();
    positions_3_copy.render();
    positions_4_copy.render();
  }

  
  for (var i = 0; i < 8; ++i) {
    console.log(env.particles.best_error_value);
    run();
  }

  var avg_v  = 0.0;
  for(var i = 0; i < particles_width*particles_height*4; i+=4)
  {
    avg_v += velocities_texture_1.value[i];
  }
  console.log("Avg vel:\t"+avg_v / (particles_width * particles_height * 4));


  var bestArr = [];
  for(var i = 0; i < 4; i++)
  {
    for(var j= 0; j < 4; j++)
    {
      bestArr.push(env.particles.global_bests[i][j]);
    }
  }
  console.log(bestArr);
  console.log(env.particles.best_error_value);

});
