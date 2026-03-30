define('scripts/pso', [
  'scripts/gl_helper',
  'text!shaders/copy.frag',
  'text!shaders/copy_uint_texture.frag',
  'text!shaders/default.vert',
  'text!shaders/mitchell-schaeffer.frag',
  'text!shaders/modified-mitchell-schaeffer.frag',
  'text!shaders/expand_error.frag',
  'text!shaders/fenton_karma.frag',
  'text!shaders/reduce_error_s1.frag',
  'text!shaders/reduce_error_s2.frag',
  'text!shaders/update_velocities.frag',
  'text!shaders/update_particles.frag',
  'text!shaders/update_local_bests.frag',
  'text!shaders/update_topological_best_complete.frag',
  'text!shaders/update_topological_best_grid.frag',
  'text!shaders/update_topological_best_ring.frag',
  'text!shaders/hector_fhn.frag',
  'text!shaders/bueno_4v.frag',
  'text!shaders/bueno_brugada.frag',
  'text!shaders/update_global_best.frag',
  'text!shaders/table_ortp.frag',
  'text!shaders/table_ovvr.frag',
  'text!shaders/table_tnnp2006.frag',
  'text!shaders/tnnp2006.frag',
  'text!shaders/ovvr.frag',
  'text!shaders/ortp.frag',
  'text!shaders/aliev-paniflov.frag',
  'text!shaders/set_ics.frag',
  'text!shaders/init_ortp.frag',
  'text!shaders/init_ovvr.frag',
], function(
  GlHelper,
  CopyShader,
  CopyUintShader,
  DefaultVertexShader,
  MitchellSchaefferShader,
  ModifiedMsShader,
  ExpandErrorShader,
  FentonKarmaShader,
  ReduceErrorS1Shader,
  ReduceErrorS2Shader,
  UpdateVelocitiesShader,
  UpdateParticlesShader,
  UpdateLocalBestsShader,
  UpdateTopologicalBestsCompleteShader,
  UpdateTopologicalBestsGridShader,
  UpdateTopologicalBestsRingShader,
  HectorFHNShader,
  Bueno4vShader,
  BuenoBrugadaShader,
  UpdateGlobalBestShader,
  TableOrtpShader,
  TableOvvrShader,
  TableTnnp2006Shader,
  Tnnp2006Shader,
  OvvrShader,
  OrtpShader,
  APShader,
  SetIcsShader,
  InitOrtpShader,
  InitOvvrShader,
) {
  'use strict';

  return class Pso {
    constructor(particle_count) {
      let dim;
      if (particle_count <= 64)
        dim = 8;
      else if (particle_count <= 256)
        dim = 16;
      else if (particle_count <= 1024)
        dim = 32;
      else if (particle_count <= 4096)
        dim = 64;
      else if (particle_count <= 16384)
        dim = 128;
      else
        dim = 256;

      this.particles_width = dim;
      this.particles_height = dim;

      this.tex_width = 2 * dim;
      this.tex_height = 2 * dim;

      const canvas = document.createElement('canvas');
      canvas.width = this.particles_width;
      canvas.height = this.particles_height;

      this.gl_helper = new GlHelper(canvas);
    }

    static data_type_map = {
      'voltage': 0,
      'apd': 1,
      'calcium': 2,
    };

    static err_type_map = {
      'square': 0,
      'abs': 1,
    };

    static model_table_shader_map = {
      'ortp': TableOrtpShader,
      'ovvr': TableOvvrShader,
      'tnnp2006': TableTnnp2006Shader,
    };

    static getEnv() {
      const env = {
        simulation: {
          model: 'fk',
          dt: 0.02,
          period: [],
          num_beats: 1,
          pre_beats: 4,
          v_init: 1.0,
          w_init: 1.0,
          h_init: 1.0,
          align_thresh: [],
          trimmed_data: [],
          data_arrays: [],
          datatypes: [],
          apd_threshs: [],
          weights: [],
          full_normalized_data: [],
          sample_interval: 1.0,
          normalize: true,
          normalization_max: 1.0,
          normalization_min: 0.0,
          normalized_align_threshold: 0.15,
          normalized_ca_align_threshold: 0.15,
          auto_normalize: false,
          err_type: 'abs',
          ics: [],
        },
        stimulus: {
          stim_dur: 10.0,
          stim_mag: 0.4,
          stim_biphasic: true,
          stim_offset_1: 7.0,
          stim_offset_2: 6.72,
          stim_t_scale: 0.725,
        },
        particles: {
          phi_local: 2.05,
          phi_global: 2.05,
          global_bests: new Float32Array(),
          best_error_value: 1e10,
          lower_bounds: [],
          upper_bounds: [],
          learning_rate: 0.05,
          omega: 1.0,
          chi: 0.73,
          parameter_textures: 1,
          topology: 'complete',
        },
        tables: {
          table_width: 512,
          table_height: 512,
          table_shift: 9,
          npoints: 20000,
          vmin: -100.0,
          vmax: 100.0,
          vekmin: -10.0,
          vekmax: 190.0,
        },
        fk_bounds: [
          [25, 10, 50, 0.15, 1, 10, 500, 5, 5, 1, 0.2, 0.05, 0.005],
          [200, 300, 900, 0.4, 20, 50, 1500, 100, 50, 15, 0.9, 0.3, 0.06],
        ],
        ms_bounds: [
          [0.15, 3.0, 75, 60, 0.065],
          [0.6, 12.0, 300, 240, 0.26],
        ],
        mms_bounds: [
          [0.05, 0.5, 60, 60, 0.065],
          [0.6, 12.0, 300, 240, 0.26],
        ],
        fhn_bounds: [
          [0.05, 0.2, 0.001, 0.2, 0.01, -0.1, 0.5],
          [0.6, 2.0, 1.0, 2.0, 1.0, 0.1, 1.5],
        ],
        b4v_bounds: [
          [0.1, 0.5, 1.0, 1.0, 0.01, 15.0, 1.8, 2.0, 1.0, 1.5, 5.0, 5.0, 100.0, 0.05, 5.0, 5.0, 0.1, 0.6, 0.8, 1.0, .1, .005, 0.004, 50.0, .01, 0.4, 1.45],
          [0.35, 300.0, 1500.0, 15.0, 0.04, 100.0, 2.2, 3.0, 20.0, 3.0, 150.0, 150.0, 1000.0, 0.5, 500.0, 10.0, 1.5, 0.7, 1.0, 4.0, 0.15, 0.25, 0.008, 250.0, .2, 1.0, 1.61],
        ],
        bb_bounds: [
          [ // min bounds
            // 1.0,   // tv1p
            // 0.5,   // tv1m
            // 1.0,   // tv2m
            // 5.8024,   // tv1p fixed
            // 60.0,   // tv1m fixed
            // 50.0,   // tv2m fixed

            5.0,   // tv1p
            60.0,   // tv1m fixed
            50.0,   // tv2m fixed
            40,    // tw1p
            150,   // tw2p
            10,    // tw1m
            20,    // tw2m
            5,     // ts1
            50,    // ts2

            // 0.05,  // tfi
            // 5.0,   // to1
            // 5.0,   // to2
            0.05,  // tfi fixed
            400.0,   // to1 fixed
            20.0,   // to2 fixed

            150,   // tso1
            1,     // tso2
            10,    // tsi1
            2,     // tsi2

            // 0.001, // twinf
            // 0.1,   // thv
            // 0.005, // thvm
            // 0.001, // thvinf
            // 0.1,   // thw
            // 0.1,   // thwinf
            // 0.1,   // thso
            // 0.1,   // thsi
            // 0.004, // tho
            // 0.1,   // ths
            0.12, // twinf fixed
            0.13,   // thv fixed
            0.006, // thvm fixed
            2.0, // thvinf fixed
            0.13,   // thw fixed
            0.12,   // thwinf fixed
            0.2,   // thso fixed
            0.13,   // thsi fixed
            0.006, // tho fixed
            0.36,   // ths fixed

            5,     // kwp
            100,   // kwm
            5,     // ks
            1.5,   // kso
            10,    // ksi
            0.02,  // uwm
            0.25,  // us

            // 0.0,   // uo
            // 1.45,  // uu
            0.0,   // uo fixed
            1.0,  // uu fixed

            0.3,   // uso
            0.6,   // sc
            0.2,   // wcp

            // 0.4,   // winfstar
            0.94,   // winfstar fixed
          ],
          [ // max bounds
            // 4.0,    // tv1p
            // 300.0,  // tv1m
            // 1500.0, // tv2m
            10.0,   // tv1p
            60.0,   // tv1m fixed
            50.0,   // tv2m fixed

            80,     // tw1p
            300,    // tw2p
            500,    // tw1m
            40,     // tw2m
            15,     // ts1
            90,     // ts2

            // 0.5,    // tfi
            // 500.0,  // to1
            // 10.0,   // to2
            0.1,  // tfi fixed
            500.0,   // to1 fixed
            35.0,   // to2 fixed

            200,    // tso1
            3,      // tso2
            20,     // tsi1
            10,      // tsi2

            // 0.2,    // twinf
            // 0.35,   // thv
            // 0.25,   // thvm
            // 5.0,    // thvinf
            // 0.15,   // thw
            // 0.2,    // thwinf
            // 0.5,    // thso
            // 0.2,    // thsi
            // 0.008,  // tho
            // 0.5,    // ths
            0.12, // twinf fixed
            0.13,   // thv fixed
            0.006, // thvm fixed
            2.0, // thvinf fixed
            0.13,   // thw fixed
            0.12,   // thwinf fixed
            0.2,   // thso fixed
            0.13,   // thsi fixed
            0.006, // tho fixed
            0.36,   // ths fixed

            10,     // kwp
            150,    // kwm
            25,     // ks
            4,      // kso
            70,     // ksi
            0.12,   // uwm
            0.4,    // us

            // 1.0,    // uo
            // 1.61,   // uu
            0.0,   // uo fixed
            1.0,  // uu fixed

            0.75,   // uso
            0.9,    // sc
            0.3,    // wcp

            // 1.0,    // winfstar
            0.94,   // winfstar fixed
          ],
        ],
        tnnp2006_bounds: [
          // GNa     GK1     Gto     GKr     GKs    GCaL     GpK     GpCa     GbNa     GbCa      pNaK   kNaCa
          [  7.419,  2.7025, 0.0365, 0.0765, 0.196, 1.99e-5, 0.0073, 0.0619,  1.45e-4, 2.96e-4,  1.362, 500.0],
          [  29.676, 10.81,  0.146, 0.306,   0.784, 7.96e-5, 0.0292, 0.2476,  5.8e-4,  0.001184, 5.448, 2000.0],
        ],
        ovvr_bounds: [
          // gnafast gnalate  gto   pca      pcana    pcak       pcacamk  pcanacamk pcakcamk   gkr    gks     gk1     gnaca   gnak  pnab       pcab     gkb     gpca
          [  37.5,   0.00375, 0.01, 0.00005, 6.25e-8, 1.787e-8, 5.5e-5,  6.875e-8, 1.9657e-8, 0.023, 0.0017, 0.0954, 0.0004, 15.0, 1.875e-10, 1.25e-8, 0.0015, 0.00025],
          [  150.0,  0.015,   0.04, 0.0002,  2.5e-7,  7.148e-8, 0.00022, 2.75e-7,  7.8628e-8, 0.092, 0.0068, 0.3816, 0.0016, 60.0, 7.5e-10,   5.0e-8,  0.006,  0.001],

          // [129.40371704101562, 0.007744332309812307, 0.025889599695801735, 0.00005450226672110148, 1.582629209906372e-7, 6.074812830547671e-8, 0.00010538144852034748, 1.4562455419309117e-7, 5.1735945305608766e-8, 0.028930267319083214, 0.0030572263058274984, 0.2827710211277008, 0.001538076321594417, 39.12760925292969, 5.102352984565073e-10, 3.937519821306523e-8, 0.002562799723818898, 0.0005523095023818314],
          // [129.40371704101562, 0.007744332309812307, 0.025889599695801735, 0.00005450226672110148, 1.582629209906372e-7, 6.074812830547671e-8, 0.00010538144852034748, 1.4562455419309117e-7, 5.1735945305608766e-8, 0.028930267319083214, 0.0030572263058274984, 0.2827710211277008, 0.001538076321594417, 39.12760925292969, 5.102352984565073e-10, 3.937519821306523e-8, 0.002562799723818898, 0.0005523095023818314],
        ],
        ortp_bounds: [
          // gna     gto   pca      pcana    pcak       pcacamk  pcanacamk pcakcamk   gkr    gks     gk1     gnaca   gnak  pnab       pcab     gkb     gpca
          [  7.419,  0.01, 0.00005, 6.25e-8, 1.787e-8, 5.5e-5,  6.875e-8, 1.9657e-8, 0.023, 0.0017, 0.0954, 0.0004, 15.0, 1.875e-10, 1.25e-8, 0.0015, 0.00025],
          [  29.676, 0.04, 0.0002,  2.5e-7,  7.148e-8, 0.00022, 2.75e-7,  7.8628e-8, 0.092, 0.0068, 0.3816, 0.0016, 60.0, 7.5e-10,   5.0e-8,  0.006,  0.001],
        ],
        ap_bounds: [
          // k    a      b      eps0     mu1    mu2
          [1.0,   0.01, 0.08, 0.0005,   0.05,   0.05], // min parameter bounds
          [20.0,  0.5,   0.5,   0.2,   0.5,   1] // max parameter bounds
        ],
        velocity_update: {},
        ics: {
          ms: { 0: [0.0, 1.0] },
          mms: { 0: [0.0, 1.0] },
          fhn: { 0: [0.0, 0.0] },
          ap: { 0: [0.0, 0.0] },
          fk: { 0: [0.0, 1.0, 1.0] },
          b4v: { 0: [0.0, 1.0, 1.0, 0.0] },
          bb: { 0: [0.0, 1.0, 1.0, 0.0] },
          tnnp2006: {
            1000: [-84.7, 0.9891, 9.413, 136.1, 0.0001021, 0.002111, 3.385, 0.001634, 0.7512, 0.7508, 0.003213, 3.27e-5, 0.9771, 0.9995, 1.0, 0.0, 0.0, 0.0, 0.0],
            900: [-8.53385313e+01, 9.83510276e-01, 1.01423012e+01, 1.35321308e+02, 1.08655721e-04, 2.30991931e-04, 3.68552855e+00, 1.67633258e-03, 7.47857371e-01, 7.46921348e-01, 3.31005462e-03, 3.32375621e-05, 9.62810937e-01, 9.99495577e-01, 9.99933564e-01, 2.37346651e-08, 5.86365549e-01, 2.12793118e-04, 4.72303646e-01],
            800: [-8.52374482e+01, 9.74470572e-01, 1.06345405e+01, 1.34790260e+02, 1.14568511e-04, 2.55486242e-04, 3.91983411e+00, 1.71270719e-03, 7.45065490e-01, 7.42309821e-01, 3.55867616e-03, 3.36887664e-05, 9.39313134e-01, 9.99488078e-01, 9.99812149e-01, 2.41399037e-08, 5.25452408e-01, 2.52337347e-04, 4.71253088e-01],
            700: [-8.51205174e+01, 9.60013106e-01, 1.11442549e+01, 1.34246204e+02, 1.22005089e-04, 2.89023615e-04, 4.12576672e+00, 1.75574881e-03, 7.41762898e-01, 7.33121483e-01, 4.29903681e-03, 3.42186794e-05, 9.01543453e-01, 9.99479011e-01, 9.99437115e-01, 2.46202923e-08, 4.57376186e-01, 5.67118518e-04, 4.70036737e-01],
            600: [-8.49654090e+01, 9.37389442e-01, 1.16117592e+01, 1.33752701e+02, 1.33295259e-04, 3.36930334e-04, 4.27163886e+00, 1.81448009e-03, 7.37218987e-01, 7.10686687e-01, 6.47963453e-03, 3.49352876e-05, 8.42226950e-01, 9.99466088e-01, 9.98336638e-01, 2.52797892e-08, 3.82995195e-01, 2.91207674e-03, 4.68420382e-01],
            500: [-8.47166025e+01, 9.03217458e-01, 1.19264953e+01, 1.33427036e+02, 1.53304064e-04, 4.07298404e-04, 4.30229870e+00, 1.91272606e-03, 7.29530270e-01, 6.55780869e-01, 1.24172752e-02, 3.61180006e-05, 7.53453969e-01, 9.99443012e-01, 9.95383122e-01, 2.63918208e-08, 3.05087667e-01, 1.71192198e-02, 4.65821656e-01],
            400: [-8.42687847e+01, 8.54058053e-01, 1.19539901e+01, 1.33413331e+02, 1.88802334e-04, 5.11161598e-04, 4.14339144e+00, 2.10281653e-03, 7.14690111e-01, 5.43209978e-01, 2.63307127e-02, 3.83514169e-05, 6.32477747e-01, 9.99375555e-01, 9.88254335e-01, 2.85491065e-08, 2.29430938e-01, 7.82784051e-02, 4.61136893e-01],
            300: [-8.35246206e+01, 7.86468866e-01, 1.18313678e+01, 1.33579840e+02, 2.44381789e-04, 6.66183072e-04, 3.76351121e+00, 2.46020717e-03, 6.87481921e-01, 3.68145464e-01, 5.19800199e-02, 4.23729121e-05, 4.93137820e-01, 9.98186389e-01, 9.71460393e-01, 3.25428810e-08, 1.64576483e-01, 2.39147017e-01, 4.53363048e-01],
            250: [-8.31203885e+01, 7.44617375e-01, 1.15898878e+01, 1.33858314e+02, 2.74475675e-04, 7.72744978e-04, 3.48322700e+00, 2.67846877e-03, 6.64502529e-01, 2.66674833e-01, 6.97759161e-02, 4.47305836e-05, 4.24144592e-01, 9.93550784e-01, 9.54530521e-01, 3.51771074e-08, 1.38173744e-01, 3.49351925e-01, 4.49151085e-01],
          },
          ovvr: {
            1000: [-87.84, 7.23, 7.23, 143.79, 143.79, 8.54e-5, 8.43e-5, 1.61, 1.56, 0.0074621, 0.692591, 0.692574, 0.692477, 0.448501, 0.692413, 0.000194015, 0.496116, 0.265885, 0.00101185, 0.999542, 0.589579, 0.000515567, 0.999542, 0.641861, 2.43015e-9, 1.0, 0.910671, 1.0, 0.99982, 0.999977, 0.00267171, 1.0, 1.0, 8.26608e-6, 0.453268, 0.270492, 0.0001963, 0.996801, 2.53943e-5, 3.17262e-7, 0.0124065],
            900: [-8.78304025e+01, 7.37895397e+00, 7.37904521e+00, 1.43586623e+02, 1.43586594e+02, 8.83052869e-05, 8.74574654e-05, 1.65479771e+00, 1.58371506e+00, 7.47187452e-03, 6.92135133e-01, 6.92109510e-01, 6.91965978e-01, 4.47941013e-01, 6.91868385e-01, 1.94495176e-04, 4.89116500e-01, 2.56510812e-01, 1.01274045e-03, 9.99540563e-01, 5.40736522e-01, 5.16021547e-04, 9.99540574e-01, 5.92843288e-01, 2.43766417e-09, 9.99999990e-01, 8.98268899e-01, 9.99999990e-01, 9.99552689e-01, 9.99917139e-01, 3.07729401e-03, 9.99999990e-01, 9.99999990e-01, 8.71550726e-06, 4.99901222e-01, 2.95185596e-01, 1.96616950e-04, 9.96804702e-01, 2.65009209e-07, 3.31018558e-07, 1.52386152e-02],
            800: [-8.78129934e+01, 7.55661169e+00, 7.55671249e+00, 1.43344500e+02, 1.43344468e+02, 9.25359954e-05, 9.20770019e-05, 1.72462361e+00, 1.61390366e+00, 7.48496658e-03, 6.91525073e-01, 6.91486053e-01, 6.91265319e-01, 4.47183833e-01, 6.91110618e-01, 1.95139414e-04, 4.78279329e-01, 2.44938505e-01, 1.01393340e-03, 9.99539115e-01, 4.86420034e-01, 5.16629693e-04, 9.99539133e-01, 5.37332408e-01, 2.44774394e-09, 9.99999990e-01, 8.81952630e-01, 9.99999990e-01, 9.98876518e-01, 9.99709765e-01, 3.74842708e-03, 9.99999990e-01, 9.99999990e-01, 1.23597740e-05, 5.51615861e-01, 3.23644829e-01, 1.97050250e-04, 9.96809502e-01, 2.79998184e-07, 3.49638745e-07, 1.96074668e-02],
            700: [-8.77905604e+01, 7.75743441e+00, 7.75754832e+00, 1.43066545e+02, 1.43066511e+02, 9.82747271e-05, 9.84440877e-05, 1.82579478e+00, 1.65025601e+00, 7.50187066e-03, 6.90737914e-01, 6.90677101e-01, 6.90327810e-01, 4.46191210e-01, 6.90059942e-01, 1.95972729e-04, 4.61627849e-01, 2.30468752e-01, 1.01547393e-03, 9.99537228e-01, 4.26506120e-01, 5.17415029e-04, 9.99537255e-01, 4.75085303e-01, 2.46080275e-09, 9.99999990e-01, 8.61483839e-01, 9.99999990e-01, 9.97219343e-01, 9.98999541e-01, 4.84159553e-03, 9.99999990e-01, 9.99999990e-01, 4.18597416e-05, 6.07838540e-01, 3.56681843e-01, 1.97637923e-04, 9.96815775e-01, 2.96881473e-07, 3.70570493e-07, 2.65867423e-02, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
            600: [-8.77602382e+01, 7.98792788e+00, 7.98806029e+00, 1.42738468e+02, 1.42738430e+02, 1.06098411e-04, 1.07350198e-04, 1.98609653e+00, 1.69951146e+00, 7.52478060e-03, 6.89672051e-01, 6.89573943e-01, 6.88994059e-01, 4.44819805e-01, 6.88398518e-01, 1.97104801e-04, 4.35910105e-01, 2.11925319e-01, 1.01756209e-03, 9.99534637e-01, 3.60605000e-01, 5.18479545e-04, 9.99534681e-01, 4.05403709e-01, 2.47858004e-09, 9.99999990e-01, 8.34940031e-01, 9.99999990e-01, 9.93227278e-01, 9.96609816e-01, 6.74423977e-03, 9.99999989e-01, 9.99999990e-01, 2.74082340e-04, 6.68712209e-01, 3.95978532e-01, 1.98787047e-04, 9.96824407e-01, 3.17809718e-07, 3.96478498e-07, 3.91387648e-02],
            500: [-8.77165275e+01, 8.25716945e+00, 8.25732909e+00, 1.42335600e+02, 1.42335559e+02, 1.16940516e-04, 1.20244320e-04, 2.26395658e+00, 1.77315913e+00, 7.55792977e-03, 6.88131807e-01, 6.87964921e-01, 6.86862702e-01, 4.42784945e-01, 6.84726394e-01, 1.98748312e-04, 3.95707373e-01, 1.87272453e-01, 1.02058387e-03, 9.99530825e-01, 2.88332411e-01, 5.20020002e-04, 9.99530901e-01, 3.27460739e-01, 2.50446160e-09, 9.99999990e-01, 7.98948640e-01, 9.99999990e-01, 9.83774109e-01, 9.88708530e-01, 1.04293616e-02, 9.99999870e-01, 9.99999982e-01, 2.03682229e-03, 7.34234058e-01, 4.44268707e-01, 2.10145392e-04, 9.96837144e-01, 3.45063699e-07, 4.30195743e-07, 6.55670866e-02],
            400: [-8.76389766e+01, 8.58061070e+00, 8.58081763e+00, 1.41821157e+02, 1.41821111e+02, 1.34922275e-04, 1.42589220e-04, 2.75226738e+00, 1.86536335e+00, 7.61710198e-03, 6.85388523e-01, 6.85075389e-01, 6.81286479e-01, 4.39073885e-01, 6.70224516e-01, 2.01698115e-04, 3.31328347e-01, 1.52696134e-01, 1.02597340e-03, 9.99523901e-01, 2.09619546e-01, 5.22767526e-04, 9.99524046e-01, 2.40584407e-01, 2.55109181e-09, 9.99999990e-01, 7.47387183e-01, 9.99999990e-01, 9.61815700e-01, 9.63163878e-01, 2.01873938e-02, 9.99984065e-01, 9.99998612e-01, 1.47565924e-02, 8.03913502e-01, 5.05473252e-01, 4.92209768e-04, 9.96860164e-01, 3.73449614e-07, 4.65162017e-07, 1.27474405e-01],
            300: [-8.74415799e+01, 8.98030224e+00, 8.98062158e+00, 1.41180230e+02, 1.41180178e+02, 1.78659941e-04, 1.97567828e-04, 3.47116021e+00, 1.86056193e+00, 7.76980193e-03, 6.78347190e-01, 6.77649508e-01, 6.41357432e-01, 4.27747481e-01, 5.89446973e-01, 2.09405596e-04, 2.25115138e-01, 1.00931933e-01, 1.03982186e-03, 9.99505606e-01, 1.26142700e-01, 5.29827380e-04, 9.99506078e-01, 1.45899030e-01, 2.67374510e-09, 9.99999865e-01, 6.69845077e-01, 9.99999989e-01, 9.13204106e-01, 8.85440130e-01, 7.00516497e-02, 9.98293507e-01, 9.99812779e-01, 9.73264951e-02, 8.75025898e-01, 5.82252543e-01, 6.95297397e-03, 9.96924228e-01, 3.76055953e-07, 4.66513698e-07, 2.46966291e-01],
          },
          ortp: {
            1000: [-87.852948063675726, 7.0392546150493178, 7.0393347733953027, 143.96135641252130, 143.96133202798799, 8.3881023438298657e-5, 8.2772564492747630e-5, 1.5672795239503670, 1.5261659336262898, 9.7991163152215356e-4, 0.80949981914040681, 0.80934296248967219, 1.0111996752124317e-3, 0.99954240640082770, 0.60310052987525253, 5.1523608441461936e-4, 0.99954241327569771, 0.65841214647941604, 2.4246864435205784e-9, 0.99999999049967647, 0.92312278514718482, 0.99999999049973498, 0.99987265676003101, 0.99998338629400885, 2.4911629583080604e-3, 0.99999999049692279, 0.99999999049708876, 8.2273738196810401e-6, 0.43975667031918381, 0.24476268201089038, 1.9608490236542564e-4, 0.99679863494646970, 2.3426825206602286e-7, 2.9268612848908905e-7, 1.1108986413883167e-2],
          },
        },
      };

      return env;
    }

    setupEnv(model, bounds, stimulus_params, pre_beats, num_beats, sample_interval, normalize, normalization_max, normalization_min, auto_normalize, hyperparams) {
      this.env = Pso.getEnv();
      const env = this.env;

      env.simulation.model = model;

      if (Number(num_beats)) {
        env.simulation.num_beats = Number(num_beats);
      }

      if (!Number.isNaN(Number(pre_beats))) {
        env.simulation.pre_beats = Number(pre_beats);
      }

      if (Number(sample_interval)) {
        env.simulation.sample_interval = Number(sample_interval);
      }

      env.simulation.auto_normalize = auto_normalize;
      if (auto_normalize) {
        normalize = true;
        normalization_max = 1;
        normalization_min = 0;
      }

      env.simulation.normalize = normalize;
      env.simulation.normalization_max = normalization_max;
      env.simulation.normalization_min = normalization_min;

      env.stimulus = stimulus_params;

      // 16 parameters per texture for now
      env.particles.parameter_textures = Math.ceil(bounds[0].length/16);

      env.particles.phi_local = hyperparams.phi1;
      env.particles.phi_global = hyperparams.phi2;
      // TODO Elizabeth's Brugada code scales the standard chi value by 0.25, which is worth investigating
      env.particles.chi = hyperparams.chi;

      env.particles.lower_bounds = bounds[0];
      env.particles.upper_bounds = bounds[1];

      // Pad out these arrays so chunks of 16 can always be used as uniforms
      while (env.particles.lower_bounds.length % 16 !== 0) {
        env.particles.lower_bounds.push(0);
        env.particles.upper_bounds.push(0);
      }

      // Initialize all global best values to 0
      env.particles.global_bests = new Float32Array(env.particles.lower_bounds.length);

      env.particles.particle_count = this.particles_width * this.particles_height;
      env.particles.iteration_count = hyperparams.iteration_count;
    }

    normalizeData(parsed_data, normalization_max, normalization_min) {
      const min = Math.min(...parsed_data);
      const max = Math.max(...parsed_data);

      const normalized_data = parsed_data.map(x => normalization_min + (x-min)*((normalization_max-normalization_min)/(max-min)));

      return normalized_data;
    }

    readData(input_data) {
      const raw_input_data = [];
      const input_cls = [];
      const datatypes = [];
      const apd_threshs = [];
      const weights = [];
      for (const obj of input_data) {
        raw_input_data.push(obj.data);
        input_cls.push(obj.cl);
        datatypes.push(obj.datatype);
        apd_threshs.push(obj.apd_thresh || 0);
        weights.push(obj.weight);
      }

      this.env.simulation.period = input_cls;
      this.env.simulation.datatypes = datatypes;
      this.env.simulation.apd_threshs = apd_threshs;

      const trimmed_data = [];
      const data_arrays = [];
      const align_thresh = [];
      const all_full_normalized_data = [];
      this.env.simulation.ics = [];

      for (let i = 0; i < raw_input_data.length; ++i) {
        if (datatypes[i] === 'apd') {
          const apd_data = raw_input_data[i];

          const data_array = new Float32Array(4 * apd_data.length);
          for (let j = 0; j < apd_data.length; ++j) {
            data_array[4*j] = apd_data[j];
          }

          trimmed_data.push(apd_data);
          data_arrays.push(data_array);
          align_thresh.push(0);
          all_full_normalized_data.push(apd_data);
        } else {
          const delta = datatypes[i] === 'calcium' ? 1e-7 : 0.001;
          const nthresh = datatypes[i] === 'calcium' ? this.env.simulation.normalized_ca_align_threshold : this.env.simulation.normalized_align_threshold;
          const raw_text = raw_input_data[i];

          const split_data = raw_text.split('\n');
          const actual_data = split_data.filter(x => !(x.trim() === ""));

          const full_parsed_data = actual_data.map(x => parseFloat(x.trim()));
          const full_normalized_data = this.env.simulation.normalize ? this.normalizeData(full_parsed_data, this.env.simulation.normalization_max, this.env.simulation.normalization_min) : full_parsed_data;

          const data_max = Math.max(...full_normalized_data);
          const data_min = Math.min(...full_normalized_data);
          const actual_align_thresh = data_min + nthresh * (data_max-data_min);
          const first_compare_index = full_normalized_data.findIndex(number => number > actual_align_thresh);
          const last_compare_index = this.env.simulation.num_beats * input_cls[i];
          const curr_trimmed_data = full_normalized_data.slice(first_compare_index, last_compare_index);

          // Pad out the extra pixel values. The data could be stored more densely by using the full pixel
          // value and by using a two-dimensional texture, but for now there is not enough to require that.
          const data_length = curr_trimmed_data.length;
          const data_array = new Float32Array(4 * data_length);
          for (let j = 0; j < data_length; ++j) {
            data_array[4*j] = curr_trimmed_data[j];
          }

          trimmed_data.push(curr_trimmed_data);
          data_arrays.push(data_array);
          align_thresh.push(curr_trimmed_data[0] - delta);
          all_full_normalized_data.push(full_normalized_data);
        }

        // Set up the initial conditions for the data
        const ics_obj = this.env.ics[this.env.simulation.model];
        const ics_cls = Object.keys(ics_obj);
        ics_cls.sort((a, b) => a - b);
        let closest_cl = ics_cls.find(cl => input_cls[i] <= cl);
        if (closest_cl === undefined) {
          closest_cl = ics_cls[ics_cls.length-1];
        }

        const ics = ics_obj[closest_cl];

        // Pad out the array so chunks of 4 can always be used as uniforms
        while (ics.length % 4 !== 0) {
          ics.push(0);
        }

        this.env.simulation.ics.push(ics);
      }

      this.env.simulation.data_arrays = data_arrays;
      this.env.simulation.trimmed_data = trimmed_data;
      this.env.simulation.align_thresh = align_thresh;
      this.env.simulation.full_normalized_data = all_full_normalized_data;
      this.env.simulation.weights = weights;
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
      const particle_zero_array = new Float32Array(particles_width*particles_height*4);

      const global_best_array = new Float32Array(16);

      this.simulation_lengths = [];
      this.data_textures = [];

      for (let i = 0; i < period.length; i++) {
        this.simulation_lengths.push(Math.ceil(Math.ceil(num_beats * period[i]) / sample_interval));
        this.data_textures.push(gl_helper.loadFloatTexture(data_arrays[i].length/4, 1, data_arrays[i]));
      }

      const sim_zero_array = new Float32Array(Math.max(...this.simulation_lengths)*4);

      this.particles_textures = [];
      this.velocities_textures = [];
      this.bests_textures = [];
      this.global_best_textures = [];
      this.particles_out_textures = [];
      this.velocities_out_textures = [];
      this.bests_out_textures = [];
      this.global_best_out_textures = [];

      for (const init_array of init_arrays) {
        this.particles_textures.push(gl_helper.loadFloatTexture(tex_width, tex_height, init_array));
        this.velocities_textures.push(gl_helper.loadFloatTexture(tex_width, tex_height, zero_array));
        this.bests_textures.push(gl_helper.loadFloatTexture(tex_width, tex_height, zero_array));
        this.global_best_textures.push(gl_helper.loadFloatTexture(2, 2, global_best_array));
        this.particles_out_textures.push(gl_helper.loadFloatTexture(tex_width, tex_height, null));
        this.velocities_out_textures.push(gl_helper.loadFloatTexture(tex_width, tex_height, null));
        this.bests_out_textures.push(gl_helper.loadFloatTexture(tex_width, tex_height, null));
        this.global_best_out_textures.push(gl_helper.loadFloatTexture(2, 2, null));
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

      this.reduced_error_1_texture = gl_helper.loadFloatTexture(particles_width, 1, null);
      this.reduced_error_2_texture = gl_helper.loadFloatTexture(1, 1, null);

      this.topological_best_idx_texture = gl_helper.loadUintTexture(particles_width, particles_height, null);

      // These need to be 2x2 to match the global best texture
      const best_error_value_array = new Float32Array(16);
      best_error_value_array[0] = this.env.particles.best_error_value;
      best_error_value_array[4] = this.env.particles.best_error_value;
      best_error_value_array[8] = this.env.particles.best_error_value;
      best_error_value_array[12] = this.env.particles.best_error_value;
      this.best_error_value_texture = gl_helper.loadFloatTexture(2, 2, best_error_value_array);
      this.best_error_value_out_texture = gl_helper.loadFloatTexture(2, 2, null);

      // Table textures
      this.table_texture = gl_helper.loadFloatTexture(this.env.tables.table_width, this.env.tables.table_height, null);

      this.state_textures = [];
      this.final_state_textures = [];
      this.state_out_textures = [];
      this.final_state_out_textures = [];

      let num_state_textures;
      if (this.env.simulation.model === 'ovvr') {
        num_state_textures = 8;
      } else if (this.env.simulation.model === 'ortp') {
        num_state_textures = 7;
      } else {
        num_state_textures = this.env.simulation.ics[0].length/4;
      }

      for (let i = 0; i < Math.floor(num_state_textures); ++i) {
        this.state_textures.push([]);
        this.state_out_textures.push([]);

        this.final_state_textures.push([]);
        this.final_state_out_textures.push([]);

        for (let cl = 0; cl < period.length; ++cl) {
          this.state_textures[i].push(gl_helper.loadFloatTexture(particles_width, particles_height, null));
          this.state_out_textures[i].push(gl_helper.loadFloatTexture(particles_width, particles_height, null));

          this.final_state_textures[i].push(gl_helper.loadFloatTexture(Math.max(...this.simulation_lengths), 1, null));
          this.final_state_out_textures[i].push(gl_helper.loadFloatTexture(Math.max(...this.simulation_lengths), 1, null));
        }
      }

      this.normalize_textures = [];
      this.final_normalize_textures = [];
      for (let cl = 0; cl < period.length; ++cl) {
        this.normalize_textures.push(gl_helper.loadFloatTexture(particles_width, particles_height, particle_zero_array));
        this.final_normalize_textures.push(gl_helper.loadFloatTexture(Math.max(...this.simulation_lengths), 1, sim_zero_array));
      }
      this.dummy_normalize_texture = gl_helper.loadFloatTexture(1, 1, new Float32Array([0, 0, 0, 0]));

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
          run: this.gl_helper.runProgram,
          dims: [this.tex_width, this.tex_height],
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
            ['topological_best_idx_texture', 'tex', () => this.topological_best_idx_texture],
            ['itinymtState', 'tex', () => this.env.velocity_update.ftinymtState],
            ['itinymtMat', 'tex', () => this.env.velocity_update.itinymtMat],
            ['phi_local', '1f', () => this.env.particles.phi_local],
            ['phi_global', '1f', () => this.env.particles.phi_global],
            ['omega', '1f', () => this.env.particles.omega],
            ['chi', '1f', () => this.env.particles.chi],
          ],
          out: [this.velocities_out_textures[num], this.env.velocity_update.stinymtState],
          run: this.gl_helper.runProgram,
          dims: [this.tex_width, this.tex_height],
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
          run: this.gl_helper.runProgram,
          dims: [this.tex_width, this.tex_height],
        };
      };

      const makeGlobalBestUpdateSolver = (num) => {
        return {
          vert: DefaultVertexShader,
          frag: UpdateGlobalBestShader,
          uniforms: [
            ['positions_texture', 'tex', () => this.particles_textures[num]],
            ['reduced_error_2_texture', 'tex', () => this.reduced_error_2_texture],
            ['best_error_value_texture', 'tex', () => this.best_error_value_texture],
            ['global_best_texture', 'tex', () => this.global_best_textures[num]],
          ],
          out: [this.global_best_out_textures[num], this.best_error_value_out_texture],
          run: this.gl_helper.runProgram,
          dims: [2, 2],
        };
      };

      const makeCopySolver = (original, copy, idx, width, height) => {
        return {
          vert: DefaultVertexShader,
          frag: CopyShader,
          uniforms: [
            ['original', 'tex', idx === undefined ? () => this[original] : () => this[original][idx]],
          ],
          out: [idx === undefined ? this[copy] : this[copy][idx]],
          run: this.gl_helper.runProgram,
          dims: [width || this.tex_width, height || this.tex_height],
        };
      };

      const makeInitStateSolver = (num, cl_idx, final) => {
        return {
          vert: DefaultVertexShader,
          frag: SetIcsShader,
          uniforms: [
            ['ics', '4fv_a', () => [this.env.simulation.ics[cl_idx], num*4, 4]],
          ],
          out: [final ? this.final_state_textures[num][cl_idx] : this.state_textures[num][cl_idx]],
          run: this.gl_helper.runProgram,
          dims: [final ? Math.max(...this.simulation_lengths) : this.particles_width, final ? 1 : this.particles_height],
        }
      };

      const makeInitOrtpSolver = (cl_idx, final) => {
        const solver = {
          vert: DefaultVertexShader,
          frag: InitOrtpShader,
          uniforms: [
            ['V',           '1f', () => this.env.simulation.ics[cl_idx][0]],
            ['Na_i',        '1f', () => this.env.simulation.ics[cl_idx][1]],
            ['Na_ss',       '1f', () => this.env.simulation.ics[cl_idx][2]],
            ['K_i',         '1f', () => this.env.simulation.ics[cl_idx][3]],
            ['K_ss',        '1f', () => this.env.simulation.ics[cl_idx][4]],
            ['Ca_i',        '1f', () => this.env.simulation.ics[cl_idx][5]],
            ['Ca_ss',       '1f', () => this.env.simulation.ics[cl_idx][6]],
            ['Ca_nsr',      '1f', () => this.env.simulation.ics[cl_idx][7]],
            ['Ca_jsr',      '1f', () => this.env.simulation.ics[cl_idx][8]],
            ['m',           '1f', () => this.env.simulation.ics[cl_idx][9]],
            ['h',           '1f', () => this.env.simulation.ics[cl_idx][10]],
            ['j',           '1f', () => this.env.simulation.ics[cl_idx][11]],
            ['a',           '1f', () => this.env.simulation.ics[cl_idx][12]],
            ['ifast',       '1f', () => this.env.simulation.ics[cl_idx][13]],
            ['islow',       '1f', () => this.env.simulation.ics[cl_idx][14]],
            ['aCaMK',       '1f', () => this.env.simulation.ics[cl_idx][15]],
            ['iCaMKfast',   '1f', () => this.env.simulation.ics[cl_idx][16]],
            ['iCaMKslow',   '1f', () => this.env.simulation.ics[cl_idx][17]],
            ['d',           '1f', () => this.env.simulation.ics[cl_idx][18]],
            ['ffast',       '1f', () => this.env.simulation.ics[cl_idx][19]],
            ['fslow',       '1f', () => this.env.simulation.ics[cl_idx][20]],
            ['fCafast',     '1f', () => this.env.simulation.ics[cl_idx][21]],
            ['fCaslow',     '1f', () => this.env.simulation.ics[cl_idx][22]],
            ['jCa',         '1f', () => this.env.simulation.ics[cl_idx][23]],
            ['n',           '1f', () => this.env.simulation.ics[cl_idx][24]],
            ['fCaMKfast',   '1f', () => this.env.simulation.ics[cl_idx][25]],
            ['fCaCaMKfast', '1f', () => this.env.simulation.ics[cl_idx][26]],
            ['xrfast',      '1f', () => this.env.simulation.ics[cl_idx][27]],
            ['xrslow',      '1f', () => this.env.simulation.ics[cl_idx][28]],
            ['xs1',         '1f', () => this.env.simulation.ics[cl_idx][29]],
            ['xs2',         '1f', () => this.env.simulation.ics[cl_idx][30]],
            ['xK1',         '1f', () => this.env.simulation.ics[cl_idx][31]],
            ['JrelNP',      '1f', () => this.env.simulation.ics[cl_idx][32]],
            ['JrelCaMK',    '1f', () => this.env.simulation.ics[cl_idx][33]],
            ['CaMKtrap',    '1f', () => this.env.simulation.ics[cl_idx][34]],
          ],
          out: [],
          run: this.gl_helper.runProgram,
          dims: [final ? Math.max(...this.simulation_lengths) : this.particles_width, final ? 1 : this.particles_height],
        }

        const out_textures = final ? this.final_state_textures : this.state_textures;

        for (let i = 0; i < out_textures.length; ++i) {
          solver.out.push(out_textures[i][cl_idx]);
        }

        return solver;
      };

      const makeInitOvvrSolver = (cl_idx, final) => {
        const solver = {
          vert: DefaultVertexShader,
          frag: InitOvvrShader,
          uniforms: [
            ['V',           '1f', () => this.env.simulation.ics[cl_idx][0]],
            ['Na_i',        '1f', () => this.env.simulation.ics[cl_idx][1]],
            ['Na_ss',       '1f', () => this.env.simulation.ics[cl_idx][2]],
            ['K_i',         '1f', () => this.env.simulation.ics[cl_idx][3]],
            ['K_ss',        '1f', () => this.env.simulation.ics[cl_idx][4]],
            ['Ca_i',        '1f', () => this.env.simulation.ics[cl_idx][5]],
            ['Ca_ss',       '1f', () => this.env.simulation.ics[cl_idx][6]],
            ['Ca_nsr',      '1f', () => this.env.simulation.ics[cl_idx][7]],
            ['Ca_jsr',      '1f', () => this.env.simulation.ics[cl_idx][8]],
            ['m',           '1f', () => this.env.simulation.ics[cl_idx][9]],
            ['hfast',       '1f', () => this.env.simulation.ics[cl_idx][10]],
            ['hslow',       '1f', () => this.env.simulation.ics[cl_idx][11]],
            ['j',           '1f', () => this.env.simulation.ics[cl_idx][12]],
            ['hCaMKslow',   '1f', () => this.env.simulation.ics[cl_idx][13]],
            ['jCaMK',       '1f', () => this.env.simulation.ics[cl_idx][14]],
            ['mL',          '1f', () => this.env.simulation.ics[cl_idx][15]],
            ['hL',          '1f', () => this.env.simulation.ics[cl_idx][16]],
            ['hLCaMK',      '1f', () => this.env.simulation.ics[cl_idx][17]],
            ['a',           '1f', () => this.env.simulation.ics[cl_idx][18]],
            ['ifast',       '1f', () => this.env.simulation.ics[cl_idx][19]],
            ['islow',       '1f', () => this.env.simulation.ics[cl_idx][20]],
            ['aCaMK',       '1f', () => this.env.simulation.ics[cl_idx][21]],
            ['iCaMKfast',   '1f', () => this.env.simulation.ics[cl_idx][22]],
            ['iCaMKslow',   '1f', () => this.env.simulation.ics[cl_idx][23]],
            ['d',           '1f', () => this.env.simulation.ics[cl_idx][24]],
            ['ffast',       '1f', () => this.env.simulation.ics[cl_idx][25]],
            ['fslow',       '1f', () => this.env.simulation.ics[cl_idx][26]],
            ['fCafast',     '1f', () => this.env.simulation.ics[cl_idx][27]],
            ['fCaslow',     '1f', () => this.env.simulation.ics[cl_idx][28]],
            ['jCa',         '1f', () => this.env.simulation.ics[cl_idx][29]],
            ['n',           '1f', () => this.env.simulation.ics[cl_idx][30]],
            ['fCaMKfast',   '1f', () => this.env.simulation.ics[cl_idx][31]],
            ['fCaCaMKfast', '1f', () => this.env.simulation.ics[cl_idx][32]],
            ['xrfast',      '1f', () => this.env.simulation.ics[cl_idx][33]],
            ['xrslow',      '1f', () => this.env.simulation.ics[cl_idx][34]],
            ['xs1',         '1f', () => this.env.simulation.ics[cl_idx][35]],
            ['xs2',         '1f', () => this.env.simulation.ics[cl_idx][36]],
            ['xK1',         '1f', () => this.env.simulation.ics[cl_idx][37]],
            ['JrelNP',      '1f', () => this.env.simulation.ics[cl_idx][38]],
            ['JrelCaMK',    '1f', () => this.env.simulation.ics[cl_idx][39]],
            ['CaMKtrap',    '1f', () => this.env.simulation.ics[cl_idx][40]],
          ],
          out: [],
          run: this.gl_helper.runProgram,
          dims: [final ? Math.max(...this.simulation_lengths) : this.particles_width, final ? 1 : this.particles_height],
        }

        const out_textures = final ? this.final_state_textures : this.state_textures;

        for (let i = 0; i < out_textures.length; ++i) {
          solver.out.push(out_textures[i][cl_idx]);
        }

        return solver;
      };

      const model_shader_map = {
        'fk': FentonKarmaShader,
        'ms': MitchellSchaefferShader,
        'mms': ModifiedMsShader,
        'fhn': HectorFHNShader,
        'b4v': Bueno4vShader,
        'bb': BuenoBrugadaShader,
        'tnnp2006': Tnnp2006Shader,
        'ovvr': OvvrShader,
        'ortp': OrtpShader,
        'ap': APShader,
      };

      const makeRunSimulationSolver = (final, prepace, normalize) => {
        const solver = {
          vert: DefaultVertexShader,
          frag: model_shader_map[this.env.simulation.model],
          uniforms: [
            ['data_texture', 'tex', (cl_idx) => this.data_textures[cl_idx]],
            ['dt', '1f', () => this.env.simulation.dt],
            ['period', '1f', (cl_idx) => this.env.simulation.period[cl_idx]],
            ['stim_dur', '1f', () => this.env.stimulus.stim_dur],
            ['stim_mag', '1f', () => this.env.stimulus.stim_mag],
            ['stim_biphasic', '1i', () => this.env.stimulus.stim_biphasic],
            ['stim_offset_1', '1f', () => this.env.stimulus.stim_offset_1],
            ['stim_offset_2', '1f', () => this.env.stimulus.stim_offset_2],
            ['stim_t_scale', '1f', () => this.env.stimulus.stim_t_scale],
            ['num_beats', '1i', () => this.env.simulation.num_beats],
            ['prepacing', '1i', () => prepace],
            ['align_thresh', '1f', (cl_idx) => this.env.simulation.align_thresh[cl_idx]],
            ['sample_interval', '1f', () => this.env.simulation.sample_interval],
            ['data_type', '1i', (cl_idx) => Pso.data_type_map[this.env.simulation.datatypes[cl_idx]]],
            ['apd_thresh', '1f', (cl_idx) => this.env.simulation.apd_threshs[cl_idx]],
            ['weight', '1f', (cl_idx) => this.env.simulation.weights[cl_idx]],
            ['err_type', '1i', () => Pso.err_type_map[this.env.simulation.err_type]],
            ['normalizing', '1i', () => normalize],
            ['auto_normalize', '1i', () => this.env.simulation.auto_normalize],
          ],
          out: [],
          run: final ? this.gl_helper.runFinal : this.gl_helper.runSimulation,
        };

        for (let i = 0; i < this.particles_textures.length; ++i) {
          if (final) {
            solver.uniforms.push(['in_particles_' + (i+1), 'tex', () => this.final_particles_textures[i]]);
          } else {
            solver.uniforms.push(['in_particles_' + (i+1), 'tex', () => this.particles_textures[i]]);
          }
        }

        if (prepace) {
          if (final) {
            solver.out.push((cl_idx) => this.final_state_out_textures[0][cl_idx]);
          } else {
            solver.out.push((cl_idx) => this.state_out_textures[0][cl_idx]);
          }
        } else if (normalize) {
          if (final) {
            solver.out.push((cl_idx) => this.final_normalize_textures[cl_idx]);
          } else {
            solver.out.push((cl_idx) => this.normalize_textures[cl_idx]);
          }
        } else {
          if (final) {
            solver.out.push(this.simulation_texture);
          } else {
            solver.out.push(this.error_texture);
          }
        }

        for (let i = 0; i < this.state_textures.length; ++i) {
          if (final) {
            solver.uniforms.push(['state_textures_' + i, 'tex', (cl_idx) => this.final_state_textures[i][cl_idx]]);
            if (i > 0) {
              solver.out.push((cl_idx) => this.final_state_out_textures[i][cl_idx]);
            }
          } else {
            solver.uniforms.push(['state_textures_' + i, 'tex', (cl_idx) => this.state_textures[i][cl_idx]]);
            if (i > 0) {
              solver.out.push((cl_idx) => this.state_out_textures[i][cl_idx]);
            }
          }
        }

        if (normalize) {
          solver.uniforms.push(['normalize_texture', 'tex', () => this.dummy_normalize_texture]);
        } else {
          if (final) {
            solver.uniforms.push(['normalize_texture', 'tex', (cl_idx) => this.final_normalize_textures[cl_idx]]);
          } else {
            solver.uniforms.push(['normalize_texture', 'tex', (cl_idx) => this.normalize_textures[cl_idx]]);
          }
        }

        if (Object.keys(Pso.model_table_shader_map).includes(this.env.simulation.model)) {
          solver.uniforms.push(
            ['table', 'tex', () => this.table_texture],
            ['table_shift', '1i', () => this.env.tables.table_shift],
            ['table_npoints', '1i', () => this.env.tables.npoints],
            ['table_vmin', '1f', () => this.env.tables.vmin],
            ['table_vmax', '1f', () => this.env.tables.vmax],
            ['table_vekmin', '1f', () => this.env.tables.vekmin],
            ['table_vekmax', '1f', () => this.env.tables.vekmax],
          );
        }

        return solver;
      };

      const shader_map = {
        run_simulation: makeRunSimulationSolver(false, false, false),
        run_final_simulation: makeRunSimulationSolver(true, false, false),

        run_simulation_prepace: makeRunSimulationSolver(false, true, false),
        run_final_simulation_prepace: makeRunSimulationSolver(true, true, false),

        run_simulation_normalize: makeRunSimulationSolver(false, false, true),
        run_final_simulation_normalize: makeRunSimulationSolver(true, false, true),

        reduce_error_1: {
          vert: DefaultVertexShader,
          frag: ReduceErrorS1Shader,
          uniforms: [
            ['error_texture', 'tex', () => this.error_texture],
          ],
          out: [this.reduced_error_1_texture],
          run: this.gl_helper.runProgram,
          dims: [this.particles_width, 1],
        },

        reduce_error_2: {
          vert: DefaultVertexShader,
          frag: ReduceErrorS2Shader,
          uniforms: [
            ['reduced_error_1', 'tex', () => this.reduced_error_1_texture],
          ],
          out: [this.reduced_error_2_texture],
          run: this.gl_helper.runProgram,
          dims: [1, 1],
        },

        expand_error: {
          vert: DefaultVertexShader,
          frag: ExpandErrorShader,
          uniforms: [
            ['error_texture', 'tex', () => this.error_texture],
          ],
          out: [this.expanded_error_texture],
          run: this.gl_helper.runProgram,
          dims: [this.tex_width, this.tex_height],
        },

        tinymt_copy: {
          vert: DefaultVertexShader,
          frag: CopyUintShader,
          uniforms: [
            ['original', 'tex', () => this.env.velocity_update.stinymtState],
          ],
          out: [this.env.velocity_update.ftinymtState],
          run: this.gl_helper.runProgram,
          dims: [this.tex_width, this.tex_height],
        },

        local_error_copy: makeCopySolver('local_bests_error_texture_out', 'local_bests_error_texture'),
        best_error_value_copy: makeCopySolver('best_error_value_out_texture', 'best_error_value_texture', undefined, 2, 2),

        update_topological_best_complete: {
          vert: DefaultVertexShader,
          frag: UpdateTopologicalBestsCompleteShader,
          uniforms: [
            ['best_error_value_texture', 'tex', () => this.best_error_value_texture],
          ],
          out: [this.topological_best_idx_texture],
          run: this.gl_helper.runProgram,
          dims: [this.particles_width, this.particles_height],
        },

        update_topological_best_ring: {
          vert: DefaultVertexShader,
          frag: UpdateTopologicalBestsRingShader,
          uniforms: [
            ['local_bests_error_texture', 'tex', () => this.local_bests_error_texture],
          ],
          out: [this.topological_best_idx_texture],
          run: this.gl_helper.runProgram,
          dims: [this.particles_width, this.particles_height],
        },

        update_topological_best_grid: {
          vert: DefaultVertexShader,
          frag: UpdateTopologicalBestsGridShader,
          uniforms: [
            ['local_bests_error_texture', 'tex', () => this.local_bests_error_texture],
          ],
          out: [this.topological_best_idx_texture],
          run: this.gl_helper.runProgram,
          dims: [this.particles_width, this.particles_height],
        },
      };

      if (Object.keys(Pso.model_table_shader_map).includes(this.env.simulation.model)) {
        shader_map.table = {
          vert: DefaultVertexShader,
          frag: Pso.model_table_shader_map[this.env.simulation.model],
          uniforms: [
            ['width', '1i', () => this.env.tables.table_width],
            ['height', '1i', () => this.env.tables.table_height],
            ['npoints', '1i', () => this.env.tables.npoints],
            ['vmin', '1f', () => this.env.tables.vmin],
            ['vmax', '1f', () => this.env.tables.vmax],
            ['vekmin', '1f', () => this.env.tables.vekmin],
            ['vekmax', '1f', () => this.env.tables.vekmax],
            ['dt', '1f', () => this.env.simulation.dt],
          ],
          out: [this.table_texture],
          run: this.gl_helper.runProgram,
          dims: [this.env.tables.table_width, this.env.tables.table_height],
        };
      }

      for (let i = 0; i < this.env.particles.parameter_textures; ++i) {
        shader_map['local_best_update_' + i] = makeUpdateLocalBestsSolver(i);
        shader_map['velocity_' + i] = makeVelocityUpdateSolver(i);
        shader_map['position_' + i] = makeParticleUpdateSolver(i);
        shader_map['local_bests_copy_' + i] = makeCopySolver('bests_out_textures', 'bests_textures', i);
        shader_map['positions_copy_' + i] = makeCopySolver('particles_out_textures', 'particles_textures', i);
        shader_map['velocities_copy_' + i] = makeCopySolver('velocities_out_textures', 'velocities_textures', i);
        shader_map['global_best_update_' + i] = makeGlobalBestUpdateSolver(i);
        shader_map['global_best_copy_' + i] = makeCopySolver('global_best_out_textures', 'global_best_textures', i, 2, 2);
      }

      for (let i = 0; i < this.state_textures.length; ++i) {
        for (let j = 0; j < this.env.simulation.period.length; ++j) {
          shader_map['state_textures_copy_' + i + '_' + j] = {
            vert: DefaultVertexShader,
            frag: CopyShader,
            uniforms: [
              ['original', 'tex', () => this.state_out_textures[i][j]],
            ],
            out: [this.state_textures[i][j]],
            run: this.gl_helper.runProgram,
            dims: [this.particles_width, this.particles_height],
          };

          shader_map['final_state_textures_copy_' + i + '_' + j] = {
            vert: DefaultVertexShader,
            frag: CopyShader,
            uniforms: [
              ['original', 'tex', () => this.final_state_out_textures[i][j]],
            ],
            out: [this.final_state_textures[i][j]],
            run: this.gl_helper.runProgram,
            dims: [Math.max(...this.simulation_lengths), 1],
          };
        }
      }

      if (this.env.simulation.model === 'ovvr') {
        for (let j = 0; j < this.env.simulation.period.length; ++j) {
          shader_map['state_textures_init_' + j] = makeInitOvvrSolver(j, false);
          shader_map['final_state_textures_init_' + j] = makeInitOvvrSolver(j, true);
        }
      } else if (this.env.simulation.model === 'ortp') {
        for (let j = 0; j < this.env.simulation.period.length; ++j) {
          shader_map['state_textures_init_' + j] = makeInitOrtpSolver(j, false);
          shader_map['final_state_textures_init_' + j] = makeInitOrtpSolver(j, true);
        }
      } else {
        for (let i = 0; i < this.state_textures.length; ++i) {
          for (let j = 0; j < this.env.simulation.period.length; ++j) {
            shader_map['state_textures_init_' + i + '_' + j] = makeInitStateSolver(i, j, false);
            shader_map['final_state_textures_init_' + i + '_' + j] = makeInitStateSolver(i, j, true);
          }
        }
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

      const out_array = new Float32Array(16);
      this.gl_helper.getFloatTextureArray(this.best_error_value_texture, 2, 2, out_array);
      env.particles.best_error_value = out_array[0];

      for (let t = 0; t < this.particles_textures.length; ++t) {
        this.gl_helper.getFloatTextureArray(this.global_best_textures[t], 2, 2, out_array);
        for (let i = 0; i < out_array.length; ++i) {
          env.particles.global_bests[out_array.length*t+i] = out_array[i];
        }
      }
    }

    // I am doing this like, just about as inefficiently as possible
    capture_parameters(param_arrays_map, itr) {
      const tex_width = this.tex_width;
      const tex_height = this.tex_height;
      const asize = 4 * tex_width * tex_height;

      param_arrays_map[itr] = {};

      for (let i = 0; i < this.env.particles.parameter_textures; ++i) {
        const cur_params = new Float32Array(asize);
        this.gl_helper.getFloatTextureArray(this.particles_textures[i], tex_width, tex_height, cur_params);
        param_arrays_map[itr][i] = cur_params;
      }
    }

    capture_velocities(vel_arrays_map, itr) {
      const tex_width = this.tex_width;
      const tex_height = this.tex_height;
      const asize = 4 * tex_width * tex_height;

      vel_arrays_map[itr] = {};

      for (let i = 0; i < this.env.particles.parameter_textures; ++i) {
        const cur_vels = new Float32Array(asize);
        this.gl_helper.getFloatTextureArray(this.velocities_textures[i], tex_width, tex_height, cur_vels);
        vel_arrays_map[itr][i] = cur_vels;
      }
    }

    capture_error(param_arrays_map, itr) {
      const tex_width = this.tex_width;
      const tex_height = this.tex_height;
      const asize = 4 * tex_width * tex_height;

      param_arrays_map[itr] = {};

      const cur_error = new Float32Array(asize);
      this.gl_helper.getFloatTextureArray(this.expanded_error_texture, tex_width,tex_height, cur_error);
      param_arrays_map[itr] = cur_error;
    }

    async initializeTables() {
      const program_map = this.program_map;
      const nextframe = () => new Promise(resolve => requestAnimationFrame(resolve));

      if (Object.keys(Pso.model_table_shader_map).includes(this.env.simulation.model)) {
        await nextframe();
        program_map.table();
      }
    }

    async runPrepacingIterations() {
      const program_map = this.program_map;
      const nextframe = () => new Promise(resolve => requestAnimationFrame(resolve));

      if (this.env.simulation.model === 'ovvr' || this.env.simulation.model === 'ortp') {
          for (let j = 0; j < this.env.simulation.period.length; ++j) {
            await nextframe();
            program_map['state_textures_init_' + j]();
          }
      } else {
        for (let i = 0; i < this.state_textures.length; ++i) {
          for (let j = 0; j < this.env.simulation.period.length; ++j) {
            await nextframe();
            program_map['state_textures_init_' + i + '_' + j]();
          }
        }
      }

      for (let ppb = 0; ppb < this.env.simulation.pre_beats; ++ppb) {
        for (let i = 0; i < this.env.simulation.period.length; ++i) {
          await nextframe();
          // The second argument disables clearing the output textures (only necessary for the
          // blending) and the third argument disables texture blending, which is needed for the
          // error accumulation but will not work for storing the state.
          program_map.run_simulation_prepace(i, false, true);
        }

        for (let i = 0; i < this.state_textures.length; ++i) {
          for (let j = 0; j < this.env.simulation.period.length; ++j) {
            await nextframe();
            program_map['state_textures_copy_' + i + '_' + j]();
          }
        }
      }
    }

    async runOneIteration() {
      const program_map = this.program_map;
      const nextframe = () => new Promise(resolve => requestAnimationFrame(resolve));

      await this.runPrepacingIterations();

      if (this.env.simulation.auto_normalize) {
        for (let i = 0; i < this.env.simulation.period.length; ++i) {
          await program_map.run_simulation_normalize(i, false, true);
        }
      }

      for (let i = 0; i < this.env.simulation.period.length; ++i) {
        await nextframe();
        program_map.run_simulation(i, i === 0, false);
      }

      await nextframe();
      program_map.reduce_error_1();
      await nextframe();
      program_map.reduce_error_2();

      await nextframe();
      program_map.expand_error();

      for (let i = 0; i < this.env.particles.parameter_textures; ++i) {
        await nextframe();
        program_map['global_best_update_' + i]();
      }

      for (let i = 0; i < this.env.particles.parameter_textures; ++i) {
        await nextframe();
        program_map['global_best_copy_' + i]();
      }

      await nextframe();
      program_map.best_error_value_copy();

      this.updateGlobalBest();

      for (let i = 0; i < this.env.particles.parameter_textures; ++i) {
        await nextframe();
        program_map['local_best_update_' + i]();
      }

      for (let i = 0; i < this.env.particles.parameter_textures; ++i) {
        await nextframe();
        program_map['local_bests_copy_' + i]();
      }

      await nextframe();
      program_map.local_error_copy();

      if (this.env.particles.topology === 'ring') {
        await nextframe();
        program_map.update_topological_best_ring();
      } else if (this.env.particles.topology === 'grid') {
        await nextframe();
        program_map.update_topological_best_grid();
      } else {
        await nextframe();
        program_map.update_topological_best_complete();
      }

      for (let i = 0; i < this.env.particles.parameter_textures; ++i) {
        await nextframe();
        program_map['velocity_' + i]();
        await nextframe();
        program_map.tinymt_copy();
      }

      for (let i = 0; i < this.env.particles.parameter_textures; ++i) {
        await nextframe();
        program_map['velocities_copy_' + i]();
      }

      for (let i = 0; i < this.env.particles.parameter_textures; ++i) {
        await nextframe();
        program_map['position_' + i]();
        await nextframe();
        program_map.tinymt_copy();
      }

      for (let i = 0; i < this.env.particles.parameter_textures; ++i) {
        await nextframe();
        program_map['positions_copy_' + i]();
      }
    }

    setFinalPosition(sim_length, values) {
      const num_textures = Math.ceil(values.length / 16);

      this.final_particles_textures = [];

      for (let t = 0; t < num_textures; ++t) {
        const particles_array = new Float32Array(values.slice(16*t, 16*(t+1)));
        this.final_particles_textures.push(this.gl_helper.loadFloatTexture(2, 2, particles_array));
      }
    }

    async runFinalPrepacingIterations(cl_idx) {
      const program_map = this.program_map;
      const nextframe = () => new Promise(resolve => requestAnimationFrame(resolve));

      if (this.env.simulation.model === 'ovvr' || this.env.simulation.model === 'ortp') {
          await nextframe();
          program_map['final_state_textures_init_' + cl_idx]();
      } else {
        for (let i = 0; i < this.state_textures.length; ++i) {
          await nextframe();
          program_map['final_state_textures_init_' + i + '_' + cl_idx]();
        }
      }

      for (let ppb = 0; ppb < this.env.simulation.pre_beats; ++ppb) {
        await nextframe();
        program_map.run_final_simulation_prepace(cl_idx, Math.max(...this.simulation_lengths));

        for (let i = 0; i < this.state_textures.length; ++i) {
          await nextframe();
          program_map['final_state_textures_copy_' + i + '_' + cl_idx]();
        }
      }
    }

    async runFinalSimulationSolver(cl_idx, values) {
      const simsize = this.simulation_lengths[cl_idx];
      const texsize = simsize*4;

      const parameter_values = values || this.env.particles.global_bests;

      // In the case of the global_bests array, it should already have length divisible by 16, since
      // it is a Float32Array and cannot be pushed to.
      while (parameter_values.length % 16 !== 0) {
        parameter_values.push(0);
      }

      this.setFinalPosition(simsize, parameter_values);

      await this.runFinalPrepacingIterations(cl_idx);

      if (this.env.simulation.auto_normalize) {
        await this.program_map.run_final_simulation_normalize(cl_idx, Math.max(...this.simulation_lengths));
      }

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
