/* global require */
require([
  'scripts/graph',
  'scripts/interface',
  'scripts/pso',
], function(
  Graph,
  PsoInterface,
  Pso,
) {
  'use strict';

  const particles_width = 32;
  const particles_height = 32;

  let pso;
  const graph = new Graph();
  const pso_interface = new PsoInterface();

  pso_interface.displayBounds(Pso.getEnv());

  document.getElementById('reset_bounds').onclick = () => {
    pso_interface.displayBounds(Pso.getEnv());
  };

  pso_interface.add_button.onclick = () => pso_interface.addInput();
  pso_interface.remove_button.onclick = () => pso_interface.removeInput();
  pso_interface.fit_all_button.onclick = () => pso_interface.setFitCheckboxes(true);
  pso_interface.fit_none_button.onclick = () => pso_interface.setFitCheckboxes(false);

  pso_interface.data_section.onclick = (e) => {
    if (e.target.getAttribute('class') === 'plot-data-button') {
      const idx = Array.from(pso_interface.data_section.children).indexOf(e.target.parentElement);
      if (idx !== -1 && pso) {
        displayGraph(idx);
      }
    }
  };

  // Start with a single data file input
  pso_interface.addInput();

  const outer_display = () => {
    pso_interface.display_all_params(pso_interface);
  };

  const run_pso = async () => {
    pso = new Pso(particles_width, particles_height);
    const input_data = await pso_interface.getAllInputData();

    const raw_input_data = [];
    const input_cls = [];
    for (const [raw, cl] of input_data) {
      raw_input_data.push(raw);
      input_cls.push(cl);
    }

    const start_time = Date.now();

    pso.setupEnv(
      pso_interface.getBounds(),
      input_cls,
      pso_interface.data_pre_beats.value,
      pso_interface.data_num_beats.value,
      pso_interface.data_sample_interval.value,
    );

    pso.readData(raw_input_data, pso_interface.normalization.value);

    pso.initializeTextures();
    pso.setupAllSolvers();

    for (let i = 0; i < 32; ++i) {
      console.log(pso.env.particles.best_error_value);
      pso.runOneIteration();
    }

    const bestArr = pso.getGlobalBests();
    pso_interface.displayResults(bestArr);
    pso_interface.displayError(pso.env.particles.best_error_value);

    displayGraph(0);

    console.log("Execution time (ms):");
    console.log(Date.now() - start_time);
  };

  function displayGraph(cl_idx) {
    const simulation_data = pso.runFinalSimulationSolver(cl_idx);
    const actual_data = pso.env.simulation.trimmed_data[cl_idx];

    const align_index = simulation_data.findIndex(number => number > 0.15);
    const plotting_sim_data = simulation_data.slice(align_index);

    const scale = [
      Math.min(...actual_data, ...plotting_sim_data),
      Math.max(...actual_data, ...plotting_sim_data),
    ];

    const num_points = Math.max(actual_data.length, plotting_sim_data.length);
    const interval = Number(pso_interface.data_sample_interval.value);

    graph.clearGraph();
    graph.runGraph(actual_data, [0, 0, 0], num_points, scale);
    graph.runGraph(plotting_sim_data, [1, 0, 0], num_points, scale);

    pso_interface.setAxes(0, num_points * interval, scale[0], scale[1]);
  }

  document.querySelector('button#PSO_button').onclick = run_pso;
  document.querySelector('button#disp_params_button').onclick = outer_display;
});
