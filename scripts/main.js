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

  const graph_canvas = document.getElementById('graph_canvas');
  const error_canvas = document.getElementById('error_canvas');

  let pso;
  const graph = new Graph(graph_canvas);
  const error_graph = new Graph(error_canvas);
  const pso_interface = new PsoInterface();

  const dump_convergence = false;
  const save_error = false;

  pso_interface.displayBounds(Pso.getEnv());
  pso_interface.displayModelParameters();
  pso_interface.displayStimulusParameters();

  document.getElementById('reset_bounds').onclick = () => {
    pso_interface.displayBounds(Pso.getEnv());
  };


  pso_interface.add_button.onclick = () => pso_interface.addInput();
  pso_interface.remove_button.onclick = () => pso_interface.removeInput();
  pso_interface.fit_all_button.onclick = () => pso_interface.setFitCheckboxes(true);
  pso_interface.fit_none_button.onclick = () => pso_interface.setFitCheckboxes(false);
  pso_interface.normalize.onclick = () => pso_interface.updateNormalizationDisplay();
  pso_interface.plot_from_vals_button.onclick = () => {
    if (!pso) {
      alert('A fit must be created before modifying the parameters');
      return;
    }

    displayGraph(pso_interface.plotting_idx, pso_interface.get_current_values());
  };

  pso_interface.model_select.addEventListener('change', () => pso_interface.displayModelParameters());
  pso_interface.stim_biphasic_checkbox.addEventListener('change', () => pso_interface.displayStimulusParameters());

  const save_output = (data, filename) => {
    const file = new Blob(data, { type: 'text/plain' });
    const e = document.createElement('a');
    const url = URL.createObjectURL(file);
    e.href = url;
    e.download = filename;
    document.body.appendChild(e);
    e.click();
    setTimeout(() => {
      document.body.removeChild(e);
      window.URL.revokeObjectURL(url);
    }, 0);
  };

  pso_interface.save_params_button.onclick = () => {
    const model = pso.env.simulation.model;
    const global_bests = pso.env.particles.global_bests;
    const params = PsoInterface.param_lists[model];

    const params_obj = {};
    for (let i = 0; i < params.length; ++i) {
      params_obj[params[i]] = global_bests[i];
    }

    save_output([JSON.stringify(params_obj, null, 2)], `pso_${model}_params_${Date.now()}.json`);
  };

  pso_interface.save_run_button.onclick = () => {
    const details_obj = {
      'simulation': (({
        model,
        dt,
        period,
        num_beats,
        pre_beats,
        datatypes,
        weights,
        sample_interval,
        normalize,
        normalization_max,
        normalization_min,
        normalized_align_threshold,
        normalized_ca_align_threshold,
      }) => ({
        model,
        dt,
        period,
        num_beats,
        pre_beats,
        datatypes,
        weights,
        sample_interval,
        normalize,
        normalization_max,
        normalization_min,
        normalized_align_threshold,
        normalized_ca_align_threshold,
      }))(pso.env.simulation),
      'stimulus': structuredClone(pso.env.stimulus),
      'particles': structuredClone(pso.env.particles),
    };

    const params = PsoInterface.param_lists[pso.env.simulation.model];

    const params_obj = {};
    const lower_bounds_obj = {};
    const upper_bounds_obj = {};

    for (let i = 0; i < params.length; ++i) {
      const param = params[i];
      params_obj[param] = pso.env.particles.global_bests[i];
      lower_bounds_obj[param] = pso.env.particles.lower_bounds[i];
      upper_bounds_obj[param] = pso.env.particles.upper_bounds[i];
    }

    details_obj.particles.global_bests = params_obj;
    details_obj.particles.lower_bounds = lower_bounds_obj;
    details_obj.particles.upper_bounds = upper_bounds_obj;

    save_output([JSON.stringify(details_obj, null, 2)], `pso_run_${Date.now()}.json`);
  };

  pso_interface.data_section.onclick = async (e) => {
    if (e.target.getAttribute('class') === 'plot-data-button') {
      const idx = Array.from(pso_interface.data_section.children).indexOf(e.target.parentElement);

      if (idx !== -1) {
        pso_interface.update_plot_idx(idx);
        if (pso && idx < pso.env.simulation.period.length) {
          displayGraph(idx);
        } else {
          displayDataGraph(idx).catch(err => alert(err));
        }
      }
    }
  };

  // Start with a single data file input
  pso_interface.addInput();

  const run_pso = async () => {
    const hyperparams = pso_interface.getHyperparams();
    pso = new Pso(hyperparams.particle_count);
    const input_data = await pso_interface.getAllInputData();

    pso.setupEnv(
      pso_interface.model_select.value,
      pso_interface.getBounds(),
      pso_interface.getStimulusParameters(),
      pso_interface.data_pre_beats.value,
      pso_interface.data_num_beats.value,
      pso_interface.data_sample_interval.value,
      pso_interface.normalize.checked,
      Number(pso_interface.normalization_max.value),
      Number(pso_interface.normalization_min.value),
      hyperparams,
    );

    // TODO: fail if calcium data is provided for non-calcium model
    pso.readData(input_data);

    pso.initializeTextures();
    pso.setupAllSolvers();

    runPsoIterations(hyperparams.iteration_count);
  };

  async function runPsoIterations(iteration_count) {
    const start_time = Date.now();
    const best_error_list = [];

    const iteration_params = {};
    const iteration_error = {};
    const iteration_velocities = {};

    for (let iter = 0; iter < iteration_count; ++iter) {
      pso_interface.updateStatusDisplay(iter, iteration_count);
      await pso.runOneIteration();

      if (dump_convergence) {
        pso.capture_parameters(iteration_params,iter);
        pso.capture_error(iteration_error,iter);
        pso.capture_velocities(iteration_velocities,iter);
      }

      best_error_list.push(pso.env.particles.best_error_value);
    }

    pso_interface.updateStatusDisplay(iteration_count, iteration_count);

    finalizePso(start_time, best_error_list);

    if (save_error) {
      const filename = `pso_error_${pso.env.simulation.model}_${pso.particles_width*pso.particles_height}_${iteration_count}_${Date.now()}.txt`;
      save_output([best_error_list.join('\n')], filename);
    }

    if (dump_convergence) {
      save_output([JSON.stringify(iteration_params)], `convergence_data.json`);
      save_output([JSON.stringify(iteration_error)], `error_data.json`);
      save_output([JSON.stringify(iteration_velocities)], `velocity_data.json`);
    }
  }

  async function finalizePso(start_time, best_error_list) {
    const bestArr = pso.env.particles.global_bests;
    pso_interface.displayResults(bestArr);
    pso_interface.displayError(pso.env.particles.best_error_value);

    const errmin = Math.min(...best_error_list);
    const errmax = Math.max(...best_error_list);
    error_graph.clearGraph();
    error_graph.runGraph(best_error_list, [0, 0, 0], best_error_list.length, [errmin, errmax]);
    pso_interface.setErrorAxes(1, best_error_list.length, errmin, errmax);

    displayGraph(pso_interface.plotting_idx);

    console.log(`Execution time (ms): ${Date.now()-start_time}`);
  }

  async function displayDataGraph(cl_idx) {
    const input_data = await pso_interface.getAllInputData();
    if (input_data[cl_idx].datatype == 'apd') {
      return;
    }
    const raw_text = input_data[cl_idx].data;
    const actual_data = raw_text.split('\n').filter(x => !(x.trim() === ""));

    const scale = [Math.min(...actual_data), Math.max(...actual_data)];

    graph.clearGraph();
    graph.runGraph(actual_data, [0, 0, 0], actual_data.length, scale);
  }

  const apd_start_indices = (data, apd_thresh) => {
    let in_ap = false;
    const idxs = [];

    for (let i = 0; i < data.length; ++i) {
      if (!in_ap && data[i] > apd_thresh) {
        in_ap = true;
        idxs.push(i);
      } else if (in_ap && data[i] < apd_thresh) {
        in_ap = false;
      }
    }

    return idxs;
  };

  function displayGraph(cl_idx, current_values) {
    const simulation_data = current_values ?
      pso.runFinalSimulationSolver(cl_idx, current_values) :
      pso.runFinalSimulationSolver(cl_idx);

    const sim_length = (pso.env.simulation.period[cl_idx] * pso.env.simulation.num_beats) / pso.env.simulation.sample_interval;
    const plotting_sim_data = simulation_data.slice(-sim_length);

    let actual_data = [];
    let apd_starts = [];
    let apd_ends = [];
    let align_index = 0;

    if (pso.env.simulation.datatypes[cl_idx] === 'apd') {
      apd_starts = apd_start_indices(plotting_sim_data, pso.env.simulation.apd_threshs[cl_idx]);
      const apds = pso.env.simulation.trimmed_data[cl_idx];
      apd_ends = apd_starts.map((x, idx) => x + (apds[idx] || 0) / pso.env.simulation.sample_interval);
    } else {
      actual_data = pso.env.simulation.trimmed_data[cl_idx];
      const align_thresh = pso.env.simulation.align_thresh[cl_idx];
      align_index = plotting_sim_data.findIndex(x => x > align_thresh);
    }

    const scale = [
      Math.min(...actual_data, ...plotting_sim_data),
      Math.max(...actual_data, ...plotting_sim_data),
    ];

    const interval = Number(pso_interface.data_sample_interval.value);

    graph.clearGraph();
    if (pso.env.simulation.datatypes[cl_idx] === 'apd') {
      for (let i = 0; i < apd_starts.length; ++i) {
        graph.runApdGraph(apd_starts[i], apd_ends[i], pso.env.simulation.apd_threshs[cl_idx], [0, 0, 0], sim_length, scale, 0);
      }
    } else {
      graph.runGraph(actual_data, [0, 0, 0], sim_length, scale, align_index);
    }
    graph.runGraph(plotting_sim_data, [1, 0, 0], sim_length, scale, 0);

    pso_interface.setAxes(0, sim_length * interval, scale[0], scale[1]);
  }

  document.querySelector('button#PSO_button').onclick = () => run_pso();
});
