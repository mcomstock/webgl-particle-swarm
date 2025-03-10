define('scripts/interface', [
], function(
) {
  'use strict';

  return class PsoInterface {
    static param_lists = {
      fk: ['tr', 'tsi', 'twp', 'td', 'tvp', 'tv1m', 'tv2m', 'twm', 'to', 'xk', 'ucsi', 'uc', 'uv'],
      ms: ['tin', 'tout', 'tclose', 'topen', 'vgate'],
      mms: ['mms_tin', 'mms_tout', 'mms_tclose', 'mms_topen', 'mms_vgate'],
      fhn: ['alpha', 'beta', 'eps', 'mu', 'gamma', 'theta', 'delta'],
      b4v: ['b4v_thv',  'b4v_tv1m',  'b4v_tv2m',  'b4v_tvp',  'b4v_uwm',  'b4v_tso1',  'b4v_kso',  'b4v_ts1',  'b4v_ts2',  'b4v_ks',  'b4v_tw1m',  'b4v_tw2m',  'b4v_tw1p',  'b4v_tfi',  'b4v_to1',  'b4v_to2',  'b4v_tso2',  'b4v_uso',  'b4v_us' ,  'b4v_tsi1', 'b4v_thw', 'b4v_thvm', 'b4v_tho', 'b4v_kwm', 'b4v_twinf', 'b4v_winfstar', 'b4v_uu'],
      bb: ['bb_tv1p', 'bb_tv1m', 'bb_tv2m', 'bb_tw1p', 'bb_tw2p', 'bb_tw1m', 'bb_tw2m', 'bb_ts1', 'bb_ts2', 'bb_tfi', 'bb_to1', 'bb_to2', 'bb_tso1', 'bb_tso2', 'bb_tsi1', 'bb_tsi2', 'bb_twinf', 'bb_thv', 'bb_thvm', 'bb_thvinf', 'bb_thw', 'bb_thwinf', 'bb_thso', 'bb_thsi', 'bb_tho', 'bb_ths', 'bb_kwp', 'bb_kwm', 'bb_ks', 'bb_kso', 'bb_ksi', 'bb_uwm', 'bb_us', 'bb_uo', 'bb_uu', 'bb_uso', 'bb_sc', 'bb_wcp', 'bb_winfstar'],
      tnnp2006: ['tnnp2006_gna', 'tnnp2006_gk1', 'tnnp2006_gto', 'tnnp2006_gkr', 'tnnp2006_gks', 'tnnp2006_gcal', 'tnnp2006_gpk', 'tnnp2006_gpca', 'tnnp2006_gbna', 'tnnp2006_gbca', 'tnnp2006_pnak', 'tnnp2006_knaca'],
      ovvr: ['ovvr_gnafast', 'ovvr_gnalate', 'ovvr_gto', 'ovvr_pca', 'ovvr_pcana', 'ovvr_pcak', 'ovvr_pcacamk', 'ovvr_pcanacamk', 'ovvr_pcakcamk', 'ovvr_gkr', 'ovvr_gks', 'ovvr_gk1', 'ovvr_gnaca', 'ovvr_gnak', 'ovvr_pnab', 'ovvr_pcab', 'ovvr_gkb', 'ovvr_gpca'],
      ortp: ['ortp_gna', 'ortp_gto', 'ortp_pca', 'ortp_pcana', 'ortp_pcak', 'ortp_pcacamk', 'ortp_pcanacamk', 'ortp_pcakcamk', 'ortp_gkr', 'ortp_gks', 'ortp_gk1', 'ortp_gnaca', 'ortp_gnak', 'ortp_pnab', 'ortp_pcab', 'ortp_gkb', 'ortp_gpca'],
    };

    static default_normalization = {
      fk: {
        normalize: true,
        max: 1,
        min: 0,
      },
      ms: {
        normalize: true,
        max: 1,
        min: 0,
      },
      mms: {
        normalize: true,
        max: 1,
        min: 0,
      },
      fhn: {
        normalize: true,
        max: 1,
        min: 0,
      },
      b4v: {
        normalize: true,
        max: 1.2,
        min: 0,
      },
      bb: {
        normalize: true,
        max: 1.2,
        min: 0,
      },
      tnnp2006: {
        normalize: false,
        max: 36.05,
        min: -86.3,
      },
      ovvr: {
        normalize: false,
        max: 39.92,
        min: -87.84,
      },
      ortp: {
        normalize: false,
        max: 43.29,
        min: -87.85,
      },
    };

    static data_type_names = {
      'voltage': 'Voltage',
      'apd': 'APD',
      'calcium': 'Calcium',
    };

    constructor() {
      for (const [model, param_list] of Object.entries(PsoInterface.param_lists)) {
        const model_elements = {};

        for (const param of param_list) {
          model_elements[param + '_val'] = document.getElementById(param + '_val');
          model_elements[param + '_min'] = document.getElementById(param + '_min');
          model_elements[param + '_max'] = document.getElementById(param + '_max');
          model_elements[param + '_fit'] = document.getElementById(param + '_fit');

          model_elements[param + '_prev_min'] = -1;
          model_elements[param + '_prev_max'] = -1;

          model_elements[param + '_fit'].addEventListener('change', () => {
            if (model_elements[param + '_fit'].checked) {
              model_elements[param + '_val'].removeAttribute('readonly');
              model_elements[param + '_min'].removeAttribute('readonly');
              model_elements[param + '_max'].removeAttribute('readonly');

              model_elements[param + '_min'].value = model_elements[param + '_prev_min'];
              model_elements[param + '_max'].value = model_elements[param + '_prev_max'];
            } else {
              model_elements[param+'_prev_min'] = model_elements[param + '_min'].value;
              model_elements[param + '_min'].value = model_elements[param + '_val'].value;

              model_elements[param+'_prev_max'] = model_elements[param + '_max'].value;
              model_elements[param + '_max'].value = model_elements[param + '_val'].value;

              model_elements[param + '_val'].setAttribute('readonly', true);
              model_elements[param + '_min'].setAttribute('readonly', true);
              model_elements[param + '_max'].setAttribute('readonly', true);
            }
          });
        }

        this[model] = model_elements;
      }

      this.normalize = document.getElementById('normalize');
      this.normalization_max = document.getElementById('normalization-max');
      this.normalization_min = document.getElementById('normalization-min');
      this.data_num_beats = document.getElementById('data_num_beats');
      this.data_pre_beats = document.getElementById('data_pre_beats');
      this.data_sample_interval = document.getElementById('data_sample_interval');
      this.fit_error = document.getElementById('fit_error');
      this.data_section = document.getElementById('data-section');
      this.add_button = document.getElementById('add-data');
      this.remove_button = document.getElementById('remove-data');
      this.fit_all_button = document.getElementById('fit-all-button');
      this.fit_none_button = document.getElementById('fit-none-button');
      this.plot_from_vals_button = document.getElementById("plot-from-vals-button");
      this.model_select = document.getElementById('model-select');
      this.status_display = document.getElementById('status_display');
      this.save_params_button = document.getElementById('save_params_button');
      this.save_run_button = document.getElementById('save_run_button');
      this.stim_biphasic_checkbox = document.getElementById('stim-biphasic');
      this.stim_square = {
        section: document.getElementById('stim-square-section'),
        stim_dur: document.getElementById('stim-dur-square'),
        stim_mag: document.getElementById('stim-mag-square'),
      };
      this.stim_biphasic = {
        section: document.getElementById('stim-biphasic-section'),
        stim_dur: document.getElementById('stim-dur-biphasic'),
        stim_mag: document.getElementById('stim-mag-biphasic'),
        stim_offset_1: document.getElementById('stim-offset-1-biphasic'),
        stim_offset_2: document.getElementById('stim-offset-2-biphasic'),
        stim_t_scale: document.getElementById('stim-t-scale-biphasic'),
      };

      this.default_button_bg = null;
      this.plotting_idx = -1;
      this.active_plot_bg = "rgb(174,216,230)";


      this.xmin = document.getElementById('xmin');
      this.x1 = document.getElementById('x1');
      this.x2 = document.getElementById('x2');
      this.x3 = document.getElementById('x3');
      this.xmax = document.getElementById('xmax');
      this.ymin = document.getElementById('ymin');
      this.y1 = document.getElementById('y1');
      this.y2 = document.getElementById('y2');
      this.y3 = document.getElementById('y3');
      this.ymax = document.getElementById('ymax');

      this.error_xmin = document.getElementById('error_xmin');
      this.error_xmax = document.getElementById('error_xmax');
      this.error_ymin = document.getElementById('error_ymin');
      this.error_ymax = document.getElementById('error_ymax');
    }

    setDefaultNormalization() {
      const settings = PsoInterface.default_normalization[this.model_select.value];
      this.normalize.checked = settings.normalize;
      this.normalization_max.value = 'max' in settings ? settings.max : '';
      this.normalization_min.value = 'min' in settings ? settings.min : '';
      this.updateNormalizationDisplay();
    }

    updateNormalizationDisplay() {
      if (this.normalize.checked) {
        this.normalization_max.disabled = false;
        this.normalization_min.disabled = false;
      } else {
        this.normalization_max.disabled = true;
        this.normalization_min.disabled = true;
      }
    }

    updateStatusDisplay(current, total) {
      let str = `Iteration: ${current}/${total}`;
      if (current === total) {
        str += ' Done!';
      }
      this.status_display.innerHTML = str;
    }

    displayStimulusParameters() {
      if (this.stim_biphasic_checkbox.checked) {
        this.stim_square.section.classList.add('data-input-hidden');
        this.stim_biphasic.section.classList.remove('data-input-hidden');
      } else {
        this.stim_square.section.classList.remove('data-input-hidden');
        this.stim_biphasic.section.classList.add('data-input-hidden');
      }
    }

    displayModelParameters() {
      const selected_id = this.model_select.value + '-parameters';
      const model_param_divs = document.querySelectorAll('div#parameter-section div.param-div,div.param-div-hidden');

      for (const model_param_div of model_param_divs) {
        if (model_param_div.id === selected_id) {
          model_param_div.setAttribute('class', 'param-div');
        } else {
          model_param_div.setAttribute('class', 'param-div-hidden');
        }
      }

      this.setDefaultNormalization();
    }

    setFitCheckboxes(value) {
      const fit_checkboxes = document.querySelectorAll('input.fit-checkbox');

      for (const cb of fit_checkboxes) {
        cb.checked = !value;
        cb.click();
      }
    }

    async getDataFromInput(element) {
      const cl = Number(element.querySelector('.data-cl-input').value);
      const data_type = element.querySelector('.data-type-select').value;
      const weight = Number(element.querySelector('.data-weight-input').value);

      if (cl === 0) {
        throw new Error('No CL entered');
      }

      if (data_type === 'apd') {
        const apds = new Float32Array(element.querySelector('.data-apd-input').value.split(','));
        const thresh = Number(element.querySelector('.data-apd-thresh-input').value);

        if (apds.length === 0) {
          throw new Error('No APDs entered');
        }

        return {
          datatype: data_type,
          data: apds,
          cl: cl,
          apd_thresh: thresh,
          weight: weight,
        };
      }

      const file = element.querySelector('.data-file-input').files[0];

      if (!file) {
        throw new Error('No file found');
      }

      const reader = new FileReader();

      const text = await new Promise((resolve) => {
        reader.onload = () => {
          resolve(reader.result);
        };

        try {
          reader.readAsText(file);
        } catch (e) {
          throw new Error('Could not read file');
        }
      });

      return {
        datatype: data_type,
        data: text,
        cl: cl,
        weight: weight,
      };
    }

    async getAllInputData() {
      const p_result = Array.from(this.data_section.children).map(async (element) => await this.getDataFromInput(element));
      const result = await Promise.all(p_result);

      return result;
    }

    createInputElement() {
      const elem = document.createElement('div');
      elem.classList.add('data-input');

      const file_in = document.createElement('input');
      file_in.setAttribute('type', 'file');
      file_in.classList.add('data-file-input');

      const apd_label = document.createElement('label');
      apd_label.innerHTML = 'APDs';

      const apd_input = document.createElement('input');
      apd_input.setAttribute('type', 'text');
      apd_input.classList.add('data-apd-input');
      apd_label.appendChild(apd_input);

      const apd_thresh_label = document.createElement('label');
      apd_thresh_label.innerHTML = 'APD threshold';

      const apd_thresh_input = document.createElement('input');
      apd_thresh_input.setAttribute('type', 'text');
      apd_thresh_input.classList.add('data-apd-thresh-input');
      apd_thresh_label.appendChild(apd_thresh_input);
      apd_thresh_input.value = String(0.15);

      const cl_label = document.createElement('label');
      cl_label.innerHTML = 'Cycle length (ms)';

      const cl_in = document.createElement('input');
      cl_in.setAttribute('type', 'text');
      cl_in.classList.add('data-cl-input');
      cl_label.appendChild(cl_in);

      const weight_label = document.createElement('label');
      weight_label.innerHTML = 'Fitting weight';

      const weight_in = document.createElement('input');
      weight_in.setAttribute('type', 'text');
      weight_in.classList.add('data-weight-input');
      weight_label.appendChild(weight_in);
      weight_in.value = "1";

      const data_select_label = document.createElement('label');
      data_select_label.innerHTML = 'Data type';

      const data_select = document.createElement('select');
      data_select.setAttribute('name', 'data-type');
      data_select.classList.add('data-type-select');

      for (const [data_type, data_type_name] of Object.entries(PsoInterface.data_type_names)) {
        const data_select_option = document.createElement('option');
        data_select_option.setAttribute('value', data_type);
        data_select_option.innerHTML = data_type_name;
        if (data_type === 'voltage') {
          data_select_option.setAttribute('selected', '');
        }
        data_select.appendChild(data_select_option);
      }

      data_select_label.appendChild(data_select);

      const plot_button = document.createElement('button');
      plot_button.setAttribute('type', 'button');
      plot_button.classList.add('plot-data-button');
      plot_button.innerHTML = 'Plot';

      if (this.default_button_bg === null) {
        this.default_button_bg = plot_button.style.backgroundColor;
      }

      if (this.plotting_idx < 0) {
        this.plotting_idx = 0;
        plot_button.style.backgroundColor = this.active_plot_bg;
      }

      elem.appendChild(file_in);
      elem.appendChild(apd_label);
      elem.appendChild(apd_thresh_label);
      elem.appendChild(cl_label);
      elem.appendChild(weight_label);
      elem.appendChild(data_select_label);
      elem.appendChild(plot_button);

      const apd_elements = [apd_label, apd_thresh_label];
      const file_elements = [file_in];

      const set_hidden = () => {
        if (data_select.value === 'apd') {
          for (const e of apd_elements) {
            e.classList.remove('data-input-hidden');
          }
          for (const e of file_elements) {
            e.classList.add('data-input-hidden');
          }
        } else {
          for (const e of file_elements) {
            e.classList.remove('data-input-hidden');
          }
          for (const e of apd_elements) {
            e.classList.add('data-input-hidden');
          }
        }
      };

      data_select.addEventListener('change', set_hidden);

      set_hidden();

      return elem;
    }

    update_plot_idx(new_idx) {
      if (new_idx !== this.plotting_idx) {
        const tst = Array.from(this.data_section.children);
        tst[this.plotting_idx].querySelector('.plot-data-button').style.backgroundColor = this.default_button_bg;
        tst[new_idx].querySelector('.plot-data-button').style.backgroundColor = this.active_plot_bg;
        this.plotting_idx = new_idx;
      }
    }

    get_plot_idx() {
      return this.plotting_idx;
    }


    addInput() {
      this.data_section.appendChild(this.createInputElement());
    }

    removeInput() {
      const children = this.data_section.children;

      if (children.length > 1) {
        this.data_section.removeChild(children[children.length-1]);
      }
    }

    displayBounds(env) {
      for (const [model, param_list] of Object.entries(PsoInterface.param_lists)) {
        const [lb, ub] = env[model + '_bounds'];

        param_list.forEach((param, idx) => {
          if (this[model][param + '_fit'].checked) {
            this[model][param + '_min'].value = lb[idx];
            this[model][param + '_max'].value = ub[idx];
          }

        });
      }
    }

    displayResults(bestArr) {
      const model = this.model_select.value;
      PsoInterface.param_lists[model].forEach((param, idx) => {
        // this[model][param + '_val'].value = bestArr[idx];
        this[model][param + '_val'].value = bestArr[idx].toFixed(3);

      });
    }

    get_current_values() {
      const val_arr = [];
      const model = this.model_select.value;
      PsoInterface.param_lists[model].forEach((param) => {
        val_arr.push(Number(this[model][param + "_val"].value));
      });
      return val_arr;
    }


    display_all_params() {
      let builder = "";
      PsoInterface.param_lists[this.model_select.value].forEach(param => {
        builder = builder.concat(param + ":\t" + this[param + '_val'].value +"\n");
      });
      prompt("(ctrl+c, Enter) to copy:",builder);
    }

    displayError(error) {
      this.fit_error.innerHTML = (Number.parseFloat(error).toPrecision(4)).toString();
    }

    getBounds() {
      const lb = [], ub = [];
      const model = this.model_select.value;
      const param_list = PsoInterface.param_lists[model];

      for (const param of param_list) {
        lb.push(Number(this[model][param + '_min'].value));
        ub.push(Number(this[model][param + '_max'].value));
      }

      return [lb, ub];
    }

    getStimulusParameters() {
      return this.stim_biphasic_checkbox.checked ? {
        stim_dur: Number(this.stim_biphasic.stim_dur.value),
        stim_mag: Number(this.stim_biphasic.stim_mag.value),
        stim_biphasic: true,
        stim_offset_1: Number(this.stim_biphasic.stim_offset_1.value),
        stim_offset_2: Number(this.stim_biphasic.stim_offset_2.value),
        stim_t_scale: Number(this.stim_biphasic.stim_t_scale.value),
      } : {
        stim_dur: Number(this.stim_square.stim_dur.value),
        stim_mag: Number(this.stim_square.stim_mag.value),
        stim_biphasic: false,
        stim_offset_1: 0.0,
        stim_offset_2: 0.0,
        stim_t_scale: 0.0,
      };
    }

    setAxes(xmin, xmax, ymin, ymax) {
      const x_quart = (xmax - xmin) / 4;
      const y_quart = (ymax - ymin) / 4;

      this.xmin.innerHTML = this.truncateString(xmin);
      this.x1.innerHTML = this.truncateString(xmin + x_quart);
      this.x2.innerHTML = this.truncateString(xmin + 2*x_quart);
      this.x3.innerHTML = this.truncateString(xmin + 3*x_quart);
      this.xmax.innerHTML = this.truncateString(xmax);

      this.ymin.innerHTML = this.truncateString(ymin);
      this.y1.innerHTML = this.truncateString(ymin + y_quart);
      this.y2.innerHTML = this.truncateString(ymin + 2*y_quart);
      this.y3.innerHTML = this.truncateString(ymin + 3*y_quart);
      this.ymax.innerHTML = this.truncateString(ymax);
    }

    setErrorAxes(xmin, xmax, ymin, ymax) {
      this.error_xmin.innerHTML = String(xmin);
      this.error_xmax.innerHTML = String(xmax);
      this.error_ymin.innerHTML = Number.parseFloat(ymin).toPrecision(4).toString();
      this.error_ymax.innerHTML = Number.parseFloat(ymax).toPrecision(4).toString();
    }

    truncateString(num) {
      return Number.parseFloat(num).toFixed(2).toString();
    }

    getHyperparams() {
      const hyperparams = {};
      hyperparams.particle_count = Number(document.querySelector('input[name="particle_count"]:checked').value);
      hyperparams.iteration_count = Number(document.getElementById('hyperparam_iteration_count').value);
      hyperparams.phi1 = Number(document.getElementById('hyperparam_phi1').value);
      hyperparams.phi2 = Number(document.getElementById('hyperparam_phi2').value);
      hyperparams.chi = Number(document.getElementById('hyperparam_chi').value);

      return hyperparams;
    }
  };
});
