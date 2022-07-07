/* global define */
define('scripts/interface', [
], function(
) {
  'use strict';

  /* global PsoInterface */
  return class PsoInterface {
    static param_lists = {
      fk: ['tr', 'tsi', 'twp', 'td', 'tvp', 'tv1m', 'tv2m', 'twm', 'to', 'xk', 'ucsi', 'uc', 'uv'],
      ms: ['gna', 'gk', 'tclose', 'topen', 'vgate'],
    }

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

      this.normalization = document.getElementById('normalization');
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
    }

    setFitCheckboxes(value) {
      const fit_checkboxes = document.querySelectorAll('input.fit-checkbox');

      for (const cb of fit_checkboxes) {
        cb.checked = !value;
        cb.click();
      }
    }

    async getDataFromInput(element) {
      const file = element.querySelector('input[type=file]').files[0];
      const cl = Number(element.querySelector('input[type=text]').value);

      const reader = new FileReader();

      const text = await new Promise((resolve) => {
        reader.onload = () => {
          resolve(reader.result);
        };

        reader.readAsText(file);
      });

      return [text, cl];
    }

    async getAllInputData() {
      const p_result = Array.from(this.data_section.children).map(async (element) => await this.getDataFromInput(element));
      const result = await Promise.all(p_result);
      return result;
    }

    createInputElement() {
      const elem = document.createElement('div');
      elem.setAttribute('class', 'data-input');

      const file_in = document.createElement('input');
      file_in.setAttribute('type', 'file');

      const cl_label = document.createElement('span');
      cl_label.innerHTML = 'Cycle length (ms):';

      const cl_in = document.createElement('input');
      cl_in.setAttribute('type', 'text');

      const plot_button = document.createElement('button');
      plot_button.setAttribute('type', 'button');
      plot_button.setAttribute('class', 'plot-data-button');
      plot_button.innerHTML = 'Plot';

      if(this.default_button_bg == null)
      {
        this.default_button_bg = plot_button.style.backgroundColor;
      }

      if(this.plotting_idx < 0)
      {
        this.plotting_idx = 0;
        plot_button.style.backgroundColor = this.active_plot_bg;
      }

      elem.appendChild(file_in);
      elem.appendChild(cl_label);
      elem.appendChild(cl_in);
      elem.appendChild(plot_button);

      return elem;
    }

    update_plot_idx(new_idx)
    {
      if(new_idx != this.plotting_idx)
      {
        var tst = Array.from(this.data_section.children);        
        tst[this.plotting_idx].querySelector('.plot-data-button').style.backgroundColor = this.default_button_bg;
        tst[new_idx].querySelector('.plot-data-button').style.backgroundColor = this.active_plot_bg;
        this.plotting_idx = new_idx;
      }
    }

    get_plot_idx()
    {
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
        const bounds = env[model + '_bounds'].flat(1);

        param_list.forEach((param, idx) => {
          if (this[model][param + '_fit'].checked) {
            const [min, max] = bounds[idx];
            this[model][param + '_min'].value = min;
            this[model][param + '_max'].value = max;
          }
        });
      }

      this.normalization.value = 1;
    }

    displayResults(bestArr) {
      const model = this.model_select.value;
      PsoInterface.param_lists[model].forEach((param, idx) => {
        // this[model][param + '_val'].value = bestArr[idx];
        this[model][param + '_val'].value = bestArr[idx].toFixed(3);

      });
    }

    display_all_params_test()
    {
      var builder = "";
      const model = this.model_select.value;
      PsoInterface.param_lists[model].forEach((param,idx) => 
          {
            builder = builder.concat(param + ":\t"+ this[model][param + "_val"].value + "\n");
          }
        );
      prompt("(ctrl+c, Enter) to copy:",builder);
    }

    get_current_values()
    {
      var val_arr = [];
      const model = this.model_select.value;
      PsoInterface.param_lists[model].forEach((param,idx) => 
          {
            val_arr.push(Number(this[model][param + "_val"].value));
          }
        );
      return val_arr;
    }


    display_all_params(real_interface)
    {
      var builder = "";
      PsoInterface.param_lists[this.model_select.value].forEach(param => {
        builder = builder.concat(param + ":\t" + this[param + '_val'].value +"\n");
      });
      prompt("(ctrl+c, Enter) to copy:",builder);
    }

    displayError(error) {
      // this.fit_error.innerHTML = this.truncateString(error);
      this.fit_error.innerHTML = (Number.parseFloat(error).toPrecision(4)).toString();
    }

    getBounds() {
      const bounds = [];
      const model = this.model_select.value;
      const param_list = PsoInterface.param_lists[model];
      const num_colors = Math.ceil(param_list.length/4);

      let c, i;
      for (i = 0; i < param_list.length; ++i) {
        const param = param_list[i];

        if (i % 4 === 0) {
          bounds.push([]);
          c = i/4;
        }

        bounds[c].push([Number(this[model][param + '_min'].value), Number(this[model][param + '_max'].value)]);
      }

      // Fill out the remaining values
      for (; i % 4 !== 0; ++i) {
        bounds[c].push([0, 0]);
      }

      return bounds;
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

    truncateString(str) {
      const str_array = String(str).split('.');

      // If there is no decimal, do nothing
      if (str_array.length <= 1) {
        return str;
      } else {
        return str_array[0] + '.' + str_array[1].slice(0, 2);
      }
    }
  };
});
