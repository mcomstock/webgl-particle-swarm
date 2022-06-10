/* global define */
define('scripts/interface', [
], function(
) {
  'use strict';

  return class PsoInterface {
    static paramList = ['tr', 'tsi', 'twp', 'td', 'tvp', 'tv1m', 'tv2m', 'twm', 'to', 'xk', 'ucsi', 'uc', 'uv'];

    constructor() {
      PsoInterface.paramList.forEach(param => {
        this[param + '_val'] = document.getElementById(param + '_val');
        this[param + '_min'] = document.getElementById(param + '_min');
        this[param + '_max'] = document.getElementById(param + '_max');
        this[param + '_fit'] = document.getElementById(param + '_fit');

        this[param + '_prev_min'] = -1;
        this[param + '_prev_max'] = -1;

        this[param + '_fit'].addEventListener('change', () => {
          if (this[param + '_fit'].checked) {
            this[param + '_val'].removeAttribute('readonly');
            this[param + '_min'].removeAttribute('readonly');
            this[param + '_max'].removeAttribute('readonly');

            this[param + '_min'].value = this[param + '_prev_min'];
            this[param + '_max'].value = this[param + '_prev_max'];
          } else {
            this[param+'_prev_min'] = this[param + '_min'].value;
            this[param + '_min'].value = this[param + '_val'].value;

            this[param+'_prev_max'] = this[param + '_max'].value;
            this[param + '_max'].value = this[param + '_val'].value;

            this[param + '_val'].setAttribute('readonly', true);
            this[param + '_min'].setAttribute('readonly', true);
            this[param + '_max'].setAttribute('readonly', true);
          }
        });
      });

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

      elem.appendChild(file_in);
      elem.appendChild(cl_label);
      elem.appendChild(cl_in);
      elem.appendChild(plot_button);

      return elem;
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
      const bounds = env.bounds.flat(1);

      PsoInterface.paramList.forEach((param, idx) => {
        if (this[param + '_fit'].checked)
        {
          const [min, max] = bounds[idx];
          this[param + '_min'].value = min;
          this[param + '_max'].value = max;
        }
      });

      this.normalization.value = 1;
    }

    displayResults(bestArr) {
      PsoInterface.paramList.forEach((param, idx) => {
        this[param + '_val'].value = bestArr[idx];
      });
    }

    display_all_params(real_interface)
    {
      var builder = "";
      PsoInterface.paramList.forEach(param => {
        builder = builder.concat(param + ":\t" + this[param + '_val'].value +"\n");
        // builder = builder.concat(param + ":\t" + "test" +"\n");
      });
      prompt("(ctrl+c, Enter) to copy:",builder);
    }

    displayError(error) {
      // this.fit_error.innerHTML = this.truncateString(error);
      this.fit_error.innerHTML = (Number.parseFloat(error).toPrecision(4)).toString();
    }

    getBounds() {
      const bounds = [[], [], [], []];

      PsoInterface.paramList.forEach((param, idx) => {
        bounds[Math.floor(idx/4)].push([Number(this[param + '_min'].value), Number(this[param + '_max'].value)]);
      });

      // Fill out the remaining values
      for (let i = 0; i < 3; ++i) {
        bounds[3].push([0, 0]);
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
