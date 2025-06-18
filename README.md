# WebGL implementation of particle swarm optimization

## How to run

To run the code, simply run any local HTTP server and go to the URL it provides. For example, I use
Python's `http.server` by running
```bash
python -m http.server
```
and navigating to the address it prints after starting.

## Interface

The following section contains a brief description of the provided behavior and intended workflow
for the program.

### Specifying data to fit

Fitting a model to a data set requires adding the data and setting a few options. The section of the
user interface focused on data input appears directly below the graph window, beginning with the
"Normalize?" button and ending with the "Add data" and "Remove data" buttons.

Currently, three types of data are available: voltage time series data, calcium time series data,
and APD data. Any combination of data types can be fit simultaneously by adding additional data
rows.

#### Time series data (voltage and calcium)

Time series data must be provided as newline-delimited files with each line containing a single
number (which can be parsed by JavaScript) representing the voltage or calcium value. Normalization
may be enabled by checking the checkbox labeled "Normalize?", in which case the data will be
normalized to the range specified by the fields "Normalization minimum" and "Normalization
maximum". Different models have different default normalization settings, which will be set
automatically when the model is selected.

The row of fields directly below the "Normalize?" button and above the data rows applies to all
input data. The "Number of stimuli" field is used to tell the program how many pacing periods are
included in each data set. The value should be an integer; for example, even if a single action
potential was recorded for 800 milliseconds with a cycle length of 1 second, the number of stimuli
should be set to 1. The "Pre-recording stimuli" field is used to set the number of cycle lengths
that should occur in the simulation before comparison with the data. Larger numbers will lead to
longer simulations, but will reduce transient effects due to the initial conditions of the
simulation. The "Sample interval" field is used to set the time resolution of the data, which
currently must be the same for all data files. High resolution data should ideally be downsampled to
a sample interval of around 1 millisecond, as excessively large data files may slow the program and
possibly prevent it from running.

Data files are added using the "Browse" button that appears on each data row (above the "Add data"
button). In addition to the data file, the pacing period/pacing cycle length corresponding to the
data must be entered in the "Cycle length" field. If multiple data sets are being fit, the relative
weight of each data set may be adjusted using the "Fitting weight" field, which defaults to a value
of 1. The fitting weight is used as the coefficient for the error term PSO cost function for the
corresponding data set, with the total cost function value being the sum of each component. Note
that the error terms are normalized by the number of data points but not by the magnitude of the
data, so manually adjusting the weights may be particularly important when fitting voltage and
calcium at the same time. The data type may be changed between voltage, calcium, and APD data using
the "Data type" drop-down menu.

The "Plot" button at the far right of the data row can be used to visualize the data before a fit is
generated. After a fit has been generated, the functionality of the plot button changes to plot both
the data from that row and the fitted model paced to the same cycle length for comparison.

To fit to more than one data set simultaneously, use the "Add data" button to create an additional
data row. Data rows may also be removed using the "Remove data" button.

#### Action potential duration (APD) data

Adding APD data to fit is a similar process to adding time series data, with a few key
differences. To change a data row to accept APD data, the "Data type" field must be set to
"APD". The data row will then change to accept APD input. Rather than uploading a file, the APDs are
entered as comma-separated numbers (in milliseconds) in the "APDs" field. The "APD threshold" field
contains the voltage value where the APDs are measured. The "Cycle length" and "Fitting weight" have
the same behavior for APD data as for time series data.

The "Plot" button will not show the APD data without a fit being generated. Once a fit has been
generated, it can be used to compare the model fit with the data: it will plot the model paced at
the cycle length for the data set and display the APDs as vertical lines aligned with the start of
each simulated APD.

### Choosing a model and parameter ranges

The fields to the right of the graph window contain the model selection and settings for parameter
fitting. The "Model" drop-down is used to select the model to fit, which will automatically populate
the fields below with the relevant entries for the parameters associated with the selected
model. Each of the following rows contains the settings relevant to a particular parameter of the
model, with one row per parameter available for fitting with PSO.

The meanings of the columns in the parameter section are:

- "Parameter": The name of the parameter.
- "Value": The value to be used for the parameter. If the parameter is to be fit, this value will be
  replaced be the result of the fit. If the parameter is not selected to be fit, then this value
  will be used for the parameter throughout the fitting process. Once a fitting is complete, the
  best value found for each parameter will be displayed in this field.
- "Minimum": The lower bound for the parameter during fitting.
- "Maximum": The upper bound for the parameter during fitting.
- "Fit?": Checked if the parameter is to be fit during fitting. If unchecked, the parameter will be
  set to the value in its value field during the fitting process.

The "Fit all" and "Fit none" buttons provide a convenient option to check or uncheck (respectively)
the "Fit?" checkboxes for every parameter. The "Plot from values" button will display the result of
the currently selected parameter values in the "Value" fields to the graph window. The "Reset
bounds" button will reset the parameter bounds to their default values.

### Stimulus settings

The stimulus current can be set to either a square or biphasic waveform. The stimulus settings are
located below the section for inputting data sets, directly below the "Add data" and "Remove data"
buttons. The stimulus is set to a biphasic waveform if the "Biphasic stimulus" checkbox is checked,
and a square waveform otherwise. The purpose of the biphasic stimulus is to mimic the current
experienced through diffusive coupling when only simulating a single cardiac cell.

When the square stimulus is selected, only two fields are present below the "Biphasic stimulus"
checkbox for configuration. The "Stimulus duration" field sets the amount of time in milliseconds
that the stimulus should remain active. The "Stimulus magnitude" field determines the magnitude of
the stimulus. The units of magnitude match the current units of the model, so the appropriate range
of values depends on the model being fit.

When the biphasic stimulus is selected, the waveform of the stimulus is computed according to the
equation $I_\text{stim}(t) = -I_\text{mag} \times (t/a - b/(1 + (t/a - c)^4))$. The stimulus
configuration fields below the "Biphasic stimulus" checkbox have the following effects:

- "Stimulus duration": The total time the biphasic stimulus is applied in milliseconds.
- "Stimulus magnitude": $I_\text{mag}$ in the biphasic stimulus equation. As with the square
  waveform, the units of the magnitude are the same as the current units for the selected model.
- "Stimulus offset 1": $b$ in the biphasic stimulus equation.
- "Stimulus offset 2": $c$ in the biphasic stimulus equation.
- "Stimulus timescale": $a$ in the biphasic stimulus equation.

### Setting PSO hyperparameters

All hyperparameters for PSO are initialized to default values, but they may be adjusted through the
user interface. The hyperparameter settings are located directly below the section for configuring
the stimulus.

The "Particle count" option allows the choice of a few different numbers of particles. In general,
larger numbers of particles are expected to produce better results but require more processing power
and running time. Depending on the graphics hardware available to the browser and the specific
browser settings, choosing too large a particle count may cause the program to fail; if this
happens, it is recommended to refresh the page and try again with a lower particle count.

The "Iteration count" field sets the number of PSO iterations to perform before selecting the
lowest-error parameterization as the best fit. The running time of the program is expected to
increase proportionally with the iteration count when other factors remain unchanged, that is, each
iteration is expected to take the same amount of time. While too few iterations will result in a
poor quality fit, increasing the number of iterations will not necessarily lead to a noticeable
improvement. The fit error plot below the hyperparameter section (described below) may be useful for
choosing a reasonable number of iterations for fitting a specific set of data.

The "Local uniform maximum" and "Global uniform maximum" fields set the relative influences of the
local best and global best position, respectively, in the PSO algorithm. The "Constriction
coefficient" field scales the velocity update step to promote convergence in the PSO
algorithm. These parameters and their meanings are described in further detail in Section 2.3.2 of
the following reference, with "Local uniform maximum" set to $\phi_1$, "Global uniform maximum" set
to $\phi_2$, and "Constriction coefficient" set to $\chi$: Loewe, Axel, et al. "Parameter estimation
of ion current formulations requires hybrid optimization approach to be both accurate and reliable."
Frontiers in bioengineering and biotechnology 3 (2016): 209.

### Running the program and saving the results

The "Run" button at the top of the initiates the model fitting process when clicked. As the program
runs, the current iteration of the PSO algorithm is displayed along with the total number of
iterations to indicate the current progress. Once the iteration process is completed and the fit has
been found, the result of the fit compared with one of the data sets is displayed in the graph
window. The input data is displayed in black and the model fit is shown in red. The data row used
for comparison can be changed to a specific data set by clicking the "Plot" button to the right of
the corresponding data row.

The parameter values resulting from the fit can be saved using the "Save parameters" button, while a
more complete set of information including the output parameter values, hyperparameter settings and
information, about the input data can be saved using the "Save run details" button. These files are
saved as downloads in the JSON format. The values of the parameters found for the best fit are also
displayed in the "Value" field for each parameter up to three digits of precision.

In addition to the plot of the model fit and the data, an additional plot at the bottom of the page
visualizes the best fit error at each iteration of the PSO algorithm. If the error still seems to
decrease significantly with further iterations even for the last few iterations, a better fit may be
achieved by increasing the number of iterations. On the other hand, if many of the final iterations
do not improve the best fit error, reducing the number of iterations may lead to similar fit quality
while reducing the amount of computation and time required to generate a fit.

## Potential failure cases

Although support for the WebGL API is available in most modern web browsers, correctly implemented
WebGL programs are still not guaranteed to work. The WebGL API (and the browser in general) has
limited access to graphics hardware, which is ultimately controlled by the operating system, often
without providing the users direct means to configure settings related to the utilization of
graphics hardware.

WebGL programs that place a large burden on graphics hardware typically face several
limitations. The first is that the program must use whatever GPU the operating system has assigned
to render graphics for the browser. A common problem that arises from this arrangement occurs when a
computer has multiple GPUs, as the correct GPU must be assigned to the browser using either
operating system or browser settings. Another common issue arises when the load placed on the GPU by
the WebGL program is too high---the definition of too high is again determined by either the
operating system or the browser and is not controlled by the program itself. In such cases, the
program will lose access to the WebGL context used to access the GPU, and the program will
ultimately fail to run successfully. The exact manifestation of this issue can take several forms
depending on both the browser and operating system, such as an error message indicating that the
context was lost or the browser freezing or even crashing entirely. In cases where this issue is
encountered, a possible solution is to attempt to reduce the load on the GPU by stopping other
programs that use system resources while the program is running or to disable power-saving settings
(for example, by connecting a laptop to a power source). Another possible solution is to reduce the
number of particles, the complexity of the model being fit, or the simulation time required to fit
the data in the settings described above if more powerful hardware is not available. Increasing or
reducing the iteration count usually does not affect this issue, as the maximum resource usage
achieved by the program depends on the amount of work per iteration, so the program is likely to
succeed once the first few iterations have completed. Yet another potential cause of failure is due
to issues in the implementation of graphics drivers. These issues are best avoided by keeping the
system software drivers up to date, and are generally specific to the combination of operating
system and hardware.

Most operating systems and browsers are configured so that graphical applications will not run when
they are not visible on-screen. Therefore, running a PSO fit with the browser not in focus or
viewing a different tab will usually result in the program running extremely slowly or stopping
entirely. In some cases, it may be possible to configure the browser to run graphical applications
when they are not in focus, although this behavior is generally not desirable, particularly if many
browser tabs are open at once.
