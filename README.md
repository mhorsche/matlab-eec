# EQUIVALENT-ELECTRIC-CIRCUIT (EEC) model class for MATLAB®

`eec` is a MATLAB® class to easily define and model equivalent-electric-circuits, mainly for the purpose to model the electric behaviour of batteries.

Installation
============

1. Extract the ZIP file (or clone the git repository) somewhere you can easily reach it. 
2. Add the `src/` folder to your path in MATLAB/Octave: e.g. 
    - using the "Set Path" dialog in MATLAB, or 
    - by running the `addpath` function from your command window or `startup` script.

Usage
=====

The basic usage is to define an `eec` model with an individual set of `OCV`, `R`, `L`, `C`, `RC` and `RL` elements. Exemplarily the command `eec('OCV+R')` generates a simple model containing a constant open-circuit-voltage in series with a resistor.

After defining a model, an `AC` or `DC` analysis can be performed. Therefor a time and current or frequency input vector is necessary. To fit the model parameter to a given reference also the voltage response or the complex impedance is required.

Most functionality is described in the given example which can be run by the static function `eec.example()`. In this example following commands are executed:

1. Generate a new model containing `OCV` + `R` `L` and three `RC` elements:

    ```MATLAB
    eecmodel = eec('OCV,R,L,RC,RC,RC',...
    ...%                                        RC1             RC2             RC3
    ...%               OCV     R       L        R    tau  init  R     tau init  R    tau init
      'InitialValues',[3.6500  0.0300  500e-9   5e-3 1e-3 0     10e-3 0.1 0     0.02 10  0 ]);
    ```

2. Run a simple `DC` analysis for a constant current pulse of `I = -2 A` (discharge) and `t = 10 sec`:

    ```MATLAB
    dc_res = dc(eecmodel,linspace(0,10,101),-2);
    ```

3. Generate a second model containing the same elements but different initial values:

    ```MATLAB
    eecmodel_fitdc = eec('OCV,R,L,RC,RC,RC',...
      InitialValues',[3.5000  0.0200  300e-9   4.0e-3 500e-6 0 15e-3 300e-3 0  25e-3 15 0 ]);
    ```

4. Run fitting algorithm in time domain to match the reference values, also show results and manual fitting frontend:

    ```MATLAB
    dcfit(eecmodel_fitdc,dc_res.t_s,dc_res.I_A,dc_res.U_V,...
      'Method','analytical',...
    ...%               OCV     R       L        RC1             RC2             RC3
      'LowerBoundary',[2.500   0       0        0     0     0   0     0     0   0     0     NaN ],...
      'UpperBoundary',[4.200   Inf     Inf      10e-3 1.0   0   0.05  10    0   0.1   100   NaN ]);
    eecmodel_fitdc.timeseries(dc_res.t_s,dc_res.U_V,dc_res.t_s,dc_res.I_A);
    ```

![Time domain](https://github.com/mhorsche/matlab-eec/blob/master/docs/images/TimeSeries.png)

5. Run a simple `AC` analysis for a frequency range of `f = 0.01 Hz ... 10 kHz`:

    ```MATLAB
    ac_res = ac(eecmodel,logspace(-2,4,101));
    ```

6. Generate a third model containing the same elements but the fitted values from the DC model:

    ```MATLAB
    eecmodel_fitac = eec('OCV,R,L,RC,RC,RC',...
      'InitialValues',[eecmodel_fitdc.Elements.value]);
    ```

7. Run fitting algorithm in frequency domain to match the reference values, also show results and manual fitting frontend:

    ```MATLAB
    acfit(eecmodel_fitac,ac_res.f_Hz,ac_res.ZRe_Ohm,ac_res.ZIm_Ohm,...
      'Method','analytical',...
    ...%               OCV     R       L        RC1             RC2             RC3
      'LowerBoundary',[NaN     0       0        0     0     0   0     0     0   0     0    NaN ],...
      'UpperBoundary',[NaN     Inf     Inf      10e-3 1.0   0   0.05  10    0   0.1   100  NaN ]);
    eecmodel_fitac.nyquist(ac_res.f_Hz,ac_res.ZRe_Ohm,ac_res.ZIm_Ohm);
    ```

![Frequency domain](https://github.com/mhorsche/matlab-eec/blob/master/docs/images/Nyquist.png)

Remarks
-------
Most functions accept numerous name-value arguments, inspect the help:

```MATLAB
doc eec
```

To run the `DC` analysis in spice mode you need to specify the path to the LTspice® executable. To read the results back into MATLAB, the plugin [LTspice2Matlab](https://github.com/PeterFeicht/ltspice2matlab) is used.

More information
================
If you experience bugs or would like to request a feature, please visit the [issue tracker](https://github.com/mhorsche/matlab-eec/issues).