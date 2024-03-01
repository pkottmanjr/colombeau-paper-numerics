# colombeau-paper-numerics
Wolfram Script files performing numerical experiments for the Colombeau paper. The scripts in this repository should reproduce the plots in Fig. 3 in the paper.

ATM the solution is suboptimal - one script per plot. To create plots, run

```powershell
wolframscript -file linear.wls
wolframscript -file power_law_illustration.wls
wolframscript -file shear_thickening.wls
wolframscript -file shear_thinning.wls
```

A better version (not implemented yet) is one where one script iterates over viscosity models to create all the plots. For this, run

```powershell
wolframscript -file .\numerics_new.wls
```