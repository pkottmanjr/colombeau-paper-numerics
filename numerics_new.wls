#!/usr/bin/env wolframscript
(* ::Package:: *)

(*(*Define smeared Heaviside function*)*)
(*cf. https://en.wikipedia.org/wiki/Bump_function#Examples*)
f[x_] := If[
    x > 0,
    Exp[-1/x],
    0
];
g[x_] := f[x] / (f[x] + f[1-x]);
hEps[eps_, x_] := g[x/eps];
(*Plot[{hEps[1, x], hEps[.2, x], hEps[.05, x]}, {x, -1, 1},  PlotLegends->{1, 0.2, 0.05}]*)

(*(*Set global parameters*)*)
R = 1; (*interface radius*)
c = 1; (*total flux*)
specVol1 = 1; (*value of specific volume for r<R*)
specVol2 = 2; (*value of specific volume for r>R*)
viscosity1 = 1; (*value of viscosity for r<R*)
viscosity2 = 2; (*coefficient for viscosity function from r = R onward*)

(*(*Define transitional profiles for specific volume and viscosity*)*)
specVol[eps_, r_] := specVol1 + specVol2 * hEps[eps, r-R];
viscosity[eps_, r_] := viscosity1 + viscosity2 * hEps[eps, r-R];
(*viscosity[eps_, r_] := viscosity1 + viscosity2 * D[hEps[eps, r-R], r];*)

(*(*Main method for computing the solution inside the interface*)*)
getJump[params_] := Module[
    {solFunc, eps = params[[1]], color = params[[2]]}, 
    sol := NDSolve[
        {
            y'[r] - 2 * c^2/r^5 * specVol[eps, r] + viscosity[eps, r] * D[4 * c/r^3 * specVol[eps, r], r] == 0,
            y[R-eps] == 0
        },
        y, {r, R-eps, R+eps}
        ];
    solFunc = First[y /. sol];
    plt = Plot[solFunc[R+z*eps], {z, -1, 1}, PlotRange->All, PlotStyle->color]
];

(*(*==========================================================*)*)

epsValues = {0.05, 0.01, 0.005, 0.001};
pltColors = {Blue, Red, Green, Magenta};
paramsList = Table[{epsValues[[i]], pltColors[[i]]}, {i, 1, 4}];

plts = Table[getJump[params], {params, paramsList}];
legend = LineLegend[pltColors, epsValues, LegendLabel->"\[Epsilon]", LegendFunction->Panel];
fig = GraphicsRow[{Show[plts], legend}];
Export["fig.png", fig]
