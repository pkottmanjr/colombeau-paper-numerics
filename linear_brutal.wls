#!/usr/bin/env wolframscript
(* ::Package:: *)

(*(*Define smeared Heaviside function*)*)
(*cf. https://en.wikipedia.org/wiki/Bump_function#Examples*)
f[x_] := If[
    x > 0,
    Exp[-1/x],
    0
];
(*g[x_] := f[x] / (f[x] + f[1-x]);*)
g[x_] := f[(x+1)/2] / (f[(x+1)/2] + f[(1-x)/2]);
hEps[eps_, x_] := g[x/eps];
(*Plot[{hEps[1, x], hEps[.2, x], hEps[.05, x]}, {x, -1, 1},  PlotLegends->{1, 0.2, 0.05}]*)

(*(*Set global parameters*)*)
R = 1; (*interface radius*)
c = 1; (*total flux*)
specVol1 = 1; (*value of specific volume for r<R*)
specVol2 = 2; (*value of specific volume for r>R*)
viscosityModel = "linear_brutal";
viscosity1 = 1; (*value of viscosity for r<R*)
viscosity2 = 2; (*coefficient for viscosity function from r = R onward*)

(*(*Define transitional profiles for specific volume and viscosity*)*)
specVol[eps_, r_] := specVol1 + (specVol2 - specVol1) * hEps[eps, r-R];

delta = 0.1;
chiDelta[chiEpsilon_] := If[
    chiEpsilon < 0.5 - delta, 0., If[
        chiEpsilon > 0.5 + delta, 1.,
        1./2. + (chiEpsilon - 0.5)/(2. * delta) + 1./(2.*Pi) * Sin[Pi * (chiEpsilon - 0.5)/delta]
    ]
];
viscosity[eps_, r_] := viscosity1 + (viscosity2 - viscosity1) * chiDelta[hEps[eps, r-R]];

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

(* Set plot parameters *)
epsValues = {0.05, 0.01, 0.005, 0.001};
pltColors = {Blue, Red, Green, Magenta};
paramsList = Table[{epsValues[[i]], pltColors[[i]]}, {i, 1, 4}];

(* Get exact value of jump *)
exactJump[z_] := -4 * c/R^3 * (specVol2 - specVol1) * (viscosity1 + 0.5);

(* Create plots and legend *)
plts = Table[getJump[params], {params, paramsList}];
pltColors = Append[pltColors, {Dashed, Black}];
epsValues = Append[epsValues, "theor. jump"];
plts = Append[plts, Plot[exactJump[z], {z, -1, 1}, PlotStyle->{Dashed, Black}]];
legend = LineLegend[pltColors, epsValues, LegendLabel->"\[CurlyEpsilon]", LegendFunction->Panel]

pltsCombined = Show[
    plts,
    AxesLabel->{
        Style["\[Xi]", FontSize->12, FontColor->Black],
        Style[Superscript[y,"(\[CurlyEpsilon])"], FontSize->12, FontColor->Black]
    }
]
(* fig = GraphicsRow[{Show[plts], legend}]; *)
Export[viscosityModel <> "_fig.pdf", pltsCombined];
Export[viscosityModel <> "_legend.pdf", legend];
