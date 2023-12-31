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
viscosity1 = 1; (*value of viscosity for r<R*)
viscosity2 = 2; (*coefficient for viscosity function from r = R onward*)

(*(*Define transitional profiles for specific volume and viscosity*)*)
specVol[eps_, r_] := specVol1 + specVol2 * hEps[eps, r-R];
(* viscosityModel = "linear";
viscosity[eps_, r_] := viscosity1 + viscosity2 * hEps[eps, r-R]; *)
(* viscosityModel = "polynomial"; q = 2;
viscosity[eps_, r_] := viscosity1 + viscosity2 * hEps[eps, r-R]^q; *)
(* viscosityModel = "power-law";
viscosity[eps_, r_] := viscosity1 + viscosity2 * D[hEps[eps, r-R], r]; *)
viscosityModel = "full-power-law"; q = 1;
viscosity[eps_, r_] := viscosity1 * ((D[c/r^2 * specVol[eps, r], r])^2 + 2 * c^2/r^6 * specVol[eps, r]^2)^(q/2);
(* Switch[
    viscosityModel,
    "linear", viscosity[eps_, r_] := viscosity1 + viscosity2 * hEps[eps, r-R];
    "polynomial", viscosity[eps_, r_] := viscosity1 + viscosity2 * hEps[eps, r-R]^q;
    "power-law", viscosity[eps_, r_] := viscosity1 + viscosity2 * D[hEps[eps, r-R], r];
    "full-power-law", viscosity[eps_, r_] := viscosity1 * ((D[c/r^2 * specVol[eps, r], r])^2 + 2 * c^2/r^6 * specVol[eps, r]^2)^((q-2)/2);
]; *)

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
    (* If[
        viscosityModel == "power-law",
        plt = LogPlot[Abs[solFunc[R+z*eps]], {z, -1, 1}, PlotRange->All, PlotStyle->color],
        plt = Plot[solFunc[R+z*eps], {z, -1, 1}, PlotRange->All, PlotStyle->color]
    ] *)
    plt = Plot[solFunc[R+z*eps], {z, -1, 1}, PlotRange->All, PlotStyle->color]
];

(*(*==========================================================*)*)

(* Set plot parameters *)
epsValues = {0.05, 0.01, 0.005, 0.001};
pltColors = {Blue, Red, Green, Magenta};
paramsList = Table[{epsValues[[i]], pltColors[[i]]}, {i, 1, 4}];

(* Get exact value of jump *)
exactJump[z_] := Switch[
    viscosityModel,
    "linear", -4 * c/R^3 * specVol2 * (viscosity1 + viscosity2/2),
    "polynomial", -4 * c/R^3 * specVol2 * (viscosity1 + viscosity2/(q+1))
];

(* Create plots and legend *)
plts = Table[getJump[params], {params, paramsList}];
If[
    viscosityModel != "power-law" && viscosityModel != "full-power-law",
    (
        pltColors = Append[pltColors, {Dashed, Black}];
        epsValues = Append[epsValues, "theor. jump"];
        plts = Append[plts, Plot[exactJump[z], {z, -1, 1}, PlotStyle->{Dashed, Black}]];
    )
];
legend = LineLegend[pltColors, epsValues, LegendLabel->"\[CurlyEpsilon]", LegendFunction->Panel]

pltsCombined = Show[
    plts,
    (* PlotLabel->Style[
        viscosityModel <> " viscosity profile" <> If[viscosityModel == "polynomial", " with q = " <> ToString[q], ""],
        FontSize->16,
        FontColor->Black
    ], *)
    AxesLabel->{
        Style["\[Xi]", FontSize->12, FontColor->Black],
        Style[Superscript[y,"(\[CurlyEpsilon])"], FontSize->12, FontColor->Black]
    }
]
(* fig = GraphicsRow[{Show[plts], legend}]; *)
Export[viscosityModel <> If[viscosityModel == "polynomial", "_q_" <> ToString[q], ""] <> "_fig.pdf", pltsCombined];
Export[viscosityModel <> "_legend.pdf", legend];