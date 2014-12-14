(* ::Package:: *)

(* Mathematica Package *)

(* Created by the Wolfram Workbench Nov 3, 2010 *)

BeginPackage["Toolbox`"]
(* Exported symbols added here with SymbolName::usage *) 

SubsetQ::usage = "SubsetQ [list,subslist] gives True if sublist is a subset of list, otherwise false."

Logit::usage = "Logit[x] is defined as Log[x/(1-x)]."

Alogit::usage = "Alogit[x] is defined as 1/(1+Exp[-x])."

ToCorrelation::usage = "ToCorrelation[x,b] transfroms a covaraince  matrix x into a correlation matrix. The diagonal elements are standard deviations if b=True, and 1 (by default) otherwsise."

IndexByInterpolation::usage = "IndexByInterpolation[x] constructs an interpolation function f for a list x of a strictly increasing values, where f[z]=i if x[[i]]<=z<x[[i+1]] and f[x[[-1]]]=Length[x]. "

PotentialSRF::usage = "PotentialSRF[data] "

GewekeDiagostic::usage = "GewekeDiagostic[data,c1,c2]"

PosteriorPDFPlot::usage = "PosteriorPDFPlot[data,smooth]"

MCMCDiagostic::usage = "MCMCDiagostic[data,xlabels]"

Begin["`Private`"]
(* Implementation of the package *)

SubsetQ[list_, sublist_] :=
    Complement[sublist, list] === {}

Logit[x_] :=
    Log[x/(1 - x)]

Alogit[x_] :=
    1/(1 + Exp[-x])
    
ToCorrelation[V_?PositiveDefiniteMatrixQ,b_:False] :=
    Module[ {sd, M,i},
        sd = Sqrt[Diagonal[V]];
        M = Transpose[Transpose[V/sd]/sd];
        If[ b,
            Do[M[[i, i]] = sd[[i]], {i, Length[sd]}]
        ];
        M
    ]
    
IndexByInterpolation[data_?(OrderedQ[#] && VectorQ[#, NumericQ] &)] :=
    Module[ {ls, f},
        ls = Transpose[{data, Range[Length[data]] - 1}];
        f = Interpolation[ls, InterpolationOrder -> 0];
        Function[x,
         Round[f[x]] + Switch[Depth[x],
            1, Total[Boole[Thread[Rest[data] == x]]],
            2, Total[Boole[Outer[Equal, x, Rest[data]]], {2}]
           ]
        ]
    ]

Options[PosteriorPDFPlot] = Options[Plot];
PosteriorPDFPlot[data_, smooth_, options : OptionsPattern[]] :=
    Module[ {hls, xmin, xmax, f, g, x,dx},
        hls = HistogramList[data, Automatic, "PDF"];
        {xmin, xmax} = hls[[1, {1, -1}]];
        dx = 10^(-3.) (xmax - xmin);
        hls = Transpose[{Rest[hls[[1]]], hls[[2]]}];
        hls = Join[{{xmin - dx, 0}, {xmin, 0}}, hls, {{xmax + dx, 0}}];
        If[ smooth >= 1,
            hls = Transpose[{MovingAverage[hls[[All, 1]], 2], hls[[2 ;;, 2]]}];
            dx /= 2;
        ];
        f = Interpolation[hls, InterpolationOrder -> smooth];
        g = Plot[f[x], {x, xmin - dx, xmax + dx}, options]
    ]    
  
PotentialSRF[data_] :=
    Module[ {n, m, k, phj, ph, bb, sj2, ww,j},
        {n, m, k} = Dimensions[data];
        phj = Table[Mean[data[[All, j]]], {j, m}];
        ph = Mean[phj];
        bb = Total[(# - ph)^2 & /@ phj] n/(m - 1);
        sj2 = (# - phj)^2 & /@ data;
        sj2 = Table[Total[sj2[[All, j]]], {j, m}]/(n - 1);
        ww = Mean[sj2];
        Sqrt[(n - 1)/n + (bb/ ww) (m + 1)/(n m) ]
    ]
  
OverlapBatchMean[data_] :=
    Module[ {n, m, xmean, datab, vv},
        n = Length[data];
        xmean = Mean[data];
        m = Ceiling[Length[data]/2];
        datab = Partition[data, m, 1];
        vv = (Mean[#] - xmean) & /@ datab;
        vv = (n m)/((n - m + 1) (n - m))  Total[vv^2];
        vv
    ]
  
GewekeDiagostic[data_, c1_, c2_] :=
    Module[ {data1, data2, n1, n2, vv1, vv2, mm1, mm2, zz},
        data1 = data[[;; Ceiling[c1 Length[data]]]];
        data2 = data[[-Floor[c2 Length[data]] ;;]];
        n1 = Length[data1];
        n2 = Length[data2];
        vv1 = OverlapBatchMean[data1];
        vv2 = OverlapBatchMean[data2];
        mm1 = Mean[data1];
        mm2 = Mean[data2];
        zz = (mm1 - mm2)/Sqrt[vv1/n1 + vv2/n2];
        {zz, 2 CDF[NormalDistribution[], -Abs[zz]]}
    ]  
    
MCMCDiagostic[data_, labels_] :=
    Module[ {data2,psrf, zscore, pvalue,i,half,ls,g},
        data2 = data[[-Floor[Length[data]/2] ;;]];
        Print["-----------------------Start----------------------------"];
        Print["Geweke Convergence Diagnositics for the 2nd half chains:"];
        Do[
          {zscore, pvalue} = GewekeDiagostic[data2[[All, i]], 0.1, 0.5];
          Print["Chain " <> ToString[i] <> ": \n", TableForm[{zscore, pvalue},
            TableHeadings -> {{"Z-Score", "P-Value"}, labels}]], {i, 
           Dimensions[data2][[2]]}];
        Print["--------------------------------------------------------"];
        Print["Brooks, Gelman, and Rubin Convergence Diagnositics:"];
        psrf = PotentialSRF[data2];
        Print["Potential Scale Reduction Factors for the 2nd half chains:\n", 
          TableForm[{psrf}, TableHeadings -> {{"PSRF"}, labels}]];
        Print["Trace plot of PSRF:"];
        half = Range[50, Length[data], 50];
        half = Thread[Span[half/2 + 1, half]];
        ls = PotentialSRF[data[[#]]] & /@ half;
        g = ListLinePlot[Transpose[ls], Joined -> True, Frame -> True, 
        PlotRange -> All, PlotStyle -> Thick, GridLines -> Automatic, 
        DataRange -> half[[{1, -1}, 2]], 
        FrameLabel -> {"iteration", "psrf"}, 
        PlotLegends -> 
        LineLegend[labels, 
        LegendFunction -> (Framed[#, RoundingRadius -> 5] &), 
        LegendMargins -> 5]];
        Print[Show[g, ImageSize -> 600]];
        Print["----------------------The End--------------------------"];
    ]    
      
      
End[]

EndPackage[]

