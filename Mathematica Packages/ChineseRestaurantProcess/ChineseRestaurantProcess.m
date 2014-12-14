(* Mathematica Package *)

(* Created by the Wolfram Workbench Jul 3, 2010 *)

BeginPackage["ChineseRestaurantProcess`"]
(* Exported symbols added here with SymbolName::usage *) 

RandomChineseRestaurant::usage = "RandomChineseRestaurant[n,\[Alpha],\[Theta],k] generates k instances of the Chinese restaurant process with discount parameter \[Alpha] and strength parameter \[Beta]. By default, k=1.";

ReturnType::usage = "ReturnType is an option specifying the form of returning results. In case of 'Configuration', return a list of each instance {customers at table 1,customers at table 2, ...}, and in case of 'Frequency' (default), return a list of each instance {{number of customers in a table, number of such tables},...}. "

ConfigToFreq::usage = "ConfigToFreq[x] transforms Configuration x from RandomChineseRestaurant into Frequency form."

LogPDFCRP::usage = "LogPDFCRP[x,n,\[Alpha],\[Theta]] gives the log probability density function for the x generated from RandomChineseRestaurant."

LogPDFnTable::usage = "LogPDFnTable[k, n, \[Theta]] gives the log probability for the Chinese restaurant process (\[Alpha]=0) with parameters n and \[Theta]."

Begin["`Private`"]

RandomChineseRestaurant::badab = "Either \[Alpha]<0 and \[Theta]=L \[Alpha] for L=1,2,..; or that 0\<=\[Alpha]<=1 and \[Theta]>=\[Alpha]";

RandomChineseRestaurant::badtype = "The option value of ReturnType must be either 'Configuration' or 'Frequency'.";

Options[RandomChineseRestaurant] = {ReturnType -> "Configuration"};

RandomChineseRestaurant[n_ /; IntegerQ[n] && n >= 1, alpha_, theta_,k_ /; IntegerQ[k] && k >= 1, OptionsPattern[]] :=
    Module[ {r, x, res, restype},
        restype = OptionValue[ReturnType];
        If[ ! MemberQ[{"Configuration", "Frequency"}, restype],
            Message[RandomChineseRestaurant::badtype];
            Return[$Failed]
        ];
        If[ alpha < 0,
            r = theta/alpha;
            If[ Abs[Round[r] - r] > 10^(-10),
                Message[RandomChineseRestaurant::badab];
                Return[$Failed]
            ],
            If[ ! (alpha <= 1 && theta > -alpha),
                Message[RandomChineseRestaurant::badab];
                Return[$Failed]
            ]
        ];
        x = Table[{1, {{1}}}, {k}];
        res = Nest[AddOneCustomer[#, alpha, theta] &, #, n - 1] & /@ x;
        Switch[restype,
         "Configuration", Sort[#]&/@res[[All, 2]],
         "Frequency", 
         Map[ConfigToFreq,res[[All,2]]]                 
         ]
    ];
    
RandomChineseRestaurant[n_ /; IntegerQ[n] && n >= 1, alpha_, theta_,opts : OptionsPattern[]] :=
    First[RandomChineseRestaurant[n, alpha, theta, 1, opts]];    

(*
npeop: number of people; 
alpha: discount parameter;
theta: strength parameter;
ntab: number of tables;
config: assignment of customers to tables;
pnew: probability that new customer occupies a new table;
pocc: probability vector that new customer is assigned to occupied \
tables;
*)
AddOneCustomer[x_, alpha_, theta_] :=
    Module[ {npeop, config, ntab, pnew, pocc, p, tab},
        {npeop, config} = x;
        ntab = Length[config];
        pnew = (theta + ntab alpha)/(npeop + theta);
        pocc = ((Length[#] & /@ config) - alpha)/(npeop + theta);
        p = Append[pocc, pnew];
        tab = RandomChoice[p -> Range[ntab + 1]];
        If[ tab == ntab + 1,
            AppendTo[config, {npeop + 1}],
            AppendTo[config[[tab]], npeop + 1]
        ];
        {npeop + 1, config}
    ];
    
ConfigToFreq[config_] :=
    Module[ {n, temp, res},
        n = Max[config];
        temp = Tally[Map[Length, config]];
        res = ConstantArray[0, n];
        res[[temp[[All, 1]]]] = temp[[All, 2]];
        res
    ];     
    
LogPDFCRP[alist : {{{__} ..} ..}, n_Integer, alpha_?NumericQ, 
  theta_?NumericQ] := Module[{len,k},
  len = Map[Length, alist, {2}];
  k = Length[#] & /@ len;
  N[If[alpha == 0,
    LogGamma[theta] + k Log[theta] - LogGamma[theta + n] + 
     Total[LogGamma[len], {2}],
    LogGamma[theta] + k Log[theta] - LogGamma[theta + n] + 
     LogGamma[theta/alpha + k] - LogGamma[theta/alpha] + 
     Total[LogGamma[len - alpha], {2}] - LogGamma[1 - alpha]
    ]]
  ]
LogPDFCRP[a : {{___} ..}, n_Integer, alpha_?NumericQ, 
  theta_?NumericQ] := First[LogPDFCRP[{a}, n, 0, theta]]
  
LogPDFCRP[a_ /; VectorQ[a, NumericQ], n_Integer, alpha_?NumericQ, 
   theta_?NumericQ] :=
    Module[ {b, k, res},
        b = Select[Transpose[{Range[n], a}], #[[2]] > 0 &];
        k = Total[b[[All, 2]]];        
        res = If[ alpha == 0,
                  LogGamma[N[theta]] + k Log[N[theta]] - LogGamma[N[theta + n]] + 
                   Total[b[[All,2]] LogGamma[b[[All,1]]]],
                  LogGamma[N[theta]] + k Log[N[alpha]] - LogGamma[N[theta + n]] + 
                   LogGamma[N[theta/alpha + k]] - LogGamma[N[theta/alpha]] + 
                   Total[b[[All, 2]] (LogGamma[b[[All, 1]] - alpha] - 
                       LogGamma[1 - alpha])]
              ];
        res += 
         LogGamma[N[n + 1]] - 
          Total[b[[All, 2]] LogGamma[N[b[[All, 1]] + 1]] + 
            LogGamma[N[b[[All, 2]] + 1]]]
    ]/;Length[a]===n
 
LogPDFnTable[k_?IntegerQ, n_?IntegerQ, theta_?NumericQ] :=
    N[Log[Abs[StirlingS1[n, k] ]] + k Log[theta] + LogGamma[theta] - 
     LogGamma[N[theta + n]]];
        
End[]

EndPackage[]


