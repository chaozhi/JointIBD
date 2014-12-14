(* Mathematica Package *)

BeginPackage["JointIBD`Likelihood`"]
(* Exported symbols added here with SymbolName::usage *)  

IBDTransition::usage = "IBDTransition  "

IBDTransitionProb::usage = "IBDTransitionProb  "

IBDLogLikelihood::usage = "IBDLogLikelihood  "

IBDFoldTransition::usage = "IBDFoldTransition  "

IBDTransitionType1::usage = "IBDTransitionType1  "

IBDTransitionType2::usage = "IBDTransitionType2  "

Begin["`Private`"] (* Begin Private Context *) 


IBDTransition[ibd_, n_, theta_, rho_, dd_] :=
    Module[ {p},
        p = 1 - (1 - rho)^dd;
        If[ RandomReal[] < p,
            IBDTransition[ibd, n, theta],
            ibd
        ]
    ];
   
IBDTransition[ibd_, n_, theta_] :=
    Module[ {z = ibd, i1, del,res},
        i1 = RandomChoice[Append[Table[1,{n}],theta]->Range[n+1]];
        del = RandomInteger[{1, n + 1}];
        res = If[ i1==del||del==n+1,
                  z,
                  z = DeleteCases[DeleteCases[z, del, {2}], {}];
                  If[ i1==n+1,
                      Append[z,{del}],
                      Replace[z, i1->Sequence[i1, del], {2}]
                  ]
              ];
        Map[Sort,res,{0,1}]
    ]    
    
(*NestList[IBDTransition[#, n, theta] &,ibd, nfold]*)    
IBDFoldTransition[ibd_, n_, theta_, nfold_Integer] :=
    Module[ {i1list, dellist, f},
        i1list = RandomChoice[Append[Table[1, {n}], theta] -> Range[n + 1], nfold];
        dellist = RandomInteger[{1, n + 1}, nfold];
        f = Function[{x, id},
          Module[ {z, i1, del},
              {i1, del} = id;
              If[ i1 == del || del == n + 1,
                  x,
                  z = DeleteCases[DeleteCases[x, del, {2}], {}];
                  If[ i1 == n + 1,
                      Append[z, {del}],
                      Replace[z, i1 -> Sequence[i1, del], {2}]
                  ]
              ]
          ]
          ];
        FoldList[f, ibd, Transpose[{i1list, dellist}]]
    ]    
  
IBDFoldTransition[ibd_, n_, theta_, ulist_List] :=
    Module[ {i1list, dellist, f, nfold},
        nfold = Length[ulist];
        i1list = RandomChoice[Append[Table[1, {n}], theta] -> Range[n + 1], nfold];
        dellist = RandomInteger[{1, n + 1}, nfold];
        f = Function[{x, id},
          Module[ {z, i1, del, u},
              {i1, del, u} = id;
              If[ i1 == del || del == n + 1 || u == 0,
                  x,
                  z = DeleteCases[DeleteCases[x, del, {2}], {}];
                  If[ i1 == n + 1,
                      Append[z, {del}],
                      Replace[z, i1 -> Sequence[i1, del], {2}]
                  ]
              ]
          ]
          ];
        FoldList[f, ibd, Transpose[{i1list, dellist, ulist}]]
    ]
    
IBDTransition[ibdlist_, n_, theta_, rho_, xstart_Integer,xend_Integer] :=
    Module[ {nb},
        nb = RandomInteger[BinomialDistribution[xend - xstart, rho],Length[ibdlist]];
        MapThread[Transpose[{Prepend[Sort[RandomSample[Range[xstart + 1, xend], #1]], xstart],
             IBDFoldTransition[#2, n, theta, #1]}] &, {nb, ibdlist}]
    ] /; xend > xstart
    
IBDTransition[{xstart_Integer, ibdstart_}, n_, theta_, rho_, xend_Integer] :=
    Module[ {nb},
        nb = RandomInteger[BinomialDistribution[xend - xstart, rho]];
        Transpose[{Prepend[Sort[RandomSample[Range[xstart + 1, xend], nb]],xstart],
            IBDFoldTransition[ibdstart, n, theta, nb]}]
    ] /; xend > xstart
    
IBDLogLikelihood[y_?VectorQ, ibd_?(Depth[#] == 3 &), af_?NumericQ, epsilon_] :=
    First[IBDLogLikelihood[{y},{ibd},{af},epsilon]];   
IBDLogLikelihood[ydata_?MatrixQ, zibd_?(Depth[#] == 4 &), af_?VectorQ, epsilon_] :=
    Module[ {yls, n0, n1, i,temp, lhs, rhs, rule},
        temp = Transpose[{ydata, zibd}];
        temp = SplitBy[temp, #[[2]] &];
        yls = Flatten[Table[Transpose[temp[[i, All, 1, #]] & /@ temp[[i, 1, 2]]], {i,Length[temp]}], 1];
        (*yls=Table[ydata[[i,#]]&/@zibd[[i]],{i,Length[ydata]}];*)
        yls = Map[Tally, yls, {1}];
        lhs = Union[Flatten[yls[[All, All, 1]], 1]];
        n1 = Total[lhs, {2}];
        n0 = Length[#] & /@ lhs - n1;
        rhs = Transpose[{(1 - epsilon)^n0 (epsilon)^n1, (1 - epsilon)^n1 (epsilon)^n0}];
        rule = Dispatch[Thread[lhs -> rhs]];
        yls[[All, All, 1]] = Replace[yls[[All, All, 1]], rule, {2}];
        Total[yls[[All, All, 2]] Log[af yls[[All, All, 1, 1]] + (1 - af)  yls[[All, All, 1, 2]]], {2}]
    ] /; Length[ydata] == Length[zibd] == Length[af]
       
IBDLogLikelihood3[ydata_?MatrixQ, zibd_?(Depth[#] == 4 &),af_?VectorQ, epsilon_] :=
    Module[ {yls, n0, n1, temp,i},
        temp = Transpose[{ydata, zibd}];
        temp = SplitBy[temp, #[[2]] &];
        yls = Flatten[Table[Transpose[temp[[i, All, 1, #]] & /@ temp[[i, 1, 2]]], {i,Length[temp]}], 1];
        (*yls = Table[ydata[[i, #]] & /@ zibd[[i]], {i, Length[ydata]}];*)
        yls = Map[Tally, yls, {1}];
        n1 = Total[yls[[All, All, 1]], {3}];
        n0 = Map[Length, yls[[All, All, 1]], {2}] - n1;
        (*assume at most two alleles*)
        MapThread[Dot,{Log[af (1 - epsilon)^n0 (epsilon)^n1 + (1 - af) (1 - epsilon)^n1 (epsilon)^n0], yls[[All, All, 2]]}]
    ] /; Length[ydata] == Length[zibd] == Length[af]
        
IBDLogLikelihood2[ydata_?MatrixQ, zibd_?(Depth[#] == 4 &),af_?VectorQ, epsilon_] :=
    Module[ {q0, q1, yls, n1, n0, len, i,res},
        yls = Table[ydata[[i, #]] & /@ zibd[[i]], {i, Length[ydata]}];
        n1 = Map[Total, yls, {2}];
        n0 = Map[Length, yls, {2}] - n1;
        len = Length[#] & /@ zibd;
        q0 = MapThread[ConstantArray, {af, len}];
        q1 = 1 - q0;
        (*res = q0 (1 - epsilon q1)^n0 (epsilon q1)^n1 +q1 (1 - epsilon q0)^n1 (epsilon q0)^n0;*)
        (*assume at most two alleles*)
        res = q0 (1 - epsilon)^n0 (epsilon)^n1 +q1 (1 - epsilon)^n1 (epsilon)^n0;
        Total[Log[res], {2}]
    ]/; Length[ydata] == Length[zibd] == Length[af]     

IBDTransitionType1[{z1_, z2_}] :=
    Module[ {z1z2, z1c, z2c, id, len, min,A},
        z1z2 = Intersection[z1, z2];
        z1c = Complement[z1, z1z2];
        z2c = Complement[z2, z1z2];
        id = Length[z1c] 10 + Length[z2c];
        Switch[id,
         0,
         {0, 0, 0},
         12,
         len = Length[#] & /@ z2c;
         min = Min[len];
         Switch[ min,
              1,{1, id, Join[z1c, {{}, Flatten[Pick[z2c, len, min]]}]},
             _, {-1, -1, -1}
         ],
         21,
         len = Length[#] & /@ z1c;
         Switch[len,  
          {1, 1}, {1,id, Append[z1c, Flatten[z1c]]},
          {1, _}, {1,id, Append[z1c, z1c[[1]]]},
          {_, 1}, {1,id, Append[Reverse[z1c], z1c[[2]]]},
           _, {-1, -1, -1}
          ],
         22,
         A = Outer[Intersection, z1c, z2c, 1];
         Which[
             MatchQ[A, {{{__}, {_}}, {{}, {__}}}], {1, id, Join[z1c, {A[[1, 2]]}]},
             MatchQ[A, {{{_}, {__}}, {{__}, {}}}], {1, id, Join[z1c, {A[[1, 1]]}]},
             MatchQ[A, {{{__}, {}}, {{_}, {__}}}], {1, id, Join[Reverse[z1c], {A[[2, 1]]}]},
             MatchQ[A, {{{}, {__}}, {{__}, {_}}}], {1, id, Join[Reverse[z1c], {A[[2, 2]]}]},
             True, {-1, -1, -1}
          ],
         _,
         {-1, -1, -1}
         ]
    ]

    
IBDTransitionType2[{z1_, z2_}] :=
    Module[ {z1z2, z1c, z2c, id, len, min, moves, ss, pos, A, lena, sgna, 
      pos2, k, oper, i, pos12, pos21, sslen, pos1x, posx1},
        z1z2 = Intersection[z1, z2];
        z1c = Complement[z1, z1z2];
        z2c = Complement[z2, z1z2];
        id = Length[z1c] 10 + Length[z2c];
        If[ MemberQ[{22, 23, 32, 33, 24, 42, 34, 43, 44}, id],
            A = Outer[Intersection, z1c, z2c, 1];
            lena = Map[Length, A, {2}];
            sgna = Sign[lena];
        ];
        Switch[id,
         0, {0, 0, 0},
         12,
         len = Length[#] & /@ z2c;
         min = Min[len];
         moves = Flatten[Pick[z2c, len, min]];
         Switch[min,
          1, {min, 12, Join[z1c, {{}, moves}]},
          2, {min, 12, Join[z1c, {{}, {#}}] & /@ moves},
          _, {-1, -1, -1}
          ],
         21,
         len = Length[#] & /@ z1c;
         min = Min[len];
         moves = Flatten[Pick[z1c, len, min]];
         Switch[min,
          1, {min, 21, Append[z1c[[Ordering[len]]], moves]},
          2, {min, 21, 
           Append[If[ MemberQ[First[z1c], #],
                      z1c,
                      Reverse[z1c]
                  ], {#}] & /@ moves},
          _, {-1, -1, -1}
          ],
         22,
         ss = Outer[List, Total[sgna, {2}], Total[sgna]] sgna;
         pos = Position[ss, {2, 2}];
         If[ Length[pos] == 1,
             moves = Flatten[Extract[A, pos]];
             oper = If[ Length[moves] > 2,
                        {},
                        Append[If[ MemberQ[First[z1c], #],
                                   z1c,
                                   Reverse[z1c]
                               ], {#}] & /@moves
                    ];
             If[ Length[moves] > 1,
                 Switch[First[pos],
                  {1, 1},
                  If[ Extract[lena, {{1, 2}, {2, 1}}] == {1, 1},
                      oper = 
                       Join[oper, {Append[z1c, A[[1, 2]]], {z1c[[1]], {}, A[[1, 2]]},
                          Append[Reverse[z1c], A[[2, 1]]]}]
                  ],
                  {1, 2},
                  If[ Extract[lena, {{1, 1}, {2, 2}}] == {1, 1},
                      oper = 
                       Join[oper, {Append[z1c, A[[1, 1]]], {z1c[[1]], {}, A[[1, 1]]},
                          Append[Reverse[z1c], A[[2, 2]]]}]
                  ],
                  {2, 1},
                  If[ Extract[lena, {{1, 1}, {2, 2}}] == {1, 1},
                      oper = 
                       Join[oper, {Append[z1c, A[[1, 1]]], {z1c[[2]], {}, A[[2, 2]]},
                          Append[Reverse[z1c], A[[2, 2]]]}]
                  ],
                  {2, 2},
                  If[ Extract[lena, {{1, 2}, {2, 1}}] == {1, 1},
                      oper = 
                       Join[oper, {Append[z1c, A[[1, 2]]], {z1c[[2]], {}, A[[2, 1]]},
                          Append[Reverse[z1c], A[[2, 1]]]}]
                  ]
                  ]
             ];
             If[ oper === {},
                 {-1, -1, -1},
                 k = Min[2, Length[moves]];
                 {k, 22, 
                 If[ k == 1,
                     First[oper],
                     oper
                 ]}
             ],
             oper =
              If[ lena[[1, 2]] + lena[[2, 1]] === 2,
                  {Append[z1c, A[[1, 2]]], Append[Reverse[z1c], A[[2, 1]]]},
                  {}
              ];
             If[ lena[[1, 1]] + lena[[2, 2]] === 2,
                 oper = 
                  Join[oper, {Append[z1c, A[[1, 1]]], 
                    Append[Reverse[z1c], A[[2, 2]]]}]
             ];
             If[ oper === {},
                 {-1, -1, -1},
                 {2, 22, oper}
             ]
         ],
         13,
         len = Length[#] & /@ z2c;
         min = Min[len];
         moves = Flatten[Pick[z2c, len, min]];
         If[ Count[len, 1] >= 2,
             {2, 13, 
             Join[z1c, {{}, {#}}] & /@ moves},
             {-1, -1, -1}
         ],
         31,
         len = Length[#] & /@ z1c;
         moves = Pick[z1c, len, Min[len]];
         Switch[Count[len, 1],
          2, {2, 31, Join[{#}, Complement[z1c, moves], {#}] & /@ moves},
          3, {2, 31, Append[#, Flatten[#]] & /@ Subsets[moves, {2}]},
          _, {-1, -1, -1}
          ],
         23,
         Switch[Total[sgna, {2}],
          {1, 3} | {3, 1},
          i = If[ Total[sgna[[1]]] == 3,
                  1,
                  2
              ];
          k = Position[sgna[[3 - i]], 1][[1, 1]];
          pos = Delete[Range[3], k];
          moves = A[[i, k]];
          If[ Length[moves] == 1,
              oper = {z1c[[i]], {}, #} & /@ 
                Pick[A[[i, pos]], lena[[i, pos]], 1];
              AppendTo[oper, Append[If[ i == 1,
                                        z1c,
                                        Reverse[z1c]
                                    ], moves]];
              {2, 23, oper},
              {-1, -1, -1}
          ],
          {2, 2},
          k = Position[Total[sgna] , 2][[1, 1]];
          pos = Delete[Range[3], k];
          Switch[Total[lena[[All, k]]],
           2,
           oper = {{z1c[[1]], {}, 
              Flatten[Pick[A[[1]], lena[[1]], 1]]}, {z1c[[2]], {}, 
              Flatten[Pick[A[[2]], lena[[2]], 1]]}};
           If[ Max[lena[[1, pos]]] == 1,
               AppendTo[oper, Append[Reverse[z1c], A[[2, k]]]]
           ];
           If[ Max[lena[[2, pos]]] == 1,
               AppendTo[oper, Append[z1c, A[[1, k]]]]
           ];
           {2, 23, oper},
           _?(# > 2 &),
           i = If[ lena[[1, k]] > 1,
                   1,
                   2
               ];
           If[ Total[lena[[i, pos]]] == 1,
               {2, 
                23, {{z1c[[i]], {}, Flatten[A[[i, pos]]]}, 
                 Append[If[ i == 1,
                            Reverse[z1c],
                            z1c
                        ], A[[3 - i, k]]]}},
               {-1, -1, -1}
           ],
           _, {-1, -1, -1}
           ],
          _, {-1, -1, -1}
          ],
         32,
         Switch[Total[sgna],
          {1, 3} | {3, 1},
          i = If[ Total[sgna[[All, 1]]] == 3,
                  1,
                  2
              ];
          k = Position[sgna[[All, 3 - i]], 1][[1, 1]];
          pos = Delete[Range[3], k];
          If[ Length[A[[k, i]]] == 1 && Count[lena[[pos, i]], 1] >= 1,
              If[ lena[[pos, i]] == {1, 1},
                  {2, 32, 
                   Prepend[{z1c[[k]], #, A[[k, i]]} & /@ A[[pos, i]], 
                    Append[A[[pos, i]], Flatten[A[[pos, i]]]]]},
                  If[ ! OrderedQ[lena[[pos, i]]],
                      pos = Reverse[pos]
                  ];
                  {2, 
                   32, {Append[A[[pos, i]], A[[pos[[1]], i]]], {z1c[[k]], 
                     z1c[[pos[[2]]]], A[[k, i]]}}}
              ],
              {-1, -1, -1}
          ],
          {2, 2},
          k = Position[Total[sgna, {2}] , 2][[1, 1]];
          pos = Delete[Range[3], k];
          pos2 = DeleteCases[Position[sgna, 1], _?(First[#] == k &)];
          oper = {};
          Scan[If[ lena[[k, #[[2]]]] == 1,
                   AppendTo[oper, {z1c[[k]], Extract[A, #], A[[k, #[[2]]]]}]
               ] &, 
           pos2, {1}];
          Scan[If[ Extract[lena, #] == 1 && lena[[k, 3 - #[[2]]]] == 1,
                   AppendTo[oper, {Extract[A, #], z1c[[k]], Extract[A, #]}]
               ] &, 
           pos2, {1}];
          {2, 32, oper},
          _, {-1, -1, -1}
          ],
         33,
         ss = Outer[List, Total[sgna, {2}], Total[sgna]] sgna;
         pos12 = Position[ss, {1, 2}];
         pos21 = Position[ss, {2, 1}];
         Switch[{Length[pos12], Length[pos21]},
          {2, 2},
          sslen = Outer[List, Total[lena, {2}], Total[lena]] sgna;
          pos1x = Position[sslen, {1, _}];
          posx1 = Position[sslen, {_, 1}];
          If[ pos1x =!= {} && posx1 =!= {},
              {2, 33, 
               Join[{Append[Extract[A, pos12], Flatten[Extract[A, pos1x]]],
                 {z1c[[pos21[[1, 1]]]], {}, Flatten[Extract[A, posx1]]}},
                Flatten[
                 Outer[{z1c[[pos21[[1, 1]]]], {#1}, {#2}} &, 
                  Flatten[Extract[A, pos1x]], Flatten[Extract[A, posx1]]], 1]]},
              {-1, -1, -1}
          ],
          {2, 0},
          k = First[Complement[Range[3], pos12[[All, 1]]]];
          {2, 33, {z1c[[k]], Extract[A, #], A[[k, #[[2]]]]} & /@ pos12},
          {0, 2},
          pos = Sort[pos21[[All, 1]]];
          i = First[Complement[Range[3], pos]];
          k = First[Complement[Range[3], pos21[[All, 2]]]];
          {2, 33, {z1c[[#]], z1c[[i]], A[[#, k]]} & /@ pos},
          {1, 1},
          i = First[Complement[Range[3], {pos12[[1, 1]], pos21[[1, 1]]}]];
          k = First[Complement[Range[3], {pos12[[1, 2]], pos21[[1, 2]]}]];
          {2, 33, {{z1c[[pos21[[1, 1]]]], z1c[[i]], 
             A[[pos21[[1, 1]], k]]}, {z1c[[i]], z1c[[pos12[[1, 1]]]], 
             A[[i, pos12[[1, 2]]]]}}},
          _, {-1, -1, -1}
          ],
         24,
         moves = MapThread[Flatten[Pick[#1, #2, 1]] &, {A, lena}];
         If[ Positive[Length[#] & /@ moves] === {True, True},
             {2, 24, MapThread[{#1, {}, #2} &, {z1c, moves}]},
             {-1, -1, -1}
         ],
         42,
         moves = 
          MapThread[
           Flatten[Pick[#1, #2, 1]] &, {Transpose[A], Transpose[lena]}];
         If[ Positive[Length[#] & /@ moves] === {True, True},
             {2, 42, 
              MapThread[
               Append[#1, #2] &, {MapThread[
                 Pick[#1, #2, 1] &, {Transpose[A], Transpose[sgna]}], moves}]},
             {-1, -1, -1}
         ],
         34,
         ss = Outer[List, Total[sgna, {2}], Total[sgna]] sgna;
         pos = Position[ss, {2, 2}];
         pos2 = Position[ss, {1, 2}];
         If[ Extract[lena, pos] =!= {1} || Length[pos2] =!= 1,
             Return[{-1, -1, -1}]
         ];
         k = First[Complement[Range[3], pos[[All, 1]], pos2[[All, 1]]]];
         moves = Flatten[Pick[A[[k]], lena[[k]], 1]];
         If[ moves =!= {},
             {2, 34, {{z1c[[k]], {}, moves}, {z1c[[pos[[1, 1]]]], 
                z1c[[pos2[[1, 1]]]], Extract[A, First[pos]]}}},
             {-1, -1, -1}
         ],
         43,
         ss = Outer[List, Total[sgna, {2}], Total[sgna]] sgna;
         pos = Position[ss, {2, 2}];
         pos2 = Flatten[Position[Total[sgna], 2]];
         If[ Extract[lena, pos] =!= {1} || Length[pos2] =!= 2,
             Return[{-1, -1, -1}]
         ];
         k = First[Complement[pos2, pos[[All, 2]]]];
         moves = Flatten[Pick[A[[All, k]], lena[[All, k]], 1]];
         If[ moves =!= {},
             {2, 43, {Append[Pick[A[[All, k]], sgna[[All, k]], 1], moves],
               {z1c[[pos[[1, 1]]]], 
                Complement[z2c[[pos[[1, 2]]]], Extract[A, pos[[1]]]], 
                Extract[A, pos[[1]]]}}
              },
             {-1, -1, -1}
         ],
         44,
         ss = Outer[List, Total[sgna, {2}], Total[sgna]] sgna;
         pos = Position[ss, {2, 2}];
         If[ Extract[lena, pos] === {1, 1},
             {2, 44, {z1c[[First[#]]], 
                 Complement[z2c[[Last[#]]], Extract[A, #]], Extract[A, #]} & /@
                pos},
             {-1, -1, -1}
         ],
         _, {-1, -1, -1}
         ]
    ]
   
(*IBDTransitionType[{z1_, z2_}] :=
    Module[ {z1z2, z1c, z2c, len, share, i, j},
        z1z2 = Intersection[z1, z2];
        z1c = Complement[z1, z1z2];
        z2c = Complement[z2, z1z2];
        Switch[{Length[z1c], Length[z2c]},
         {0, 0}, {0, 0},
         {2, 1},
         len = Length[#] & /@ z1c;
         Switch[len,  
          {1, 1}, {21, Append[z1c, Flatten[z1c]]},
          {1, _}, {21, Append[z1c, z1c[[1]]]},
          {_, 1}, {21, Append[Reverse[z1c], z1c[[2]]]},
           _, {-1, -1}
          ],
         {1, 2},
         len = Length[#] & /@ z2c;
         If[ len === {} || Min[len] > 1,
             {-1, -1},
             {12, Join[z1c, {{},Flatten[Pick[z2c, len, 1]]}]}
         ],
         {2, 2},
         share = Table[Intersection[z1c[[i]], z2c[[j]]], {i, 2}, {j, 2}];
         Which[
          MatchQ[share, {{{__}, {_}}, {{}, {__}}}], {22, Append[z1c, {share[[1, 2, 1]]}]}, 
          MatchQ[share, {{{_}, {__}}, {{__}, {}}}], {22, Append[z1c, {share[[1, 1, 1]]}]},
          MatchQ[share, {{{__}, {}}, {{_}, {__}}}], {22, Append[Reverse[z1c], {share[[2, 1, 1]]}]}, 
          MatchQ[share, {{{}, {__}}, {{__}, {_}}}], {22, Append[Reverse[z1c], {share[[2, 2, 1]]}]},
          True, {-1, -1}
          ],
         _, {-1, -1}
         ]
    ]
*)           
    
IBDTransitionProb[zlist : {_, __}, n_, theta_] :=
    Module[ {pairs, typetran, types, pos, res, index, z1list, len},
        pairs = Partition[zlist, 2, 1];
        typetran = Rest[IBDTransitionType1[#]] & /@ pairs;
        types = typetran[[All, 1]];
        pos = Flatten[Position[types, #]] & /@ {-1, 0, 12, 21, 22};
        res = ConstantArray[0, Length[pairs]];
        (*0*)
        index = pos[[2]];
        If[ index =!= {},
            z1list = pairs[[index, 1]];
            len = Map[Length, z1list, {2}];
            res[[index]] = (theta ((Count[#, 1] & /@ len) + 1) + Total[len^2, {2}] + n)/((n + theta) (n + 1))
        ];
        (*12*)
        index = pos[[3]];
        If[ index =!= {},
            res[[index]] = (theta (Length[#] & /@ typetran[[index, 2, 3]]))/((n + theta) (n + 1))
        ];
        (*21*)
        index = pos[[4]];
        If[ index =!= {},
            len = Map[Length, typetran[[index, 2, 2 ;;]], {2}];
            len = Times @@ # & /@ len;
            res[[index]] = len/((n + theta) (n + 1))
        ];
        (*22*)
        index = pos[[5]];
        If[ index =!= {},
            len = Length[#] & /@ typetran[[index, 2, 2]];
            res[[index]] = len/((n + theta) (n + 1))
        ];
        N[res]
    ]    


End[] (* End Private Context *)

EndPackage[]