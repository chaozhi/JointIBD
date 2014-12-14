(* Mathematica Package *)

(* Created by the Wolfram Workbench Jan 23, 2011 *)

BeginPackage["Metropolis`"]
(* Exported symbols added here with SymbolName::usage *) 

DomainConstraint::usage = "DomainConstraint is an option to specify a constraint which random variable satisfies."

PrintDetails::usage = "PrintDetails is an option to specify whether print all the details for debugging."

StartLogpdf::usage = "StartLogpdf is an option to provide Log[p[x0]] where p[x] is the target distribution and x0 is the start value.";

IsReturnLast::usage = "IsReturnLast is an option to speficify whether return the last Metropolis update."

Metropolis::usage = "Metropolis[logpdf,logpdfjump,randomjump,x0,n]generates n Metropolis updates starting from x0. The target distribution is p[x] with the proposal distribution  J[xt, x] randomized by x = randomjump[xt] given current xt. Here logpdf[x] = Log[p[x]] and logpdfjump[xt,x] = Log[J[xt,x]]."

SamplingPattern::usage = "SamplingPattern={DecomposeMethod,StudentTDegree,IsSamplingComponentwise}. DecomposeMethod (\"Cholesky\" or \"Eigen\") specifies the method of decomposing the covariance matrix when sampling from multivariate student T with degree=StudentTDegree (>=3). IsSamplingComponentwise specifies whether sampling  is componentwise.";

IsDiscrete::usage = "IsDiscrete is an option to specify whether the random variable takes integers."

StudentTDisplacements::usage = "StudentTDisplacements][s,V,n] gives n random samples from student T distribution with mean 0, covaraince s V, and the default Infinity degree."

Metropolis::usage = "Metropolis[logpdf,logpdfjump,randomjump,x0,n]generates n Metropolis updates starting from x0. The target distribution is p[x] with the proposal distribution  J[xt, x] randomized by x = randomjump[xt] given current xt. Here logpdf[x] = Log[p[x]] and logpdfjump[xt,x] = Log[J[xt,x]]."

BlockMetropolis::usage = "BlockMetropolis[logpdf, logpdfjump, randomjump, x0, n] generates n updates for a set of conditional independent random variables starting from x0. The target distribution is p[x] with the proposal distribution  J[xt, x] randomized by x = randomjump[xt] given current xt. For any subset x[[indices]], logpdf[x, indices] = Log[p[x[[indices]]]] and logpdfjump[xt,x, indices] = Log[J[xt[[indices]],x[[indices]]]]."

StudentTMetropolis::usage = "StudentTMetropolis[logpdf, x0, s, V, n] generates n Metropolis udpates for the target distribution p[x]. logpdf[x] = Log[p[x]], and the symmetric proprosal distribution is the student T distribution with mean being the current value, (co)variance s^2 V, and degree infinity."

BlockStudentTMetropolis::usage = "BlockStudentTMetropolis[logpdf, x0, s, V, n] generates n updates for a set of conditional indepdent random variables starting from x0. The target distribution p[x] is provided by logpdf[x, indices] = Log[p[x, indices]] for any subset x[[indices]] of x. The symmetric proprosal distribution is the student T distribution with mean being the current value, (co)variance s^2 V, and degree infinity."

AdaptiveMetropolis::usage = "AdaptiveMetropolis[logpdf, x0, size,stepsize, adpfraction,adprate,adpset] generates size udpates starting from x0. Each update is given by the last update of  stepsize StudentTMetropolis updatings. The target distribution p[x] is provided by logpdf[x] = Log[p[x]], and the symmetric proprosal distribution is the student T distribution. The variance of the proposal distribution is adapted in the initial Floor[size adpfraction] iterations with rate adprate, and the adaptive parameter set adpset."

BlockAdaptiveMetropolis::usage = "BlockAdaptiveMetropolis[logpdf, x0, size, stepsize, adpfraction, adprate, adpset] generates size updates for a list of conditional indepdent random variables starting from x0. Each update is given by the last update of  stepsize BlockStudentTMetropolis updatings. The target distribution p[x] is provided by logpdf[x, indices] = Log[p[x, indices]] for any subset x[[indices]] of x. The symmetric proprosal distribution is the student T distribution. The variance of the proposal distribution is adapted in the initial Floor[size adpfraction] iterations with rate adprate, and the adaptive parameter set adpset."

MultiTryMetropolis::usage = "MultiTryMetropolis[logpdf,logpdfjump,randommultijump, x0,n,isindep] generates one Metropolis update starting from x0.  The target distribution is p[x] with the proposal distribution J[xt,x] given current value xt. The  randommultijump[xt, n] generates n proposal tries from J[xt,x].  Here logpdf[x] = Log[p[x]] and logpdfjump[xt,x] = Log[J[xt,x]].  The isindep specifies whether the proposal distribution is independent of current values."

InitialProposalScale::usage = "InitialProposalScale[targetratio, r,q] initializes scaleadpset with target accept ratio targetratio, and the default r=0.5 and q=0.005."

AdaptiveProposalScale::usage ="AdaptiveProposalScale[isadp, indaccept, adrate, scaleadpset] gives updated scaleadpset for one random variable or a set of (conditional) independent random variables. It calculates the recent accept rates, and the scales of proposal distributions with adpative rate adprate if isadp."

InitialProposal::usage = "InitialProposal[targetratio,x0,r, q] initializes adpative parameter set adpset  from the target accept ratio targetratio and initial value x0. The optional  r=0.5 and q=0.005. \nInitialProposal[targetratiols,xls0,r, q] initializes adpative parameter set adpset  from the target accept ratio targetratiols and initial values xls0. The optional r=0.5 and  q=0.005." 

AdaptiveProposal::usage = "AdaptiveProposal[isadp, indaccept, adprate, adpset, xls] gives updated adpset. It calculates the recent accept rate, and the scale and the covariance of proposal distributions with adpative rate adprate if isadp. The xls is a list of random sample x, and the indaccept is the average accept probability of the corresponding prposal list. \nAdaptiveProposal[isadp, indacceptls, adprate, adpset, {xls_ 1,..}] gives updated adpset for a set of random variables. It calculates the recent list of accept rates, and the scale and the covariance of proposal distributions with adpative rate adprate if isadp. The xls_i is a list of random sample x_i for i=1,.., and the indaccept is the corresponding list of average accept probabilities."

Begin["`Private`"]

(* Implementation of the package *)

(**********************************************Metropolis******************************************)

Options[Metropolis] = {
    DomainConstraint -> (True &), 
    StartLogpdf->None,
    IsReturnLast -> True,
    PrintDetails-> False
    };    
 
Metropolis::badratio :=
    "{oldlike,proplike,oldtoprop,proptoold}={`1`,`2`,`3`,`4`} are not real values, given the proposal from `5` to `6`!";
Metropolis[logpdf_Function, logpdfjump_Function, 
   randomjump_Function, start_, size_Integer?Positive, 
   opts : OptionsPattern[]] :=
    Module[ {constraint,startlogpdf,isreturnlast, printdetails,oneupdate, startpack, res},
        {constraint,startlogpdf, isreturnlast,printdetails} = 
         OptionValue@{DomainConstraint,StartLogpdf, IsReturnLast, PrintDetails};
        {isreturnlast,printdetails} = TrueQ/@{isreturnlast,printdetails};
        If[ startlogpdf===None,
            startlogpdf = logpdf[start]
        ];
        (*{x, logpdf[x],1/0 (accept or not)}*)
        oneupdate =
            Function[{x},
                Module[ {old, oldlike, prop, proplike, oldtoprop, proptoold,logratio, probaccept},
                    {old, oldlike} = x[[;;2]];
                    prop = randomjump[old];
                    Which[
                         prop === old, {old, oldlike, 1,1},
                         ! constraint[prop], {old, oldlike, 0,0},
                         True,
                         proplike = logpdf[prop];
                         oldtoprop = logpdfjump[old, prop];
                         proptoold = logpdfjump[prop, old];
                         logratio = (proplike + proptoold) - (oldlike +oldtoprop );
                         If[ printdetails,
                             Print["proprosal =",prop,"\n logpdf[old]=",oldlike,"\n logpdf[prop]=",proplike,"\n logpdfjump[old,prop]=",oldtoprop,
                                 "\n logpdfjump[prop,old]=",proptoold, "\n Log[MHratio]=",logratio];
                         ];
                         If[ !(Head[N[logratio]]===Real||logratio===-Infinity),
                             Message[Metropolis::badratio, oldlike,proplike, oldtoprop, proptoold, old, prop];
                             Return[$Failed]
                         ];
                         probaccept = Min[1,Exp[logratio]];
                         If[ RandomReal[] < probaccept,
                             {prop, proplike,1,probaccept},
                             {old, oldlike, 0,probaccept}
                         ]
                     ]
                ]
            ];
        startpack = {start, startlogpdf,0,0};
        res = Rest[NestList[oneupdate, startpack, size]];
        If[ isreturnlast,            
            res[[-1,-2;;]] = Mean[res[[All,-2;;]]];
            res = Last[res];
        ];
        res
    ]
    

(********************************************BlockMetropolis***************************************)

Options[BlockMetropolis] = Options[Metropolis];

BlockMetropolis::badratio :=
    "{oldlike,proplike, oldtoprop, proptoold}={`1`,`2`,`3`,`4`} are not real numbers, given the proposal from `5` to `6`!";
BlockMetropolis[logpdf_Function, logpdfjump_Function, 
   randomjump_Function, start_, size_Integer?Positive, 
   opts : OptionsPattern[]] :=
    Module[ {blockOneStepUpdate,constraint,startlogpdf,isreturnlast,printdetails, startpack, res,dim},
        {constraint,startlogpdf, isreturnlast,printdetails} = 
         OptionValue@{DomainConstraint,StartLogpdf, IsReturnLast,PrintDetails};
        {isreturnlast,printdetails} = TrueQ/@{isreturnlast,printdetails};
        dim = Length[start];
        If[ startlogpdf===None,
            startlogpdf = logpdf[start,Range[dim]]
        ];
        (*{x, logpdf[x],1/0 (accept or not),Min[1,MHratio]}*)
        blockOneStepUpdate =
            Function[x,
                Module[ {old, oldlike, prop, proplike, logratio, oldtoprop, proptoold,indaccept,indices,pos,probaccept},
                    {old, oldlike} = x[[;;2]];
                    prop = randomjump[First[x]];
                    If[ prop === old,
                        {old, oldlike, ConstantArray[1,dim],ConstantArray[1,dim]},
                        proplike = indaccept = probaccept = ConstantArray[0,dim];
                        indices = If[ constraint===(True&),
                                      Range[dim],
                                      Flatten[Position[constraint[prop], True]]
                                  ];
                        If[ indices==={},
                            {old,oldlike,indaccept,indaccept},
                            proplike[[indices]] = logpdf[prop,indices];
                            oldtoprop = logpdfjump[old, prop,indices];
                            proptoold = logpdfjump[prop, old,indices];
                            logratio = (proplike[[indices]] + proptoold) - (oldlike[[indices]] +oldtoprop);
                            If[ printdetails,
                                Print["proprosal =",prop,"\n logpdf[prop]=",proplike,"\n logpdfjump[old,prop]=",oldtoprop,
                                    "\n logpdfjump[prop,old]=",proptoold, "\n Log[MHratio]=",logratio];
                            ];
                            If[ !VectorQ[N[logratio],Head[#] === Real||#===-Infinity &],
                                Message[BlockMetropolis::badratio, oldlike,proplike, oldtoprop, proptoold, old, prop];
                                Return[$Failed]
                            ];
                            indaccept[[indices]] = Boole[Thread[Log[RandomReal[{0, 1}, Length[logratio]]] < logratio]];
                            pos = Flatten[Position[indaccept,1,{1},Heads->False]];
                            old[[pos]] = prop[[pos]];
                            oldlike[[pos]] = proplike[[pos]];
                            probaccept[[indices]] = Exp[logratio/.{_?Positive -> 0}];
                            {old,  oldlike, indaccept,probaccept}
                        ]
                    ]
                ]
            ];
        startpack = {start, startlogpdf, ConstantArray[0,dim],ConstantArray[0,dim]};
        res = Rest[NestList[blockOneStepUpdate, startpack, size]];
        If[ isreturnlast,            
            res[[-1,-2;;]] = Mean[res[[All,-2;;]]];
            res = Last[res];
        ];
        res
    ]
    
(**********************************************StudentTDisplacements*******************************)

Options[StudentTDisplacements] = {
    SamplingPattern->{"Eigen",Infinity,True},
    IsDiscrete->False    
}

StudentTDisplacements::badpatt = "The option SamplingPattern must be a list of length 3."
StudentTDisplacements::baddf = "StudentTDegree of  SamplingPattern->{DecomposeMethod, StudentTDegree, IsSamplingComponentwise} must be either Infinity or a positive integer no less than 3."
StudentTDisplacements::badmethod = "DecomposeMethod of  SamplingPattern->{DecomposeMethod, StudentTDegree, IsSamplingComponentwise} must be either \"Cholesky\" or \"Eigen\"."
StudentTDisplacements[scale_,cov_,size_,OptionsPattern[]] :=
    Module[ {patt,decomposemethod,deg,issc,isdiscrete,dim,dist,delta,V},
        {patt,isdiscrete} =
            OptionValue@{SamplingPattern,IsDiscrete};
        If[ Length[patt]==3,
            {decomposemethod,deg,issc} = patt;
            If[ !((deg>=3&&IntegerQ[deg])||deg===Infinity),
                Message[StudentTDisplacements::baddf];
                Return[$Failed]
            ];
            If[ !MemberQ[{"Cholesky","Eigen"},decomposemethod],
                Message[StudentTDisplacements::badmethod];
                Return[$Failed]
            ],
            Message[StudentTDisplacements::badpatt];
            Return[$Failed]
        ];
        {issc,isdiscrete} = TrueQ/@{issc,isdiscrete};
        dim = Length[cov];
        If[ dim===1,
            dist = If[ deg===Infinity,
                       NormalDistribution[],
                       StudentTDistribution[deg]
                   ];
            delta = scale RandomReal[dist,size] Sqrt[cov[[1,1]]],
            V = cov+10^(-6.) IdentityMatrix[dim];
            If[ deg!=Infinity,
                V = V (deg-2)/deg
            ];
            delta = Switch[decomposemethod,
                "Cholesky",randomMultiNormCholesky[scale,V,size,issc],
                "Eigen", randomMultiNormEigen[scale,V,size,issc]                
            ];
            If[ deg!=Infinity,
                delta = delta Sqrt[deg/RandomReal[ChiSquareDistribution[deg],Length[delta]]]
            ];
        ];
        If[ isdiscrete,
            Round[delta],
            delta
        ]
    ];

randomMultiNormCholesky[scale_,V_,n_,issc_] :=
    Module[ {up,z,dim,temp},
        dim = Length[V];
        up = CholeskyDecomposition[(V+Transpose[V])/2];
        If[ issc,
            z = RandomReal[NormalDistribution[],{Ceiling[n/dim],dim}];
            temp = Sqrt[dim] up scale # & /@ z;
            If[ Mod[n, dim] != 0,
                temp[[-1]] = RandomSample[temp[[-1]], Mod[n, dim]]
            ];
            Flatten[temp, 1],
            z = RandomReal[NormalDistribution[], {n, dim}];
            scale Transpose[up].#&/@z
        ]
    ]

randomMultiNormEigen[scale_,V_,n_,issc_] :=
    Module[ {eigenval,eigenvec,dim,z,temp},
        dim = Length[V];
        {eigenval, eigenvec} = Eigensystem[V];
        eigenval = Abs[eigenval];
        If[ issc,
            z = RandomReal[NormalDistribution[], {Ceiling[n/dim], dim}];
            temp = Re[Sqrt[eigenval] Sqrt[dim]  scale eigenvec];
            temp = (temp # & /@ z);
            If[ Mod[n, dim] != 0,
                temp[[-1]] = RandomSample[temp[[-1]], Mod[n, dim]]
            ];
            Flatten[temp, 1],
            z = RandomReal[NormalDistribution[], {n, dim}];
            scale Re[Transpose[eigenvec].DiagonalMatrix[Sqrt[eigenval]].# & /@z]
        ]
    ]
    

(**********************************************StudentTMetroplis***********************************)

Options[StudentTMetropolis] = Join[
    Options[Metropolis],
    Options[StudentTDisplacements]        
    ];

StudentTMetropolis::baddim = "Dimensions for arguments and/or options do not match!";
StudentTMetropolis[
    logpdf_Function,
    start_?(NumericQ[#]||VectorQ[#,NumericQ]&),
    scale_?(NumericQ[#]||VectorQ[#,NumericQ]&),
    cov_?(MatrixQ[#,NumericQ]&&Equal@@Dimensions[#]&),
    size_Integer?Positive,
    opts : OptionsPattern[]] :=
    Module[ {displacements,randomprop,n},
        displacements = StudentTDisplacements[scale,cov,size,FilterRules[{opts}, Options[StudentTDisplacements]]];
        n = Length[displacements];
        $index = 0;
        randomprop = Function[{xt},
                        $index = Mod[$index,n]+1;
                        xt+displacements[[$index]]
                    ];
        Metropolis[logpdf,(0&),randomprop,start,n,FilterRules[{opts}, Options[Metropolis]]]
    ]/;(Length[scale]==0||Length[scale]==Length[cov])  

(********************************************BlockStudentTMetroplis***********************************)

Options[BlockStudentTMetropolis] = Join[
    Options[BlockMetropolis],
    Options[StudentTDisplacements]        
    ]; 

BlockStudentTMetropolis::baddim = "Dimensions for arguments and/or options do not match!"; 
BlockStudentTMetropolis::isdiscrete = "The option IsDiscrete must be either boolean (True|False) or a list of boolean values with length equal to the number of varibles!"
BlockStudentTMetropolis::sampatt = "The option SamplingPattern must be either {} or a list of them with length equal to the number of varibles!"
BlockStudentTMetropolis[logpdf_Function, 
  start : {_?(NumericQ[#] || VectorQ[#, NumericQ] &) ..}, 
  scale_?(VectorQ[#, NumericQ] &), 
  cov : {_?(Equal @@ Dimensions[#] && MatrixQ[#, NumericQ] &) ..}, 
  size_Integer?Positive, opts : OptionsPattern[]] :=
    Module[ {sampatt,isdiscrete,opts2,i, displacements, randomprop, n,dim},
        {sampatt,isdiscrete} = OptionValue@{SamplingPattern,IsDiscrete};
        dim = Length[start];
        Switch[Dimensions[sampatt],
            {3}, sampatt = Table[sampatt,{dim}],
            {dim,3}, 0,
            _, Message[BlockStudentTMetropolis::sampatt];
               Return[$Failed]
        ];
        Switch[Dimensions[isdiscrete],
            {},isdiscrete = Table[isdiscrete,{dim}],
            {dim},0,
            _,Message[BlockStudentTMetropolis::isdiscrete];
              Return[$Failed]
            ];
        opts2 = FilterRules[FilterRules[{opts}, Options[StudentTDisplacements]],Except[{SamplingPattern,IsDiscrete}]];
        displacements = Transpose[
          Table[StudentTDisplacements[scale[[i]], cov[[i]], size, SamplingPattern->sampatt[[i]], IsDiscrete->isdiscrete[[i]],opts2], {i, Length[scale]}]
          ];
        n = Length[displacements];
        $blockindex = 0;
        randomprop = Function[{xt},
          $blockindex = Mod[$blockindex, n] + 1;
          xt + displacements[[$blockindex]]
          ];
        BlockMetropolis[logpdf, (0 &), randomprop, start, n, FilterRules[{opts}, Options[BlockMetropolis]]]
    ] /; (Length[scale] == Length[cov]) 


    
(**************************************AdaptiveMetropolis*****************************************)
    
Options[AdaptiveMetropolis] = Options[StudentTMetropolis];

AdaptiveMetropolis[
   logpdf_Function,
   inputstart_ ?(NumericQ[#]||VectorQ[#, NumericQ]&),
   size_Integer?Positive,
   stepsize_Integer?Positive,
   adpfraction_?(NumericQ[#]&&0<=#<=1&),
   adprate_?NumericQ,   
   inputadpset_?(Length[#]==7&),
   opts : OptionsPattern[]] :=
    Module[ {dim,adpthreshold,start,startlogpdf, isreturnlast, adpset,temp, res, 
         opts2, scale, cov,i},
        {start,adpset} = {inputstart,inputadpset};
        {isreturnlast, startlogpdf} =  OptionValue@{IsReturnLast,StartLogpdf};
        opts2 = FilterRules[FilterRules[{opts}, Options[StudentTMetropolis]], 
                    Except[{IsReturnLast,StartLogpdf}]];
        adpthreshold = Floor[adpfraction size stepsize];
        dim = Length[start];
        res = Table[
                 {scale, cov} = adpset[[{-3, -1}]];
                 If[ dim<=1,
                     cov = {{1}}
                 ];
                 temp = StudentTMetropolis[logpdf, start, scale, cov, stepsize,IsReturnLast -> False,
                  StartLogpdf -> startlogpdf, opts2];
                 adpset = AdaptiveProposal[i<=adpthreshold, Mean[temp[[All,-1]]], adprate, adpset, temp[[All,1]]];                 
                 temp[[-1,-2;;]] = Mean[temp[[All,-2;;]]];
                 temp = Last[temp];
                 {start, startlogpdf} = temp[[;; 2]];
                 Append[temp,adpset],{i,1,size stepsize,stepsize}];                       
        (*list of {x, loglike, 1/0, ProbAccept, adpset}*)        
        If[ isreturnlast,            
            res[[-1,{3,4}]] = Mean[res[[All,{3,4}]]];
            res = Last[res];
        ];
        res
    ]
        
(**************************************BlockAdaptiveMetropolis************************************)
   
Options[BlockAdaptiveMetropolis] = Options[BlockStudentTMetropolis];

BlockAdaptiveMetropolis[
    logpdf_Function, 
       inputstart : {_?(NumericQ[#] || VectorQ[#, NumericQ] &) ..}, 
       size_Integer?Positive, 
       stepsize_Integer?Positive, 
       adpfraction_?(NumericQ[#] && 0 <= # <= 1 &),
       adprate_?NumericQ, 
       inputadpset_?(Length[#] ==7&), 
       opts : OptionsPattern[]] :=
    Module[ {dim, adpthreshold, start, startlogpdf, isreturnlast, adpset, temp, res, opts2, scale, cov,  i},
        {start, adpset} = {inputstart, inputadpset};
        {isreturnlast, startlogpdf} = OptionValue@{IsReturnLast, StartLogpdf};
        opts2 = 
         FilterRules[FilterRules[{opts}, Options[BlockStudentTMetropolis]],
           Except[{IsReturnLast, StartLogpdf}]];
        adpthreshold = Floor[adpfraction size stepsize];
        dim = Length[start];
        res = Table[
          {scale, cov} = adpset[[{-3, -1}]];
          cov = cov /. {{{_}} -> {{1}}};
          temp = BlockStudentTMetropolis[logpdf, start, scale, cov, stepsize, 
             IsReturnLast -> False, StartLogpdf -> startlogpdf, opts2];
          adpset = AdaptiveProposal[i<=adpthreshold, Mean[temp[[All,-1]]], adprate, adpset, temp[[All,1]]];          
          temp[[-1,-2;;]] = Mean[temp[[All,-2;;]]];
          temp = Last[temp];
          {start, startlogpdf} = temp[[;; 2]];
          Append[temp, adpset], {i, 1,size stepsize, stepsize}];
        (*list of {x,loglike,1/0,probaccept,adpset}*)
        If[ isreturnlast,            
            res[[-1,{3,4}]] = Mean[res[[All,{3,4}]]];
            res = Last[res];
        ];
        res
    ];
        
(*****************************************MultiTryMetropolis**************************************)

Options[MultiTryMetropolis] = {
    DomainConstraint -> (True &)
    };  
    
MultiTryMetropolis[logpdf_Function, logpdfjump_Function, randommultijump_Function, 
   start_,1,isindependent:True|False,opts : OptionsPattern[]] :=
    Module[ {randomjump},
        randomjump = First[randommultijump[#,1]]&;
        Metropolis[logpdf,logpdfjump,randomjump,start,1,FilterRules[{opts}, Options[Metropolis]]]
    ];  
    
MultiTryMetropolis[logpdf_Function, logpdfjump_Function, randommultijump_Function, 
   start_,ntry_Integer/;ntry>1, False,OptionsPattern[]] :=
    Module[ {constraint,weightfun,multitry, weighttry,baseweight,pos,prop,multiref,weightref,ratio},
        constraint = OptionValue[DomainConstraint];
        weightfun = Function[{x,y},(logpdf[x]+logpdfjump[x,y])];
        multitry = randommultijump[start,ntry];
        weighttry = Map[If[ constraint[#],
                            weightfun[#,start],
                            -Infinity
                        ]&,multitry];
        baseweight = Max[weighttry];
        If[ baseweight===-Infinity,
            Return[{start,0}]
        ];
        weighttry = Exp[weighttry-baseweight];
        pos = RandomChoice[weighttry -> Range[ntry]];
        prop = multitry[[pos]];
        multiref = Append[randommultijump[prop,ntry-1],start];
        weightref = Exp[(weightfun[#,prop]& /@ multiref)-baseweight];
        ratio = Total[weighttry]/Total[weightref];
        (*Print[{{ratio,Total[weighttry],Total[weightref]},pos,multitry,multiref}];*)
        If[ RandomReal[] < Min[1, ratio],
            {prop, 1},
            {start, 0}
        ]
    ];

MultiTryMetropolis[logpdf_Function,logpdfjump_,randommultijump_Function, start_, 
    ntry_Integer/; ntry>1, True,OptionsPattern[]] :=
    Module[ {constraint, weightfun,multitry, weighttry,baseweight,pos,prop,totweight,startweight,ratio},
        constraint = OptionValue[DomainConstraint];
        weightfun = Function[{x},(logpdf[x]-logpdfjump[0,x])];
        multitry = randommultijump[0,ntry];
        weighttry = Map[If[ constraint[#],
                            weightfun[#],
                            -Infinity
                        ]&,multitry];
        baseweight = Max[weighttry];
        If[ baseweight===-Infinity,
            {start,0},
            weighttry = Exp[weighttry-baseweight];
            pos = RandomChoice[weighttry -> Range[ntry]];
            prop = multitry[[pos]];
            totweight = Total[weighttry];
            startweight = Exp[weightfun[start]-baseweight];
            If[ weighttry[[pos]]>startweight,
                {prop,1},
                ratio = totweight/(totweight-weighttry[[pos]]+startweight);
                (*Print[{ratio,Total[weighttry],pos,multitry}];*)
                If[ RandomReal[] < ratio,
                    {prop, 1},
                    {start, 0}
                ]
            ]
        ]
    ];

(*************Adaptive Proposal for a scalar or a set of independent scalars************)

(*q: greedy of adaptation,which is exponentially decay*)
InitialProposalScale[
    targetratio_?(NumericQ[#]||VectorQ[#, NumericQ]&),  
    r_: 0.5,      
    q_: 0.05] :=
    Module[ {accept,scale},
        scale = Sign[targetratio];
        accept = targetratio;
        {targetratio,r,q,accept,scale}
    ] /; NumericQ[r]&&NumericQ[q]   
      
AdaptiveProposalScale[
   isadp:True|False,
   indaccept_ ?(NumericQ[#]||VectorQ[#, NumericQ]&),
   adprate_?NumericQ,
   adpset_?(Length[#] == 5&)] :=
    Module[ {targetratio,r,q,accept, scale},
        {targetratio,r,q, accept,scale} = adpset;
        accept+=q (indaccept-accept);
        If[ isadp,
            scale*= Exp[r(accept - targetratio)];
            r*=(1 - adprate);
            q*=(1 - adprate);
        ];
        {targetratio,r,q,accept,scale}
    ]

(***************************Adaptive Proposal  vectorized random variable(s)*********************)
InitialProposal[
    targetratio_?NumericQ,
    theta_ ?(NumericQ[#]||VectorQ[#, NumericQ]&),
    r_: 0.5,
    q_: 0.05] :=
    Module[ {accept,scale, mu, cov},
        scale = 1;
        accept = targetratio;
        mu = theta;
        cov = IdentityMatrix[Max[1,Length[theta]]];
        {targetratio,r,q,accept,scale, mu, cov}
    ] /; NumericQ[r]&&NumericQ[q]

InitialProposal[
  targetratio_?(VectorQ[#, NumericQ] &), 
  theta : {_?(NumericQ[#]||VectorQ[#, NumericQ] &) ..}, 
  r_: 0.5,
  q_: 0.05] :=
    Module[ {accept, scale, mu, cov},
        scale =  Sign[targetratio];
        accept = targetratio;
        mu = theta;
        cov = IdentityMatrix[Max[1,Length[#]]]&/@theta;
        {targetratio, r,q, accept, scale, mu, cov}
    ] /; NumericQ[r]&&NumericQ[q]&&(Length[targetratio]==Length[theta])
        
AdaptiveProposal[
   isadp:True|False,
   indaccept_ ?NumericQ,
   adprate_?NumericQ,
   adpset_?(Length[#] == 7&),
   thetals_?(ArrayQ[#,1|2, NumericQ]&)] :=
    Module[ {targetratio,r,q,accept, scale, mu, cov, temp},
        {targetratio,r,q, accept,scale, mu, cov} = adpset;
        accept+=q (indaccept-accept);
        If[ isadp,
            temp = Switch[Depth[thetals],                             
                    2,Total[(thetals - mu)^2]/Length[thetals],
                    3,Total[Transpose[{# - mu}].{# - mu} & /@ thetals]/Length[thetals]
            ];
            cov += q (temp - cov);
            mu += q (Mean[thetals] - mu);
            scale*=Exp[r (accept - targetratio)];
            r*=(1 - adprate);
            q*=(1 - adprate);
        ];
        {targetratio,r,q,accept,scale,mu,cov}
    ]
                 
AdaptiveProposal[
    isadp : True | False, 
    indaccept_?(VectorQ[#, NumericQ]&), 
    adprate_?NumericQ, 
    adpset_?(Length[#] == 7&), 
    thetals : {{_?(NumericQ[#] || VectorQ[#, NumericQ] &) ..} ..}] :=
    Module[ {targetratio,r,q, accept, scale, mu, cov, thetals2,temp},
        {targetratio, r,q, accept, scale, mu, cov} = adpset;
        accept+=q (indaccept-accept);
        If[ isadp,
            thetals2 = # - mu & /@ thetals;
            thetals2 = Replace[thetals2, x_?NumericQ :> {x}, {2}];
            temp = Map[Transpose[{#}].{#} &, thetals2, {2}];
            temp = Total[#] & /@ Transpose[temp];
            mu += q (Total[thetals]/Length[thetals]- mu);
            cov += q (temp - cov);
            scale*=Exp[r (accept - targetratio)];
            r*=(1 - adprate);
            q*=(1 - adprate);
        ];
        {targetratio, r,q, accept, scale, mu, cov}
    ]         

End[]

EndPackage[]

