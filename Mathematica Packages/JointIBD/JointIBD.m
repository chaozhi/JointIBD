(* Mathematica Package *)

(* Created by the Wolfram Workbench 2013-2-23 *)

BeginPackage["JointIBD`",{"JointIBD`Likelihood`","ChineseRestaurantProcess`","Metropolis`","Toolbox`"}]
(* Exported symbols added here with SymbolName::usage *) 

initializeJIBDchain::usage = "initializeJIBDchain  "

saveJIBDchain::usage = "saveJIBDchain  "

updateJIBDchain::usage = "updateJIBDchain  "

swapJIBDchain::usage = "swapJIBDchain  "

permuteJIBDchain::usage = "permuteJIBDchain  "

testibdls::usage = "testibdls  "

testsnpls::usage = "testsnpls  ";

(*getmidoperators::usage=""

randomcondibd::usage=""
*)

Begin["`Private`"]
(* Implementation of the package *)
        
initializeJIBDchain[inputsnpls_, nsq_, nbp_, initheta_,inirho_,inieps_,maxcpden_] :=
    Module[ {inithetashape, inithetascale, inithetamin, inithetamax, inirhoshape, inirhoscale, inirhomin, inirhomax, 
        iniepsalpha, iniepsbeta, iniepsmin, iniepsmax,snpls = inputsnpls,theta, rho, epsilon, ncp,ibd1, ibdls, ls2,
     map, ibds, adpeps, adptheta, adprho,indibd},
        {inithetashape, inithetascale, inithetamin, inithetamax} = initheta;
        {inirhoshape, inirhoscale, inirhomin, inirhomax} = inirho;
        {iniepsalpha, iniepsbeta, iniepsmin, iniepsmax} = inieps;
        epsilon = RandomReal[TruncatedDistribution[{iniepsmin, iniepsmax},BetaDistribution[iniepsalpha,iniepsbeta]]];
        theta = RandomReal[If[ inithetascale == Infinity,
                               UniformDistribution[{inithetamin, inithetamax}],
                               TruncatedDistribution[{inithetamin, inithetamax},GammaDistribution[inithetashape, inithetascale]]
                           ]];
        rho = RandomReal[If[ inirhoscale == Infinity,
                             UniformDistribution[{inirhomin, inirhomax}],
                             TruncatedDistribution[{inirhomin, inirhomax},GammaDistribution[inirhoshape, inirhoscale]]
                         ]];
        ncp = Min[Ceiling[rho nbp], maxcpden nbp/10^6];
        ibd1 = RandomChineseRestaurant[nsq, 0, theta,ReturnType->"Configuration"];
        ibdls = Transpose[{Prepend[Sort[RandomSample[Range[1, nbp], ncp]], 1],IBDFoldTransition[ibd1, nsq, theta, ncp]}];
        ibdls[[All, 2]] = Map[Sort, ibdls[[All, 2]], 2];
        If[ Length[ibdls]==1,
            ibdls[[1]] = Append[ibdls[[1]],LogPDFCRP[ibd1, nsq, 0, theta]],
            ibdls = Transpose[Append[Transpose[ibdls],Prepend[Log[IBDTransitionProb[ibdls[[All, 2]], nsq, theta]], 
                LogPDFCRP[ibd1, nsq, 0, theta]]]];
        ];
        ls2 = Append[ibdls, ReplacePart[ibdls[[-1]], 1 -> nbp + 1]];
        map = IndexByInterpolation[ls2[[All, 1]]];
        ibds = ls2[[map[snpls[[All, 1]]], 2]];
        snpls[[All, -1]] = IBDLogLikelihood[snpls[[All, 2]], ibds, snpls[[All, 3]],epsilon];
        adpeps = InitialProposalScale[0.44];
        adptheta = InitialProposalScale[0.44];
        adprho = InitialProposal[0.44,rho];
        indibd = Table[0, {10}];
        {{theta, epsilon, rho}, ibdls, {adpeps, adptheta,adprho,indibd,{}}, snpls}
    ]    

saveJIBDchain[mcstate_, adpheat_,outputconcise_, outputcold_] :=
    Module[ {mcstate2, mcstate3,i,ind},
    	ind=Flatten[Transpose[Table[Range[Length[mcstate]], {Dimensions[mcstate][[2]]}]]];
        mcstate2 = Flatten[mcstate,1][[All, ;; 5]];
        Do[
            mcstate2[[i, 5]] = Prepend[mcstate2[[i, 5]], adpheat[[ind[[i]]]]];
            mcstate2[[i, 4]] = {#[[1]], Length[#[[2]]]} &/@mcstate2[[i, 4]], {i, Length[mcstate2]}];
        mcstate2[[All, 5, ;; 4]] = mcstate2[[All, 5, ;; 4, {4, 5}]];
        PutAppend[mcstate2, outputconcise];
        mcstate3 = First[Select[#, First[#] == 1 &]] & /@ mcstate;
        Do[mcstate3[[i, 5]] = Prepend[mcstate3[[i, 5]],  adpheat[[ind[[i]]]]],{i,Length[mcstate3]}];
        mcstate3[[All, 5, ;; 4]] = mcstate3[[All, 5, ;; 4, {4, 5}]];
        PutAppend[mcstate3, outputcold];
    ]   
    
permuteJIBDchain[inputgroupmcstate_, isadp_, adprate_,inputadpset_,size_] :=
    Module[ {groupmcstate = inputgroupmcstate, adpset = inputadpset, ii, ar, ind = 0},
        Do[
            ii = RandomSample[Range[Length[groupmcstate]]];
            ar = groupmcstate[[All, 1]].groupmcstate[[ii, 2]] - 
              groupmcstate[[All, 1]].groupmcstate[[All, 2]];
            ar = Min[1, Exp[ar]];
            ind += ar;
            If[ RandomReal[] < ar,
                groupmcstate[[All, 1]] = groupmcstate[[ii, 1]];
            ],{size}];
        adpset = AdaptiveProposalScale[isadp, ind/size, adprate, adpset];
        adpset[[-1]] = Min[adpset[[-1]], 1./Length[groupmcstate]];
        groupmcstate = SortBy[groupmcstate, First];
        If[ isadp,
            groupmcstate[[All, 1]] = 1 - Last[adpset]  Range[Length[groupmcstate]-1,0, -1];
        ];
        {groupmcstate, adpset}
    ]       
    
swapJIBDchain[inputgroupmcstate_, isadp_, adprate_,inputadpset_,size_] :=
    Module[ {groupmcstate = inputgroupmcstate, adpset = inputadpset, ii,k, ar, ind},
        ind = Table[0,{size}];
        Do[
         ii = RandomSample[Range[Length[groupmcstate]], 2];
         ar = groupmcstate[[ii, 1]].groupmcstate[[Reverse[ii], 2]]-groupmcstate[[ii, 1]].groupmcstate[[ii, 2]];
         ar = Min[1, Exp[ar]];
         ind[[k]] = ar;
         If[ RandomReal[] < ar,
             groupmcstate[[ii, 1]] = groupmcstate[[Reverse[ii], 1]];
         ], {k,size}];
        groupmcstate = SortBy[groupmcstate, First];
        adpset = AdaptiveProposalScale[isadp, Min[ind], adprate, adpset];
        adpset[[-1]] = Min[adpset[[-1]], 1./Length[groupmcstate]];
        If[ isadp,
            groupmcstate[[All, 1]] = 1 - Last[adpset]  Range[Length[groupmcstate]-1,0, -1];
        ];
        {groupmcstate, adpset}
    ]       
            
updateJIBDchain[chainstate_, nsq_, nbp_, prirho_,pritheta_,prieps_, maxcpden_,isphased_,isadp_, adprate_, rjfreq_,isprint_] :=
    Module[ {prithetashape, prithetascale, prithetamin, prithetamax, prirhoshape, prirhoscale, prirhomin, prirhomax, 
        priepsalpha, priepsbeta, priepsmin, priepsmax, heat, logl, theta, epsilon, rho, ibdls, 
        adpeps, adptheta, adprho,indibd, segarlist,snpls, direct},
        {heat, logl, {theta, epsilon, rho}, ibdls, {adpeps, adptheta, adprho,indibd,segarlist}, snpls} = chainstate;
        {prithetashape, prithetascale, prithetamin, prithetamax} = pritheta;
        {prirhoshape, prirhoscale, prirhomin, prirhomax} = prirho;
        {priepsalpha, priepsbeta, priepsmin, priepsmax} = prieps;
        direct = RandomChoice[{1,-1}];
        If[ direct == -1,
            {snpls, ibdls} = reversechromosome[snpls, ibdls, theta, epsilon, nbp, nsq];
        ];
        If[ isprint,
            Print["update epsilon..."]
        ];
        {snpls, epsilon, adpeps} = updateepsilon[snpls, ibdls, epsilon, nbp, priepsalpha, priepsbeta,priepsmin,
           priepsmax, isadp, adprate, adpeps,heat];
        If[ isprint,
            Print["update theta..."]
        ];
        {ibdls, theta, adptheta} = updatetheta[ibdls, theta, nsq, prithetashape,prithetascale, prithetamin,prithetamax,isadp, adprate,adptheta,heat];
        If[ isprint,
            Print["update ibdls..."]
        ];
        {snpls, ibdls, indibd,segarlist} = rjupdateibdls[snpls, ibdls, theta, epsilon, rho,nsq, nbp, prirhoshape, 
            prirhoscale,prirhomin,prirhomax, maxcpden,isphased,rjfreq,heat,isprint];
        If[ isprint,
            Print["update rho..."]
        ];
        {rho,adprho} = updaterho[ibdls, rho,nbp, prirhoshape, prirhoscale,prirhomin,prirhomax,maxcpden,isadp, adprate, adprho, heat];
        If[ direct == -1,
            {snpls, ibdls} = reversechromosome[snpls, ibdls, theta, epsilon, nbp, nsq];
        ];
        logl = Total[snpls[[All, -1]]]+Total[ibdls[[All, -1]]];
        logl += (Length[ibdls] - 1) Log[rho] - rho nbp+
                (priepsalpha - 1) Log[epsilon] + (priepsbeta - 1) Log[1 - epsilon]+
                (prithetashape-1) Log[theta]- theta/prithetascale + 
                (prirhoshape - 1) Log[rho] -rho/prirhoscale;
        If[ maxcpden<Infinity,
            logl+=-(Log[Gamma[1+Ceiling[maxcpden nbp/10^6],rho nbp]]-LogGamma[1+Ceiling[maxcpden nbp/10^6]])
        ];
        {heat, logl, {theta, epsilon, rho},ibdls, {adpeps, adptheta,adprho,indibd,segarlist}, snpls}
    ]            

reversechromosome[inputsnpls_, inputibdls_, theta_, epsilon_, nbp_,nsq_] :=
    Module[ {snpls = inputsnpls, ibdls = inputibdls, pri, test},
        test = RandomReal[] < 0.1;
        If[ test&&(! testsnpls[snpls, ibdls, epsilon, nbp]),
            Print["wrong snpls before reverse chromosome!"]
        ];
        If[ test&&(! testibdls[ibdls, theta, nsq]),
            Print["wrong ibdls before reverse chromosome!"]
        ];
        If[ snpls =!= {},
            snpls[[All, 1]] = nbp + 1 - snpls[[All, 1]];
            snpls = Reverse[snpls];
        ];
        ibdls[[All, 1]] = nbp + 1 - Append[Rest[ibdls[[All, 1]]] - 1, nbp];
        pri = LogPDFCRP[ibdls[[All, 2]], nsq, 0, theta];
        ibdls[[All, -1]] = Append[ibdls[[2 ;;, -1]] - Differences[pri], Last[pri]];
        ibdls = Reverse[ibdls];
        If[ test&&(! testsnpls[snpls, ibdls, epsilon, nbp]),
            Print["wrong snpls after reverse chromosome!"]
        ];
        If[ test&&(! testibdls[ibdls, theta, nsq]),
            Print["wrong ibdls after chromosome!"]
        ];
        {snpls, ibdls}
    ]   
    
testsnpls[snpls_, ibdls_, epsilon_,nbp_] :=
    Module[ {ls, ibdpos, temp},
        If[ snpls == {},
            True,
            ls = Append[ibdls[[All, 1]], nbp + 1];
            ibdpos = IndexByInterpolation[ls][snpls[[All, 1]]];
            temp = IBDLogLikelihood[snpls[[All, 2]], ibdls[[ibdpos, 2]],snpls[[All, 3]], epsilon];
            If[ Abs[Total[temp-snpls[[All, -1]]]]>10^(-10.),
                Print["{snplstrue,snpls}=", Transpose[{temp, snpls[[All, -1]]}]];
                False,
                True
            ]
        ]
    ]  

testibdls[ibdls_, theta_,nsq_] :=
    Module[ {temp},
        temp = If[ Length[ibdls] == 1,
                   {},
                   Log[Flatten[IBDTransitionProb[#, nsq, theta] & /@Partition[ibdls[[All, 2]], 2, 1]]]
               ];
        temp = Prepend[temp, LogPDFCRP[ibdls[[1, 2]], nsq, 0, theta]];
        If[ Abs[Total[temp-ibdls[[All, -1]]]]>10^(-10.),
            Print["{ibdlstrue,ibdls}=",Transpose[{temp, ibdls[[All, -1]]}]];
            False,
            True
        ]
    ]    
    
updateepsilon[inputsnpls_, ibdls_, inputepsilon_, nbp_, priepsalpha_, 
   priepsbeta_, priepsmin_,priepsmax_, isadpt_, adprate_, inputadpset_,heat_] :=
    Module[ {snpls = inputsnpls, epsilon = inputepsilon, 
     adpset = inputadpset, size = 1, ind = 0, logpri, ls2, map, zls, 
     deltxls, propx, propsnpls, ar,i},
        logpri = Function[{x}, (priepsalpha - 1) Log[x] + (priepsbeta - 1) Log[1 - x]];
        If[ snpls=!={},
            ls2 = Append[ibdls, ReplacePart[ibdls[[-1]], 1 -> nbp + 1]];
            map = IndexByInterpolation[ls2[[All, 1]]];
            zls = ls2[[map[snpls[[All, 1]]], 2]];
        ];
        deltxls = RandomReal[NormalDistribution[0, adpset[[-1]]], size];
        Do[
         propx = epsilon + deltxls[[i]];
         If[ priepsmin <= propx <= priepsmax,
             propsnpls = snpls;
             If[ propsnpls==={},
                 ar = 0,
                 propsnpls[[All, -1]] = IBDLogLikelihood[propsnpls[[All, 2]], zls,propsnpls[[All, 3]], propx];
                 ar = Total[propsnpls[[All, -1]]] - Total[snpls[[All, -1]]];
             ];
             ar += logpri[propx] - logpri[epsilon];
             ar = Min[1, Exp[heat ar]];
             If[ RandomReal[] < ar,
                 epsilon = propx;
                 snpls = propsnpls;
             ],
             ar = 0
         ];
         ind += ar, {i, size}];
        adpset = AdaptiveProposalScale[isadpt, ind/size, adprate, adpset];
        {snpls, epsilon, adpset}
    ]

updatetheta[ibdls_, theta_, nsq_,prithetashape_,prithetascale_, prithetamin_,prithetamax_,isadp_, adprate_,inputadpset_,heat_] :=
    Module[ {adpset = inputadpset,logpri, propx, propibdls, ar},
        logpri = Function[{x}, (prithetashape-1) x-Exp[x]/prithetascale+x];
        propx = RandomReal[NormalDistribution[Log[theta], Last[adpset]]];
        If[ prithetamax> Exp[propx] > prithetamin,
            propibdls = ibdls;
            propibdls[[All, -1]] = Prepend[If[ Length[propibdls]==1,
                                               {},
                                               Log[IBDTransitionProb[propibdls[[All, 2]], nsq, Exp[propx]]]
                                           ],
                LogPDFCRP[propibdls[[1, 2]], nsq, 0, Exp[propx]]];
            ar = Total[propibdls[[All, -1]]] - Total[ibdls[[All, -1]]];
            ar += logpri[propx] - logpri[Log[theta]];
            ar = Min[1, Exp[heat ar]],
            ar = 0
        ];
        adpset = AdaptiveProposalScale[isadp, ar, adprate, adpset];
        If[ RandomReal[] < ar,
            {propibdls, Exp[propx], adpset},
            {ibdls, theta, adpset}
        ]
    ]
         
updaterho[ibdls_, rho_, nbp_, prirhoshape_, prirhoscale_, prirhomin_, prirhomax_,maxcpden_,isadp_, adprate_, adprho_, heat_] :=
    Module[ {a, b, logpdf, res},
        a = heat (Length[ibdls] + prirhoshape - 2) + 1;
        b = heat (nbp + 1./prirhoscale);
        If[ maxcpden ===prirhomax=== Infinity&&prirhomin==0,
            {RandomReal[GammaDistribution[a, 1/b]], adprho},
            logpdf = Function[{x}, (a - 1) Log[x] - x b - heat Log[Gamma[1 + Ceiling[maxcpden nbp/10^6], x nbp]]];
            res = AdaptiveMetropolis[logpdf, rho, 1, 1, Boole[isadp], adprate, adprho, DomainConstraint -> (prirhomin<#<prirhomax&)];
            res[[{1, -1}]]
        ]
    ]    
    
    
rjupdateibdls[inputsnpls_, inputibdls_, theta_, epsilon_,rho_, nsq_, nbp_, prirhoshape_, prirhoscale_,
    prirhomin_,prirhomax_,maxcpden_,isphased_,rjfreq_,heat_,isprint_] :=
    Module[ {snpls = inputsnpls, ibdls = inputibdls, cc, weight,actls, ar,act,ind,segarlist,seglen},
        weight = {(1 - 2 cc)/3,(1-2 cc)/3,(1-2 cc)/3,cc,cc};
        ind = Table[{0, 10^(-6.)}, {Length[weight]+2}]; 
        cc=If[maxcpden<Infinity, 0.125,0.2];
        actls = RandomChoice[weight -> Range[Length[weight]],rjfreq];
        segarlist = Reap[
            Do[
                If[ isprint,
                    Print["act=",act]
                ];
                Switch[act,
                 1, {snpls, ibdls, ar} = updateibdpos[snpls, ibdls, epsilon, nsq, nbp,heat],
                 2, {snpls, ibdls, ar} = updateibdibd[snpls, ibdls, theta, epsilon, nsq, nbp,heat],
                 3, {snpls, ibdls, {seglen, ar}} = updatesegibdibd[snpls, ibdls, theta, epsilon,nsq, nbp,heat],
                 4, {snpls, ibdls, {seglen, ar}} = updatesegibdinsert[snpls, ibdls, theta, epsilon,rho, nsq, nbp,prirhoshape, prirhoscale,prirhomax,maxcpden,heat],
                 5, {snpls, ibdls, {seglen, ar}} = updatesegibddelete[snpls, ibdls, theta, epsilon,rho, nsq, nbp,prirhoshape, prirhoscale,prirhomax,maxcpden,heat]                 
                 ];
                If[ MemberQ[{3,4,5},act],
                    Sow[{seglen,ar}]
                ];
                ind[[act]] += {ar, 1}, {act, actls}];
        ];
        {snpls, ibdls, ar} = updatechrswap[snpls, ibdls, epsilon, nsq, nbp,heat];
        ind[[Length[weight]+1]] += {ar, 1};
        If[ !isphased,
            {snpls, ar} = updatephasesingle[snpls, ibdls, epsilon, nsq, nbp, heat];
            ind[[-1]] += {ar, 1};
            {snpls, ibdls} = updatephaseblock[snpls, ibdls, nsq, nbp];
        ];        
        segarlist = If[ segarlist[[2]]=!={},
                        Mean[#] & /@ SplitBy[SortBy[segarlist[[2, 1]], First], First],
                        {}
                    ];
        {snpls, ibdls, ind[[All, 1]]/ind[[All, 2]],segarlist}
    ]    
            


updateibdpos[snpls_, ibdls_, epsilon_, nsq_, nbp_,heat_] :=
    Module[ {foc,ar, s1, s2, s3, props, propibdls, propsnpls, ibd1, ibd2},
        If[ Length[ibdls] == 1,
            ar = 0,
            foc = RandomInteger[{2, Length[ibdls]}];
            {{s1, ibd1}, {s2, ibd2}} = ibdls[[foc - 1 ;; foc, {1, 2}]];
            s3 = If[ foc == Length[ibdls],
                     nbp + 1,
                     ibdls[[foc + 1, 1]]
                 ];
            props = RandomInteger[{s1 + 1, s3 - 1}];
            propibdls = ibdls;
            propibdls[[foc, 1]] = props;
            {propsnpls, ar} = If[ s2 < props,
                                  reviseloglike[snpls, ibd1, s2, props, epsilon],
                                  reviseloglike[snpls, ibd2, props, s2, epsilon]
                              ];
            ar = Min[1, Exp[heat ar]]
        ];
        If[ RandomReal[] < ar,
            {propsnpls, propibdls, ar},
            {snpls, ibdls, ar}
        ]
    ]


updateibdibd[snpls_, ibdls_, theta_, epsilon_, nsq_, nbp_,heat_] :=
    Module[ {foc,propibd, propibdls, propsnpls, ar, ar2, z1, z2,z3, ops, type},
        foc = RandomInteger[{1, Length[ibdls]}];
        Which[
         foc == Length[ibdls] == 1,
         propibd = IBDTransition[ibdls[[1, 2]], nsq, theta];
         If[ propibd==ibdls[[1,2]],
             propibdls = ibdls;
             propsnpls = snpls;
             ar = 1,
             propibdls = ibdls;
             propibdls[[1, 2 ;;]] = {propibd, LogPDFCRP[propibd, nsq, 0, theta]};
             {propsnpls, ar} = reviseloglike[snpls, propibdls[[1, 2]], ibdls[[1, 1]], nbp + 1,epsilon];
             ar += propibdls[[1, -1]] - ibdls[[1, -1]];
             ar2 = Log[IBDTransitionProb[{ibdls[[1, 2]], propibd}, nsq, theta][[1]]] - 
               Log[IBDTransitionProb[{propibd, ibdls[[1, 2]]}, nsq, theta][[1]]];
             ar = Min[1, Exp[heat ar - ar2]]
         ],
         foc == 1,
         propibd = IBDTransition[ibdls[[2, 2]], nsq, theta];
         propibdls = ibdls;
         propibdls[[1, 2 ;;]] = {propibd, LogPDFCRP[propibd, nsq, 0, theta]};
         propibdls[[2, -1]] = Log[Total[IBDTransitionProb[propibdls[[;; 2, 2]], nsq, theta]]];
         {propsnpls, ar} = reviseloglike[snpls, propibdls[[1, 2]], ibdls[[1, 1]],ibdls[[2, 1]], epsilon];
         ar += Total[propibdls[[;; 2, -1]]] - Total[ibdls[[;; 2, -1]]];
         ar2 = Log[IBDTransitionProb[{ibdls[[2, 2]], propibd}, nsq, theta][[1]]] - 
           Log[IBDTransitionProb[{ibdls[[2, 2]], ibdls[[1, 2]]}, nsq, theta][[1]]];
         ar = Min[1, Exp[heat ar - ar2]],
         foc == Length[ibdls],
         propibd = IBDTransition[ibdls[[-2, 2]], nsq, theta];
         propibdls = ibdls;
         propibdls[[-1, 2]] = propibd;
         propibdls[[-1, 3]] = Log[Total[IBDTransitionProb[propibdls[[-2 ;;, 2]], nsq, theta]]];
         {propsnpls, ar} = reviseloglike[snpls, propibdls[[-1, 2]], ibdls[[-1, 1]], nbp + 1, epsilon];
         ar += propibdls[[-1, -1]] - ibdls[[-1, -1]];
         ar2 = propibdls[[-1, -1]] - ibdls[[-1, -1]];
         ar = Min[1, Exp[heat ar-ar2]],
         True,
         {z1, z2,z3} = ibdls[[foc - 1;;foc + 1, 2]];
         If[ z1 == z3,
             propibd = IBDTransition[z1, nsq, theta];
             ar2 = Log[IBDTransitionProb[{z1, propibd}, nsq, theta][[1]]] - 
               Log[IBDTransitionProb[{z1, ibdls[[foc, 2]]}, nsq, theta][[1]]],
             type = IBDTransitionType1[{z1, z3}];
             If[ First[type] == -1,
                 propibd = randommidibd2[z1, z3];
                 ar2 = 0,
                 ops = getmidoperators[type];
                 propibd = Switch[RandomChoice[{0.25,0.25,0.5}->{1,2,3}],
                     1,z1,
                     2,z3,
                     3, randommidibd1[z1,theta,ops]                 
                 ];
                 ar2 = Log[0.25 KroneckerDelta[propibd,z1]+0.25 KroneckerDelta[propibd,z3]+0.5 pdfmidibd1[z1, propibd, nsq, theta, ops]]-
                     Log[0.25 KroneckerDelta[z2,z1]+0.25 KroneckerDelta[z2,z3]+0.5 pdfmidibd1[z1, z2, nsq, theta, ops]]
             ]
         ];
         propibdls = ibdls;
         propibdls[[foc, 2]] = propibd;
         propibdls[[foc ;; foc + 1, -1]] = Log[IBDTransitionProb[propibdls[[foc - 1 ;; foc + 1, 2]], nsq, theta]];
         {propsnpls, ar} = reviseloglike[snpls, propibdls[[foc, 2]], ibdls[[foc, 1]],ibdls[[foc + 1, 1]], epsilon];
         ar += Total[propibdls[[foc ;; foc + 1, -1]]] - Total[ibdls[[foc ;; foc + 1, -1]]];
         ar = Min[1, Exp[heat ar - ar2]]
        ];
        If[ RandomReal[] < ar,
            {propsnpls, propibdls, ar},
            {snpls, ibdls,ar}
        ]
    ]    
    

    
updateibdinsert[snpls_, ibdls_, theta_, epsilon_,rho_, nsq_, nbp_,prirhoshape_, prirhoscale_,prirhomax_,princpmax_,heat_] :=
    Module[ {foc, s1, z1, props, propibd, propibdls, propsnpls, ar, ar2,s2, z2, type, ops},
        If[ Length[ibdls]-1>Ceiling[princpmax nbp/10^6]-1,
            ar = 0,
            foc = RandomInteger[{1, Length[ibdls]}];
            {s1, z1} = ibdls[[foc, ;; 2]];
            If[ foc == Length[ibdls],
                If[ nbp - s1 < 1,
                    ar = 0,
                    props = RandomInteger[{s1 + 1, nbp}];
                    propibd = IBDTransition[z1, nsq, theta];
                    propibdls = Insert[ibdls, {props, propibd, Log[IBDTransitionProb[{z1, propibd}, nsq, theta][[1]]]},foc + 1];
                    {propsnpls, ar} = reviseloglike[snpls, propibd, props,nbp + 1, epsilon];
                    ar += propibdls[[-1, -1]];
                    ar = heat ar+logpriKx[Length[propibdls], rho,princpmax, prirhoshape, prirhoscale,prirhomax, nbp,heat]-
                            logpriKx[Length[ibdls], rho,princpmax, prirhoshape, prirhoscale,prirhomax, nbp,heat];
                    ar2 = propibdls[[-1, -1]] + Log[1./(nbp - s1)];
                    ar = Min[1, Exp[ar - ar2]]
                ],
                {s2, z2} = ibdls[[foc + 1, ;; 2]];
                If[ s2 - s1 < 2,
                    ar = 0,
                    props = RandomInteger[{s1 + 1, s2 - 1}];
                    type = IBDTransitionType1[{z1, z2}];
                    If[ First[type]==-1,
                        Print["Wrong ibds in ibdls {foc,z1,z2}=", {foc, z1, z2}];
                        ar = 0,
                        If[ First[type]==0,
                            propibd = IBDTransition[z1, nsq, theta];
                            ar2 = Log[IBDTransitionProb[{z1, propibd}, nsq, theta][[1]]] + Log[1./(s2 - s1 - 1)],
                            ops = getmidoperators[type];
                            propibd = randommidibd1[z1,theta,ops];
                            ar2 = Log[pdfmidibd1[z1, propibd, nsq, theta, ops]] +Log[1./(s2 - s1 - 1)];
                        ];
                        propibdls = Insert[ibdls, {props, propibd, Log[IBDTransitionProb[{z1, propibd}, nsq, theta][[1]]]},foc + 1];
                        propibdls[[foc + 2, -1]] = Log[IBDTransitionProb[{propibd, z2}, nsq, theta][[1]]];
                        {propsnpls, ar} = reviseloglike[snpls, propibd, props, s2, epsilon];
                        ar += Total[propibdls[[foc + 1 ;; foc + 2, -1]]]-ibdls[[foc + 1, -1]];
                        ar = heat ar+ logpriKx[Length[propibdls], rho,princpmax, prirhoshape, prirhoscale,prirhomax, nbp,heat]-
                            logpriKx[Length[ibdls], rho,princpmax, prirhoshape, prirhoscale,prirhomax, nbp,heat];
                        ar = Min[1, Exp[ar - ar2]]
                    ];
                ];
            ];
        ];
        If[ RandomReal[] < ar,
            {propsnpls, propibdls, ar},
            {snpls, ibdls, ar}
        ]
    ]   
  
updateibddelete[snpls_, ibdls_, theta_, epsilon_, rho_,nsq_, nbp_, prirhoshape_, prirhoscale_,prirhomax_,princpmax_,heat_,isprint_:False] :=
    Module[ {ar, foc, propibdls, propsnpls, ar2, s1, s2, s3, z1, z2, z3,type},
        If[ Length[ibdls] == 1,
            ar = 0,
            foc = RandomInteger[{2, Length[ibdls]}];
            (*If[ isprint,
                Print["{Length[ibdls],foc}=",{Length[ibdls],foc}]
            ];*)
            If[ foc == Length[ibdls],
                propibdls = Delete[ibdls, foc];
                {propsnpls, ar} = reviseloglike[snpls, ibdls[[-2, 2]], ibdls[[-1, 1]], nbp + 1, epsilon];
                ar += 0 - ibdls[[-1, -1]];
                ar = heat ar+logpriKx[Length[propibdls], rho,princpmax, prirhoshape, prirhoscale,prirhomax, nbp,heat] - 
                  logpriKx[Length[ibdls], rho,princpmax, prirhoshape, prirhoscale,prirhomax, nbp,heat];
                (*If[ isprint,
                    Print["{Length[ibdls],foc,ibdls[[-1,1]]}=",{Length[ibdls],foc,ibdls[[-1,1]]}]
                ];*)
                ar2 = 0 - (ibdls[[-1, -1]] + Log[1./(nbp - ibdls[[-2, 1]])]);
                ar = Min[1, Exp[ar - ar2]],
                {{s1, z1}, {s2, z2}, {s3, z3}} = ibdls[[foc - 1 ;; foc + 1, ;; 2]];
                type = IBDTransitionType1[{z1, z3}];
                If[ First[type]==-1,
                    ar = 0,
                    ar2 = If[ First[type]==0,
                              -(ibdls[[foc,-1]]+Log[1./(s3 - s1-1)]),
                              -(Log[pdfmidibd1[z1,z2, nsq, theta, getmidoperators[type]]]+ Log[1./(s3 - s1-1)])
                          ];
                    (*If[ isprint,
                        Print["{s1,s3,ar2}=",{s1,s3,ar2}];
                        Print["{z1,z2,z3,nsq,theta,type}=",{z1,z2,z3,nsq,theta,type}];
                    ];*)
                    propibdls = Delete[ibdls, foc];
                    propibdls[[foc, -1]] = Log[IBDTransitionProb[{z1, z3}, nsq, theta][[1]]];
                    {propsnpls, ar} = reviseloglike[snpls, z1, s2, s3, epsilon];
                    ar += propibdls[[foc, -1]] - Total[ibdls[[foc;;foc+1, -1]]];
                    ar = heat ar+logpriKx[Length[propibdls], rho,princpmax, prirhoshape, prirhoscale,prirhomax, nbp,heat] - 
                      logpriKx[Length[ibdls], rho,princpmax, prirhoshape, prirhoscale,prirhomax, nbp,heat];
                    ar = Min[1, Exp[ar - ar2]]
                ];
            ];
        ];
        (*Print[{heat,ar,ar2}];*)
        If[ RandomReal[] < ar,
            {propsnpls, propibdls, ar},
            {snpls, ibdls,ar}
        ]
    ]  


    
updatesegibdibd[snpls_, ibdls_, theta_, epsilon_, nsq_, nbp_,heat_] :=
    Module[ {foc, propibd, propibdls, propsnpls, ar, ar2, z1, seglen = 1, 
      propsub, ar2sub, maxpos},
        foc = RandomChoice[ReplacePart[Table[1, {Length[ibdls]}], 1 -> nsq + 1] ->Range[Length[ibdls]]];
        (*foc = RandomInteger[{1, Length[ibdls]}];*)
        Which[
         Length[ibdls] == foc == 1,
         propibd = IBDTransition[ibdls[[1, 2]], nsq, theta];
         If[ propibd == ibdls[[1, 2]],
             propibdls = ibdls;
             propsnpls = snpls;
             ar = 1,
             propibdls = ibdls;
             propibdls[[1, 2 ;;]] = {propibd, LogPDFCRP[propibd, nsq, 0, theta]};
             {propsnpls, ar} = reviseloglike[snpls, propibdls[[1, 2]], ibdls[[1, 1]], nbp + 1, epsilon];
             ar += propibdls[[1, -1]] - ibdls[[1, -1]];
             ar2 = Log[IBDTransitionProb[{ibdls[[1, 2]], propibd}, nsq, theta][[1]]] -
                Log[IBDTransitionProb[{propibd, ibdls[[1, 2]]}, nsq, theta][[1]]];
             ar = Min[1, Exp[heat ar - ar2]]
         ],
         foc == Length[ibdls],
         propibd = IBDTransition[ibdls[[-2, 2]], nsq, theta];
         propibdls = ibdls;
         propibdls[[-1, 2 ;;]] = {propibd, Log[IBDTransitionProb[{ibdls[[-2, 2]], propibd}, nsq, theta][[1]]]};
         {propsnpls, ar} = reviseloglike[snpls, propibdls[[-1, 2]], ibdls[[-1, 1]], nbp + 1,epsilon];
         ar += propibdls[[-1, -1]] - ibdls[[-1, -1]];
         ar2 = propibdls[[-1, -1]] - ibdls[[-1, -1]];
         ar = Min[1, Exp[heat ar-ar2]],
         True,
         If[ foc == 1,
             z1 = ibdls[[1, 2]];
             propibd = IBDTransition[z1, nsq, theta];
             ar2 = If[ propibd == z1,
                       0,
                       Log[IBDTransitionProb[{z1, propibd}, nsq, theta][[1]]] - 
                        Log[IBDTransitionProb[{propibd, z1}, nsq, theta][[1]]]
                   ],
             z1 = ibdls[[foc, 2]];
             {ar2, propibd} = randomcondibd[ibdls[[foc - 1, 2]], z1, nsq, theta];
         ];
         If[ propibd == z1,
             propibdls = ibdls;
             propsnpls = snpls;
             ar = 1,
             propibdls = ibdls;
             propibdls[[foc, 2 ;;]] = {propibd, If[ foc == 1,
                                                    LogPDFCRP[propibd, nsq, 0, theta],
                                                    Log[IBDTransitionProb[{propibdls[[foc - 1, 2]], propibd}, nsq,theta][[1]]]
                                                ]};
             {ar2sub, propsub} = randomsegcondibd[ibdls[[foc ;;, 2]], propibdls[[foc, 2]]];
             If[ ar2sub === None,
                 seglen = -1;
                 ar = 0,
                 ar2 += ar2sub;
                 seglen = Length[propsub];
                 maxpos = foc + seglen - 1;
                 (*Print["seglen=",seglen,"; last==?",propsub[[-1]]==ibdls[[maxpos,2]]];*)
                 propibdls[[foc ;; maxpos, 2]] = propsub;
                 propibdls[[foc + 1 ;; maxpos, -1]] = Log[IBDTransitionProb[propsub, nsq, theta]];
                 {propsnpls, ar} = reviseloglike2[snpls, propibdls, foc, maxpos, epsilon, nbp];
                 ar += Total[propibdls[[foc ;; maxpos, -1]]]-Total[ibdls[[foc ;; maxpos, -1]]];
                 ar = Min[1, Exp[heat ar - ar2]];
             ];
         ];
         ];
        If[ RandomReal[] < ar,
            {propsnpls, propibdls, {seglen, ar}},
            {snpls,ibdls, {seglen, ar}}
        ]
    ]           
    
updatesegibdibd2[snpls_, ibdls_, theta_, epsilon_, nsq_, nbp_,heat_] :=
    Module[ {foc, propibd, propibdls, propsnpls, ar, ar2,seglen = 1},
        (*foc = RandomChoice[ReplacePart[Table[1, {Length[ibdls]}], 1 -> nsq + 1] ->Range[Length[ibdls]]];*)
        foc = RandomInteger[{1, Length[ibdls]}];
        Which[
         Length[ibdls] == foc == 1,
         propibd = IBDTransition[ibdls[[1, 2]], nsq, theta];
         If[ propibd == ibdls[[1, 2]],
             propibdls = ibdls;
             propsnpls = snpls;
             ar = 1,
             propibdls = ibdls;
             propibdls[[1, 2 ;;]] = {propibd, LogPDFCRP[propibd, nsq, 0, theta]};
             {propsnpls, ar} = reviseloglike[snpls, propibdls[[1, 2]], ibdls[[1, 1]], nbp + 1, epsilon];
             ar += propibdls[[1, -1]] - ibdls[[1, -1]];
             ar2 = Log[IBDTransitionProb[{ibdls[[1, 2]], propibd}, nsq, theta][[1]]] -
                Log[IBDTransitionProb[{propibd, ibdls[[1, 2]]}, nsq, theta][[1]]];
             ar = Min[1, Exp[heat ar - ar2]]
         ],
         foc == Length[ibdls],
         propibd = IBDTransition[ibdls[[-2, 2]], nsq, theta];
         propibdls = ibdls;
         propibdls[[-1, 2 ;;]] = {propibd, Log[IBDTransitionProb[{ibdls[[-2, 2]], propibd}, nsq, theta][[1]]]};
         {propsnpls, ar} = reviseloglike[snpls, propibdls[[-1, 2]], ibdls[[-1, 1]], nbp + 1,epsilon];
         ar += propibdls[[-1, -1]] - ibdls[[-1, -1]];
         ar2 = propibdls[[-1, -1]] - ibdls[[-1, -1]];
         ar = Min[1, Exp[heat ar-ar2]],
         True,
         {ar,seglen,propsnpls,propibdls} = randomsegprop[foc, snpls, ibdls, theta, epsilon, nsq, nbp,heat];        
        ];
        If[ RandomReal[] < ar,
            {propsnpls, propibdls, {seglen, ar}},
            {snpls,ibdls, {seglen, ar}}
        ]
    ] 

    
updatesegibdinsert[snpls_, ibdls_, theta_, epsilon_,rho_, nsq_, nbp_,prirhoshape_, prirhoscale_,prirhomax_,princpmax_,heat_] :=
    Module[ {seglen = 1, foc, seg0, props, ar, propibd, propibdls,propsnpls, ar2, ar2sub, propsub, maxpos},
        If[ Length[ibdls]-1>Ceiling[princpmax nbp/10^6]-1,
            ar = 0,
            foc = RandomInteger[{1, Length[ibdls]}];
            seg0 = {ibdls[[foc, 1]] + 1, If[ foc == Length[ibdls],
                                             nbp,
                                             ibdls[[foc + 1, 1]] - 1
                                         ]};
            props = RandomInteger[seg0];
            Which[
             Greater @@ seg0,
             ar = 0,
             foc == Length[ibdls],
             propibd = IBDTransition[ibdls[[foc, 2]], nsq, theta];
             propibdls = Insert[ibdls, {props, propibd,Log[IBDTransitionProb[{ibdls[[foc, 2]], propibd}, nsq,theta][[1]]]}, foc + 1];
             {propsnpls, ar} = reviseloglike[snpls, propibd, props, nbp + 1, epsilon];
             ar += propibdls[[foc + 1, -1]];
             ar = heat ar+logpriKx[Length[propibdls], rho,princpmax,prirhoshape, prirhoscale,prirhomax, nbp,heat] - 
               logpriKx[Length[ibdls], rho,princpmax, prirhoshape, prirhoscale,prirhomax, nbp,heat];
             ar2 = propibdls[[foc + 1, -1]] + Log[1./(seg0[[2]] - seg0[[1]] + 1)];
             ar = Min[1, Exp[ar - ar2]],
             True,
             propibd = IBDTransition[ibdls[[foc, 2]], nsq, theta];
             propibdls = Insert[ibdls, {props, propibd, Log[IBDTransitionProb[{ibdls[[foc, 2]], propibd}, nsq,theta][[1]]]}, foc + 1];
             ar2 = propibdls[[foc + 1, -1]] + Log[1./(seg0[[2]] - seg0[[1]] + 1)];
             If[ propibd == ibdls[[foc, 2]],
                 propsnpls = snpls;
                 ar = propibdls[[foc+1, -1]];
                 ar = heat ar+logpriKx[Length[propibdls], rho,princpmax, prirhoshape, prirhoscale,prirhomax, nbp,heat] - 
                   logpriKx[Length[ibdls], rho,princpmax, prirhoshape, prirhoscale,prirhomax, nbp,heat];
                 ar = Min[1, Exp[ar - ar2]],
                 {ar2sub, propsub} = randomsegcondibd[ibdls[[foc ;;, 2]], propibdls[[foc + 1, 2]]];
                 If[ ar2sub === None,
                     seglen = -1;
                     ar = 0,
                     ar2 += ar2sub;
                     seglen = Length[propsub];
                     maxpos = foc + seglen;
                     propibdls[[foc + 1 ;; maxpos, 2]] = propsub;
                     propibdls[[foc + 2 ;; maxpos, -1]] = Log[IBDTransitionProb[propsub, nsq, theta]];
                     {propsnpls, ar} = reviseloglike2[snpls, propibdls, foc + 1, maxpos, epsilon, nbp];
                     ar += Total[propibdls[[foc + 1 ;; maxpos, -1]]]-Total[ibdls[[foc + 1 ;; maxpos - 1, -1]]];
                     ar = heat ar+logpriKx[Length[propibdls], rho,princpmax, prirhoshape, prirhoscale,prirhomax, nbp,heat] - 
                       logpriKx[Length[ibdls], rho,princpmax, prirhoshape, prirhoscale,prirhomax, nbp,heat];
                     ar = Min[1, Exp[ar - ar2]];
                 ];
             ];
            ];
        ];
        If[ RandomReal[] < ar,
            {propsnpls, propibdls, {seglen, ar}},
            {snpls, ibdls, {seglen, ar}}
        ]
    ]    

updatesegibddelete[snpls_, ibdls_, theta_, epsilon_,rho_, nsq_, nbp_, prirhoshape_, prirhoscale_,prirhomax_, princpmax_,heat_] :=
    Module[ {seglen = 1, ar, foc, propibdls, seg0, propsnpls, ar2, ar2sub, propsub, maxpos},
        If[ Length[ibdls] == 1,
            ar = 0,
            foc = RandomInteger[{2, Length[ibdls]}];
            propibdls = Delete[ibdls, foc];
            seg0 = {ibdls[[foc - 1, 1]] + 1, If[ foc == Length[ibdls],
                                                 nbp,
                                                 ibdls[[foc + 1, 1]] - 1
                                             ]};
            ar2 = 0 - (ibdls[[foc, -1]] + Log[1/(seg0[[2]] - seg0[[1]] + 1)]);
            If[ foc == Length[ibdls],
                {propsnpls, ar} = reviseloglike[snpls, ibdls[[foc - 1, 2]], ibdls[[foc, 1]],nbp + 1, epsilon];
                ar += 0 - ibdls[[foc, -1]];
                ar = heat ar+logpriKx[Length[propibdls], rho,princpmax, prirhoshape, prirhoscale,prirhomax, nbp,heat] - 
                  logpriKx[Length[ibdls], rho,princpmax, prirhoshape, prirhoscale,prirhomax, nbp,heat];
                ar = Min[1, Exp[ar - ar2]],
                If[ ibdls[[foc - 1, 2]] == ibdls[[foc, 2]],
                    propsnpls = snpls;
                    ar = 0 - ibdls[[foc, -1]];
                    ar = heat ar+logpriKx[Length[propibdls], rho,princpmax, prirhoshape, prirhoscale,prirhomax, nbp,heat] - 
                      logpriKx[Length[ibdls], rho,princpmax, prirhoshape, prirhoscale,prirhomax, nbp,heat];
                    ar = Min[1, Exp[ar - ar2]],
                    {ar2sub, propsub} = randomsegcondibd[ibdls[[foc ;;, 2]], propibdls[[foc - 1, 2]]];
                    If[ ar2sub === None,
                        seglen = -1;
                        ar = 0,
                        ar2 += ar2sub;
                        seglen = Length[propsub];
                        maxpos = foc + seglen - 2;
                        propibdls[[foc - 1 ;; maxpos, 2]] = propsub;
                        propibdls[[foc ;; maxpos, -1]] = Log[IBDTransitionProb[propsub, nsq, theta]];
                        {propsnpls, ar} = reviseloglike2[snpls, propibdls, foc - 1, maxpos, epsilon,nbp];
                        ar += Total[propibdls[[foc ;; maxpos, -1]]] -Total[ibdls[[foc ;; maxpos + 1, -1]]];
                        ar = heat ar+ logpriKx[Length[propibdls], rho,princpmax, prirhoshape, prirhoscale,prirhomax, nbp,heat] - 
                          logpriKx[Length[ibdls], rho,princpmax, prirhoshape, prirhoscale,prirhomax, nbp,heat];
                        ar = Min[1, Exp[ar - ar2]];
                    ];
                ];
            ];
        ];
        If[ RandomReal[] < ar,
            {propsnpls, propibdls, {seglen, ar}},
            {snpls, ibdls, {seglen, ar}}
        ]
    ]

updatechrswap[inputsnpls_, inputibdls_, epsilon_, nsq_, nbp_,heat_] :=
    Module[ {ind = {0, 10^(-6.)}, snpls = inputsnpls, ibdls = inputibdls,ls,ibdpos,seg,size,
       pair,  bool, segpos, i1, i2, propibdls, propsnpls, xleft, xright, pos, ar},
        size = Ceiling[nsq/2.];
        ls = Append[ibdls[[All, 1]], nbp + 1];
        ibdpos = IndexByInterpolation[ls][snpls[[All, 1]]];
        Do[
            pair = Sort[RandomSample[Range[nsq], 2]];
            bool = (MemberQ[#, {pair[[1]]}]&&MemberQ[#, {pair[[2]]}])||MemberQ[#, {___,pair[[1]], ___, pair[[2]], ___}]& /@ ibdls[[All, 2]];
            segpos = SplitBy[Transpose[{Range[Length[bool]], bool}],Last][[All, {1, -1}]];
            segpos = {#[[All, 1]], #[[1, 2]]} & /@ segpos;
            segpos = Cases[segpos, {_, False}][[All, 1]];
            Do[
                {i1, i2} = seg;
                propibdls = ibdls;
                propibdls[[i1 ;; i2, 2]] = Map[Sort,propibdls[[i1 ;; i2, 2]] /. Thread[pair -> Reverse[pair]], 2];
                propsnpls = snpls;
                If[ propsnpls === {},
                    ar = 1,
                    xleft = ibdls[[i1, 1]];
                    xright = If[ i2 == Length[ibdls],
                                 nbp + 1,
                                 ibdls[[i2 + 1, 1]]
                             ];
                    pos = Flatten[Position[propsnpls[[All, 1]], _?(xleft <= # < xright &), {1},Heads -> False]];
                    pos = Pick[pos, Unequal @@ # & /@ propsnpls[[pos, 2, pair]]];
                    If[ pos === {},
                        ar = 1,
                        propsnpls[[pos, -1]] = IBDLogLikelihood[propsnpls[[pos, 2]],propibdls[[ibdpos[[pos]], 2]], 
                                                    propsnpls[[pos, 3]], epsilon];
                        ar = Total[propsnpls[[pos, -1]]] - Total[snpls[[pos, -1]]];
                        ar = Min[1, Exp[heat ar]];
                    ]
                ];
                ind += {ar, 1};
                If[ RandomReal[] <ar,
                    {snpls, ibdls} = {propsnpls, propibdls}
                ], {seg,segpos}],{size}];
        {snpls, ibdls, If[ segpos==={},
                           1,
                           Divide @@ ind
                       ]}
    ]    
         
updatephasesingle[inputsnpls_, ibdls_, epsilon_, nsq_,nbp_,heat_] :=
    Module[ {snpls = inputsnpls, temp, heterogeno, ls2, map, zls, 
      apt = {0, 10^(-6.)}, pair, zibd, chrs, ar,i,pairls,ils},
        temp = Partition[#, 2] & /@ snpls[[All, 2]];
        temp = Map[Equal @@ # &, temp, {2}];
        temp = Flatten[Position[#, False, {1}, Heads -> False]] & /@ temp;
        heterogeno = Transpose[{2 # - 1, 2 #}] & /@ temp;
        If[ snpls =!= {},
            ls2 = Append[ibdls, ReplacePart[ibdls[[-1]], 1 -> nbp + 1]];
            map = IndexByInterpolation[ls2[[All, 1]]];
            zls = ls2[[map[snpls[[All, 1]]], 2]];
        ];
        ils = Sort[RandomSample[Range[Length[snpls]], Round[Length[snpls]/2]]];
        (*ils=RandomSample[Range[Length[snpls]]];*)
        Do[
         pairls = Pick[heterogeno[[i]], RandomChoice[{True, False}, Length[heterogeno[[i]]]]];
         Do[
             zibd = Select[zls[[i]], Intersection[#, pair] != {} &];
             chrs = Union @@ zibd;
             If[ !(Length[chrs] == 2 || Length[zibd] == 1),
                 zibd = zibd /. Thread[chrs -> Range[Length[chrs]]];
                 ar = IBDLogLikelihood[snpls[[i, 2, chrs /. Thread[pair -> Reverse[pair]]]], zibd,snpls[[i, 3]], epsilon] -
                   IBDLogLikelihood[snpls[[i, 2, chrs]], zibd, snpls[[i, 3]], epsilon];
                 ar = Min[1, Exp[heat ar]];
                 If[ RandomReal[] < ar,
                     snpls[[i, 2, pair]] = Reverse[snpls[[i, 2, pair]]];
                     apt += {ar, 1};
                 ]
             ],{pair,pairls}], {i, ils}];
        snpls[[ils, -1]] = IBDLogLikelihood[snpls[[ils, 2]], zls[[ils]], snpls[[ils, 3]], epsilon];
        {snpls,Divide @@ apt}
    ]         

updatephaseblock[inputsnpls_, inputibdls_, nsq_,nbp_] :=
    Module[ {snpls = inputsnpls,ibdls = inputibdls,pairls,bool,segls,xleft,xright,pos,pair},
        pairls = Transpose[{Range[1, nsq, 2], Range[2, nsq, 2]}];
        pairls = Pick[pairls, RandomChoice[{True, False}, Length[pairls]]];
        (*pairls=RandomSample[pairls];*)
        Do[
          bool = (MemberQ[#, {pair[[1]]}] && MemberQ[#, {pair[[2]]}]) || 
              MemberQ[#, {___, pair[[1]], ___, pair[[2]], ___}] & /@ ibdls[[All, 2]];
          segls = SplitBy[Transpose[{Range[Length[bool]], bool}], Last][[All, {1, -1}]];
          segls = {#[[All, 1]], #[[1, 2]]} & /@ segls;
          segls = Pick[segls, RandomChoice[{True, False}, Length[segls]]];
          If[ segls =!= {},
              pos = Cases[segls, {_, False}][[All, 1]];
              pos = Flatten[Range @@ # & /@ pos];
              ibdls[[pos, 2]] = Map[Sort, ibdls[[pos, 2]] /. Thread[pair -> Reverse[pair]], 2];
              xleft = ibdls[[segls[[All, 1, 1]], 1]];
              xright = ibdls[[segls[[;; -2, 1, 2]] + 1, 1]];
              AppendTo[xright, If[ segls[[-1, 1, 2]] == Length[ibdls],nbp + 1,ibdls[[segls[[-1, 1, 2]] + 1, 1]]]];
              pos = Flatten[Table[Flatten[Position[snpls[[All, 1]], _?(xleft[[i]] <= # < xright[[i]] &), {1},Heads -> False]], {i, Length[xleft]}]];
              snpls[[pos, 2, pair]] = snpls[[pos, 2, Reverse[pair]]];
          ], {pair, pairls}];
        {snpls,ibdls}
    ]

logpriKx0[k_, prirhoshape_, prirhoscale_, prirhomax_,nbp_,heat_] :=
    Module[ {a,b},
        a = heat (k + prirhoshape-2)+1.;
        b = heat (nbp + 1./prirhoscale);
        -a Log[b]+
        If[ prirhomax==Infinity,
            LogGamma[a],
            If[ prirhomax> a/b,
                Log[Gamma[a]-Gamma[a,prirhomax b]],
                Log[Gamma[a,0,prirhomax b]]
            ]
        ]
    ]
    
logpriKx[k_, rho_, princpmax_, prirhoshape_, prirhoscale_, prirhomax_,nbp_, heat_] :=
    If[ princpmax == Infinity,
        logpriKx0[k, prirhoshape, prirhoscale, prirhomax, nbp, heat],
        heat ((k + prirhoshape - 2) Log[rho] - rho (nbp + 1/prirhoscale) - 
           Log[Gamma[1 + Ceiling[princpmax nbp/10^6], rho nbp]])
    ]    

reviseloglike[inputsnpls_, ibd_, xstart_, xend_, epsilon_] :=
    Module[ {pos, newsnpls = inputsnpls},
        pos = Flatten[Position[inputsnpls[[All, 1]], _?(xstart <= # < xend &), {1},Heads -> False]];
        If[ pos === {},
            {newsnpls, 0},
            newsnpls[[pos, -1]] = IBDLogLikelihood[newsnpls[[pos, 2]], Table[ibd, {Length[pos]}],newsnpls[[pos, 3]], epsilon];
            {newsnpls, Total[newsnpls[[pos, -1]]] - Total[inputsnpls[[pos, -1]]]}
        ]
    ]
    
(*update inputsnpls[[All,-1]] when ibdls[[istart;;iend]] were revised*)
reviseloglike2[inputsnpls_, ibdls_, istart_, iend_, epsilon_, nbp_] :=
    Module[ {snpls = inputsnpls, xstart, xend, pos, ls, ibds},
        xstart = ibdls[[istart, 1]];
        xend = If[ iend == Length[ibdls],
                   nbp + 1,
                   ibdls[[iend + 1, 1]]
               ];
        pos = Flatten[Position[snpls[[All, 1]], _?(xstart <= # < xend &), {1},Heads -> False]];
        If[ pos === {},
            {snpls, 0},
            ls = If[ iend == Length[ibdls],
                     Append[ibdls[[istart ;;]],ReplacePart[ibdls[[-1]], 1 -> nbp + 1]],
                     ibdls[[istart ;; iend + 1]]
                 ];
            ibds = ls[[IndexByInterpolation[ls[[All, 1]]][inputsnpls[[pos, 1]]], 2]];
            snpls[[pos, -1]] = IBDLogLikelihood[snpls[[pos, 2]], ibds, snpls[[pos, 3]],epsilon];
        ];
        {snpls, Total[snpls[[pos, -1]]] - Total[inputsnpls[[pos, -1]]]}
    ]    

(*type=IBDTransitionType1[{z1,z2}] statisfying type[[1]]=1 (d=1)*)
getmidoperators[type_] :=
    Module[ {ops},
        Switch[type[[2]],
         12,
         ops = If[ Length[type[[-1, 1]]] ==2,
                   {{type[[-1, 1]], -1, #} & /@Partition[type[[-1, 1]], 1]},
                   {{{type[[-1, 1]], -1, type[[-1, -1]]}}}
               ];
         If[ Length[type[[-1, 1]]] == 3,
             ops = Join[ops, {type[[-1, 1]], {}, #} & /@Partition[Complement @@ type[[-1, {1, -1}]], 1]]
         ],
         21,
         ops = If[ Length[type[[-1, -1]]] ==2,
                   {{#, -1, #} & /@Partition[type[[-1, -1]], 1]},
                   {{{type[[-1, 1]], -1,type[[-1, -1]]}}}
               ];
         If[ Length[type[[-1, 1]]] == 1 && Length[type[[-1, 2]]] == 2,
             ops = Join[ops, {type[[-1, 2]], type[[-1, 1]], #} & /@Partition[type[[-1, 2]], 1]]
         ],
         22,
         ops = {{{type[[-1, 1]], -1, type[[-1, -1]]}}};
         If[ Length[type[[-1, 1]]] == 2 && Length[type[[-1, 2]]] == 1,
             ops = Join[ops, {{type[[-1, 2]], type[[-1, 1]], type[[-1, 2]]},
                             {type[[-1,1]],type[[-1,2]],Complement @@ type[[-1, {1, -1}]]}}]
         ],
         _, Print["Wrong type in getmidoperators!", "type=", type]
         ];
        ops
    ]

(*return z by operating op on z1*)
transformibd[z1_, opp : {{_Integer ..}, {_Integer ...}, {_Integer..}}] :=
    Module[ {z = z1,op = opp},
        If[ op[[2]] === {},
            PrependTo[z, {}]
        ];
        Which[
            Length[op[[1]]] == 2&&op[[2]]==={},op[[-1]] = {op[[1,1]]},
            Length[op[[1]]] == 1&&Length[op[[2]]]==1,op[[-1]] = op[[1]]
        ];
        z = Replace[z, {op[[1]] -> Complement[op[[1]], op[[3]]], 
           op[[2]] -> Join[op[[2]], op[[3]]]}, {1}];
        Map[Sort, DeleteCases[z, {}], {0, 1}]
    ]

(*randomize middle tree given |z1-z2|=1,ops=getmidoperators[type], \
and type=IBDTransitionType1[{z1,z2}]*)
randommidibd1[z1_, theta_, ops_] :=
    Module[ {ind, op, z, w, pos},
        ind = RandomInteger[{1, Length[ops]}];
        op = If[ ind == 1,
                 op = RandomChoice[ops[[ind]]],
                 ops[[ind]]
             ];
        If[ op[[2]] === -1,
            z = Prepend[z1, {}];
            w = Length[#] & /@ z;
            w[[1]] = theta;
            pos = RandomChoice[w -> Range[Length[w]]];
            z = DeleteCases[z, op[[-1, 1]], {2}];
            z[[pos]] = Join[z[[pos]], op[[-1]]];
            Map[Sort, DeleteCases[z, {}], {0, 1}],
            transformibd[z1, op]
        ]
    ]

(*randomize middle tree given |z1-z2|=2*)
randommidibd2[z1_, z2_] :=
    Module[ {ops},
        ops = Last[IBDTransitionType2[{z1, z2}]];
        transformibd[z1, RandomChoice[ops]]
    ]
  
pdfmidibd1[z1_, z_, nsq_, theta_, ops_] :=
    Module[ {type, len},
        type = IBDTransitionType1[{z1, z}];
        Switch[First[type],
         0,
         len = Length[ops[[1, 1, 1]]];
         If[ len == 1,
             theta + 1,
             len
         ]/(nsq + theta)/Length[ops],
         1,
         If[ MemberQ[ops, type[[-1]]],
             1/Length[ops],
             If[ type[[2]] == 12,
                 theta/(nsq +theta),
                 (Length[type[[-1, 2]]]/(nsq + theta))* 
                 (Length[Intersection[type[[-1, -1]], ops[[1, All, -1, 1]]]]/Length[ops[[1]]])
             ]/Length[ops]
         ],
         _,
         Print["Wrong middle ibd in pdfmidibd!", 
          "{z1,z,ops}=", {z1, z, ops}];
         1
        ]
    ]

(*z1---z2 along chromosomes;
  |    |
 z1---(to sample z22)*)    
(*return {ar2,z22} if exist, otherwise {None,None} if failed*)
randomcondibd[z1_, z2_, nsq_, theta_] :=
    Module[ {type, z, ar2,ops,temp},
        type = IBDTransitionType1[{z1, z2}];
        Switch[First[type],
             0,
             (*Print["{z1,z2,theta}=",{z1,z2,theta}];*)
             z = IBDTransition[z1, nsq, theta];
             ar2 = If[ z==z1,
                       0,
                       temp = getmidoperators[IBDTransitionType1[{z1, z}]];
                       Log[IBDTransitionProb[{z1, z}, nsq, theta][[1]]]-
                       Log[pdfmidibd1[z1, z2, nsq, theta,temp]]
                   ];
             {ar2, z},
             1,
             (*Print["{z1,z2,theta}=",{z1,z2,theta}];*)
             ops = getmidoperators[type];
             z = randommidibd1[z1, theta,ops];
             ar2 = Which[
               z == z2, 0,
               z == z1, 
               Log[pdfmidibd1[z1, z, nsq, theta, ops]]-Log[IBDTransitionProb[{z1, z2}, nsq, theta][[1]]],
               True,
               temp = getmidoperators[IBDTransitionType1[{z1, z}]];
               Log[pdfmidibd1[z1, z, nsq, theta, ops]]-Log[pdfmidibd1[z1, z2, nsq, theta,temp]]
             ];
             {ar2, z},
             _,
             Print["Wrong in randomconibd!"];
             Print["{z1,z2}=", {z1, z2}];
             {None, None}
        ]
    ]    
 
(*z1---z2 along chromosomes;
  |    |
z11---(to sample z22)*)    
(*return {ar2,z22} if exist, otherwise {None,None} if failed*)
randomcondibd[z1_, z2_, z11_] :=
    Module[ {type, type2, ops,ops22,z22,type22,msg = "Wrong in randomcondibd!"},
        type = IBDTransitionType1[{z1, z2}];
        Switch[First[type],            
             0, 
             (*{0, z11}*)
             {0,z2},
             1,
             type2 = IBDTransitionType2[{z11, z2}];
             Switch[First[type2],
                  0, 
                  (*{0, z1}*)
                  {0,z2},
                  1, {0, z2},
                  2,
                  ops = Select[type2[[-1]],Intersection[#[[-1]], type[[-1, -1]]] =!= {} &];
                  (*Print["{z1,z2,z11}=",{z1,z2,z11}];*)
                  If[ ops === {},
                      Print["1: {z1,z2,z11}=", {z1, z2, z11}];
                      {None, None},
                      z22 = transformibd[z11, RandomChoice[ops]];
                      (*Print["{z11,z22}=",{z11,z22}];*)
                      type22 = IBDTransitionType2[{z1, z22}];
                      If[ First[type22]==2,
                          ops22 = Select[type22[[-1]],Intersection[#[[-1]], IBDTransitionType1[{z11, z22}][[-1, -1]]] =!= {} &];
                          If[ ops22==={},
                              Print["2: {z1,z2,z11}=", {z1, z2, z11}];
                              {None,None},
                              {Log[1./Length[ops]]-Log[1./Length[ops22]],z22}
                          ],
                          (*Print["3: {z1,z2,z11}=", {z1, z2, z11}];*)
                          {None,None}
                      ]
                  ],
                  _, 
                  Print[msg];
                  Print["{z1,z2,z11}=", {z1, z2, z11}]
             ],
             _, 
             Print[msg];
             Print["{z1,z2,z11}=", {z1, z2, z11}]
        ]
    ]

(*z1---z2--z3--... along chromosomes;
  |    |
z11---z22--z33--..*)
(*return {ar2,{z11,z22,z33,..}} if exists, otherwise {None,None} if failed*)        
randomsegcondibd[ibds_, z11_] :=
    Module[ {logprob = 0, propibds, logp, propz,i},
        propibds = Table[0, {Length[ibds]}];
        propibds[[1]] = z11;
        Do[
         If[ propibds[[i]] == ibds[[i]],
             propibds = Take[propibds, i];
             Break[]
         ];
         {logp, propz} = randomcondibd[ibds[[i]], ibds[[i + 1]], propibds[[i]]];
         If[ logp === None,
             logprob = None;
             propibds = None;
             Break[],
             logprob += logp;
             propibds[[i + 1]] = propz
         ], {i, Length[ibds] - 1}];
        {logprob, propibds}
    ]    


randomsegprop[1, snpls_, ibdls_, theta_, epsilon_, nsq_, nbp_,heat_] :=
    Module[ {propibd, ar2sub, propsub, propibdls, seglen,ar2,propsnpls,ar},
        propibd = IBDTransition[ibdls[[1, 2]], nsq, theta];
        If[ propibd == ibdls[[1, 2]],
            {1, 0,snpls, ibdls},
            propibdls = ibdls;
            propibdls[[1, 2 ;;]] = {propibd, LogPDFCRP[propibd, nsq, 0, theta]};
            {ar2sub, propsub} = randomsegcondibd[ibdls[[All, 2]], propibd];
            If[ ar2sub === None,
                {0, -1, snpls,ibdls},
                seglen = Length[propsub];
                propibdls[[;; seglen, 2]] = propsub;
                propibdls[[2 ;; seglen, -1]] = Log[IBDTransitionProb[propibdls[[;; seglen, 2]], nsq, theta]];
                ar2 = Log[IBDTransitionProb[{ibdls[[1, 2]], propibd}, nsq, theta][[1]]] - 
                  Log[IBDTransitionProb[{propibd, ibdls[[1, 2]]}, nsq, theta][[1]]];
                ar2 += ar2sub;
                {propsnpls, ar} = reviseloglike2[snpls, propibdls,1, seglen, epsilon, nbp];
                ar += Total[propibdls[[;;seglen, -1]]]-Total[ibdls[[;;seglen, -1]]];
                ar = Min[1, Exp[heat ar - ar2]];
                {ar, seglen, propsnpls,propibdls}
            ]
        ]
    ]
                 
    
randomsegprop[foc_/;foc>=2, snpls_, ibdls_, theta_, epsilon_, nsq_, nbp_,heat_] :=
    Module[ {ii,zibd, ar2sub, propsub, ar2sub2, propsub2, propibdls, seglen,ar2,ar,propsnpls},
        ii = foc-1;
        zibd = IBDTransition[ibdls[[ii, 2]], nsq, theta];
        {ar2sub, propsub} = randomsegcondibd[ibdls[[ii ;;, 2]], zibd];
        If[ ar2sub === None,
            {0, -1, snpls,ibdls},
            propibdls = Insert[ibdls, {-1, zibd, -1}, ii + 1];
            propibdls[[ii + 1 ;; ii + Length[propsub], 2]] = propsub;
            {ar2sub2, propsub2} = randomsegcondibd[propibdls[[ii + 1 ;;, 2]], ibdls[[ii, 2]]];
            (*Print[{ar2sub, ar2sub2}, Length[#] & /@ {propsub, propsub2}];*)
            If[ ar2sub2 === None,
                {0, -1, snpls,ibdls},
                propibdls = Delete[propibdls, ii + 1];
                propibdls[[ii ;; ii + Length[propsub2] - 1, 2]] = propsub2;
                seglen = Max[Length[#] & /@ {propsub, propsub2}];
                propibdls[[ii+1 ;; ii + seglen - 1, -1]] = Log[IBDTransitionProb[propibdls[[ii ;; ii + seglen - 1, 2]],nsq, theta]];
                ar2 = ar2sub + ar2sub2;
                {propsnpls, ar} = reviseloglike2[snpls, propibdls,ii+1, ii+seglen-1, epsilon, nbp];
                ar += Total[propibdls[[ii+1;; ii + seglen - 1, -1]]]-Total[ibdls[[ii+1;; ii + seglen - 1, -1]]];
                ar = Min[1, Exp[heat ar - ar2]];
                {ar, seglen, propsnpls,propibdls}
            ]
        ]
    ]    
    
            
End[]

EndPackage[]

