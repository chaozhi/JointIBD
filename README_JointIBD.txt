Steps for infering joint IBD from SNP data:
1) There are two folders "JointIBD" and "Mathematica Packages", and one text file "README_JointIBD.txt".


2) Run the mathematica notebook "JointIBD_MathematicaCode\JointIBD Inference\Infer_JointIBD.nb" by double clicking it.

3.1) In "Infer_JointIBD.nb", the default setup for running conditions:

	isrestart = False;
	isparallel = True;
	inphased = True;
	ascertainmodel = "M1";
	ngroupofchain = 2;
	nchainofgroup = 8;
	heatingAR = 0.5;
	datafile = "Data_S-NoLD.txt";
	outcompact = "JointIBD_OutCompact_" <> If[isphased, "Phased_", "Unphased_"] <> datafile
	outcold = "JointIBD_OutCold_" <> If[isphased, "Phased_", "Unphased_"] <> datafile
	outset = "JointIBD_OutSetting_" <> If[isphased, "Phased_", "Unphased_"] <> datafile
	outstate = "JointIBD_MCMCState_" <> If[isphased, "Phased_", "Unphased_"] <> datafile

	rjfreq = Ceiling[nbp/10^5];
	adpmax = 10000;
	adprate = 3/adpmax;
	itmin = 1;
	itmax = 1000001
	printfreq = 100;
	savefreq = 10;
	swapfreq = nchainofgroup;
	isprint = False;

3.2) 
        *) If isrestart=True, continue running by using the saved mcmc states from the previous running.
	*) If isparallel=True, parallel compuation for multiple MCMC chains
	*) If isphased=True, the SNP data is regarded as pahsed. 
	*) There are in total ngroupofchain*nchainofgroup mcmc chains. They are divided into "ngroupofchain" independent groups. 
	*) Within each group, parallel temperating algorithm is used, and the target accept ratio is heatingAR.
	*) The running conditions are saved in the outset file. And the MCMC state is saved in the outstate file, which may be used to restart/continue the mcmc chains. 
	*) The results are saved in the two files. Only cold chains (temperature=1) are saved in the e.g. "SMARTree_OutCold_M1_Standard_Data.txt". All the chains are saved un the outcompact where estimated local trees are excluded.


3.3) 
	*) rjfreq is the number of reversible jump sampling per iteration
	*) adpmax is the number of iniital iterations, during which the proposal distributions are adpatived.
	*) adprate is the rate of adaptation if isadp=True.
	*) itmax-itmin: is the max number of iteations.
	*) printfreq is the interval number of iterations at which some parameter values are printed on monitor.
	*) savefreq is the interval number of iteration at which results are saved.
	*) swapfreq is the number of swapping of mcmc chains with each group
	*) isprint=True, messages for the detailed steps are printed on monitor.
