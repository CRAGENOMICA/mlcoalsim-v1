perl ./mlcoalsim_mhmcmc.pl -infile ./ex01_input_likelihood.txt -numpar 2 -mlcout ex01allsims.out -burniter 500 -mcmciter 4500 -mcmcout ex01acceptedmcmc.out -parfile ex01par_mcmc.out -os mac < rangepar_mcmc_ex01.txt > ex01res_mcmc.out