# mlcoalsim version 1.42 (20080318)
## usage: [input file] [output file]

### The application program mlcoalsim (multilocus coalescent simulations) is designed to:

	(i) Generate samples and calculate neutrality tests, and other statistics, under stationary model,several demographic models or strong positive selection by mean of coalescent theory. 
	
	(ii) Perform coalescent simulations with the mutational phase given:
	
		1. the population mutation rate θ (θ = 4Nμ, where N is the effective population size and μ is the mutational rate).
		2. a fixed number of mutations.
		3. a distribution of θ values. A prior uniform (bounded) and a gamma distributions are enabled.
		4. a fixed number of biallelic segregating sites taking into account the uncertainty of the population mutation rate (conditioning on biallelic segregating sites). A prior uniform (bounded) and a gamma distributions are enabled.
	
	(iii) Perform coalescent simulations with recombination given:

		1. the population recombination rate R (R = 4Nr, where r is the recombination rate).
		2. a distribution of r values. A prior uniform (bounded) and a gamma distributions are enabled.
		3. a fixed number of minimum recombination events (Rm) taking into account the uncer- tainty of the population recombination rate (fixing Rm). A prior uniform (bounded) and a gamma distributions are enabled.		
		4. a fixed number of minimum recombination events (Rm) and a fixed number of haplo- types, considering the uncertainty of the population recombination rate.

	(iv) Perform multilocus analyses. Linked loci and unlinked loci are enabled. Multilocus statistics for unlinked loci are the average and the variance for each statistic.
	
	(v) Include recurrent mutations (multiple hits) or not.
	
	(vi) Include heterogeneity in mutation rate across the length of the sequence. A gamma distri-
	bution is used. Also, a number of invariant positions can also be defined.
	
	(vii) Include heterogeneity in recombination rate across the length of the sequence. A gamma distribution is used. Hotspots or a constant value for all positions are possible.
	This program is based on a previous version of Hudson’s coalescent program ms (Hudson, 2002) and modified for the above purposes. The function to calculate minimum recombinant values is a modification of Wall’s code (Wall, 2000). The gamma function was partially obtained from Grassly, Adachi and Rambaut code (Grassly et al., 1997).

This program is distributed under the GNU GPL License.
