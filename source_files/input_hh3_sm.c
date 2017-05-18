/************************************************************************************/
/* Copyright (C) 2006  Sebastian E. Ramos-Onsins                                    */
/*                                                                                  */ 
/* This program is free software; you can redistribute it and/or                    */
/* modify it under the terms of the GNU General Public License                      */
/* as published by the Free Software Foundation; either version 2                   */
/* of the License, or (at your option) any later version.                           */
/*                                                                                  */
/* This program is distributed in the hope that it will be useful,                  */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of                   */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the                    */
/* GNU General Public License for more details.                                     */
/*                                                                                  */
/* You should have received a copy of the GNU General Public License                */
/* along with this program; if not, write to the Free Software                      */
/* Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.  */
/************************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <errno.h>
#include "mlsp_sm.h"

int n_var = 138;
char var_file[138][20] =
{
    { "mhits"},		/* 1 if allowed multiple hits, 0 if not */

    { "n_iterations"},	/* 'howmany' iterations for each simulation */
    { "seed1"},		/* seed for ran1 */
    { "n_loci"},	/* number of loci */

    { "n_samples"},	/* samples in each locus */
    { "recombination"},	/* recombination in each locus(4Nor), f.ex: 2-10/2, 3 is in first locus 2 4 6 8 10 and 2nd is 3 */
    { "n_sites"},	/* number of nt in each locus */
    { "f"},		/* gene conversion ... */
    { "track_len"},	/* gene conversion ... */
    { "thetaw"},	/* theta for each loci, in case speciation, split or subdivision, only for the first group (4Nou) */
    { "mutations"},	/* (segsitesin in case no mhits), this is for the complete tree in each locus */
    { "npop"},		/* the number of pops. The same in each locus */
    { "ssize_pop"},	/* sample size of each population, comma and the next locus */
    { "mig_rate"},	/* vector with different rates of migration (en 4Nom)*/
    {"nlinked_loci"},	/* if linked, indicate the number of linked loci */
    {"pos_linked"},/* the initial position and the end of the fragment for all fragments, then "comma" and next linked loci*/
    {"print_matrixpol"},/* 0 not print, 1 DNA, 2 print Hudson modified output(mhits is different),3 dna excluding mhits*/
    {"mhratio_sv"},	/* in case mhits and theta. vector with the ratio vs transitions and transversions */
    {"seed2"},		/* Seed2 for generating random numbers */
    
    {"factor_pop"},	/* vector with the Ne proportion of each population */
    {"ran_factorpop"},	/* 0 no random, 1 yes, calculated every iteration */
    {"same_factorpop"},	/* 0 no constant, 1 same Ne for each population */
    {"neutral_tests"},	/* 0 neutral tests not calculated, 1 they are calcaluted */
    {"print_neuttest"},	/* 0 no print neut tests, 1 print average and variance from total loci, 2 all vales each loci */
    {"npop_sampled"},  	/* Number of populations sampled.*/

    { "theta"},/* REPEATED: theta for each loci, in case speciation, split or subdivision, only for the first group (4Nou) */
    { "sfix_allthetas"},/*Rejection algorithm and mhmcmc for Sfix and all thetas (0/1)*/
    { "mc_jump"},	/*separation among chosen values in the markov chain in Sfix_allthetas and rmfix*/
    /*Following Braverman et al. 1995*/
    {"ifselection"},	/*if selection 1, if not 0. Put 0 or 1 for each locus. Default is zero*/
    {"pop_size"},	/*effective population size: N. (Same for all loci)*/
    {"pop_sel"},	/*4Ns for each locus*/
    {"sel_nt"},		/*nt position of the selected mutation. 0 is the left value of the sequence. negative is allowed!!*/
    {"sinit"},		/*time in 4N generations the selection process ended. negative value means sel is not finished..*/
    /*limits for theta in case Sfix_all_theta*/
    {"range_thetant"},	/*Put limits to theta/nt (in case Sfix_all theta)*/
    {"thetant_min"},	/*theta/nt min*/
    {"thetant_max"},	/*theta/nt max*/
    
    {"dist_out"},	/* in mhits. Time of divergence in 4No gen */
    {"displ"},		/*in case sliding windows, displacement*/
    {"window"},		/*in case liding windows, size of the window*/

    { "nintn"},		/* change of the population size. Number of intervals */
    { "nrec"},		/* change of the population size. Initial size (en No). First must be 1.0 */
    { "npast"},		/* change of the population size. Final size (en No) */
    { "tpast"},		/* change of the population size. Duration (en 4No) */

    {"refugia"},	
    {"time_split"},
    {"time_scoal"},
    {"factor_anc"},
    {"freq_refugia"},
    
    {"tlimit"},		/* limit of the trees with recombination*/
	{"npoprefugia"},

	{"iflogistic"},/*in case use logistic growth, if not, instantaneous*/
	{"ts"},/*the first step in a growth process can start (for example) in the middle of the exponential phase*/
	
	{"ifgamma"},/*in case using a gamma distribution for theta values, if not, a uniform*/
	{"p_gamma"},
	{"alpha_gamma"},
	{"factorn_chr"},

	{"rmfix"}, /*uncertainty in recombination 0/1*/
	
	{"rm"}, /*values of Rm for each locus*/
	{"method_samp"}, /*methodology for sampling trees, in both Sfix_allthetas and rmfix: 1 rejection algorithm; 2 mcmc. mc_jump is useful also for both*/
	
	{"range_rnt"}, /*in case using a uniform distribution for all loci*/
	{"recnt_min"}, /*minimum value of recombination per nt*/
	{"recnt_max"}, /*maximum value of recombination per nt*/
	{"ifgammar"}, /*in case using a gamma distribution*/
	{"alpha_gammar"}, /*alpha value*/
	{"p_gammar"}, /*p,also named lambda*/

	{"no_rec_males"}, /*no recombination in males*/
	{"nhapl"}, /*values of ZA for each locus*/

	/*observed values*/
	{"td_obs"},
	{"fs_obs"},
	{"fdn_obs"},
	{"ffn_obs"},
	{"fd_obs"},
	{"ff_obs"},
	{"h_obs"},
	{"b_obs"},
	{"q_obs"},
	{"za_obs"},
	{"fst_obs"},
	{"kw_obs"},
	{"hw_obs"},
	{"r2_obs"},
	{"s_obs"},
	{"pi_w_obs"},
	{"pi_b_obs"},
	{"thetawatt_obs"},
	{"thetataj_obs"},
	{"thetafw_obs"},
	{"d_dmin_obs"},
	{"hnorm_obs"},
	{"maxhap_obs"},
	{"maxhap1_obs"},
	{"rm_obs"},

	{"linked_segsites"},
	{"linked_rm"},
	{"heter_theta_alphag"},
	{"heter_rec_alphag"},
	{"linked_nhapl"},

	{"correct_gamma"},
	{"correct_gammar"},

	{"invar_mut_sites"},

	{"likelihood_line"},
	
	{"td_err"},
	{"fs_err"},
	{"fdn_err"},
	{"ffn_err"},
	{"fd_err"},
	{"ff_err"},
	{"h_err"},
	{"b_err"},
	{"q_err"},
	{"za_err"},
	{"fst_err"},
	{"kw_err"},
	{"hw_err"},
	{"r2_err"},
	{"s_err"},
	{"pi_w_err"},
	{"pi_b_err"},
	{"thetawatt_err"},
	{"thetataj_err"},
	{"thetafw_err"},
	{"d_dmin_err"},
	{"hnorm_err"},
	{"maxhap_err"},
	{"maxhap1_err"},
	{"rm_err"},

	{"thetafl_obs"},
	{"thetafl_err"},
	{"thetal_obs"},
	{"thetal_err"},
	{"zenge_obs"},
	{"zenge_err"},

	{"ew_obs"},
	{"ew_err"},
	{"fstw_obs"},
	{"fstw_err"},
	{"pwh_obs"},
	{"pwh_err"},
};

void input_data( FILE *file_input,struct var **data)
{
    char *f;
    int c;
    char name_var[20];
    char number[20];
    int x,y,z;
    int nam,numb_1,numb_2/*,val*/;
    double **numb_c;
    /*double end,interval,st;*/
	int maxn=19;
    
    void init_seed1(long);
    double ran1(void);
	
	int totalnloci;
	
	memset(name_var,0,20*sizeof(char));

    if(!(*data = (struct var *)calloc(1,sizeof(struct var)))) perror("calloc error.0");;
    
    if(!((*data)->nsam = (int *)calloc((unsigned)20,sizeof(int)))) perror("calloc error.1");
    if(!((*data)->r = (double *)calloc((unsigned)20,sizeof(double )))) perror("calloc error.2");
    if(!((*data)->nsites = (long int *)calloc((unsigned)20,sizeof(long int)))) perror("calloc error.4");    
    if(!((*data)->loci_linked = (long int **)calloc((unsigned)20,sizeof(long int *)))) perror("calloc error.4b");
    for(x=0;x<20;x++) 
        if(!((*data)->loci_linked[x] = (long int *)calloc((unsigned)20,sizeof(long int))))
            perror("calloc error.4c");
    if(!((*data)->f = (double *)calloc((unsigned)20,sizeof(double )))) perror("calloc error.5");
    if(!((*data)->track_len = (double *)calloc((unsigned)20,sizeof(double )))) perror("calloc error.7");
    if(!((*data)->theta_1 = (double *)calloc((unsigned)20,sizeof(double )))) perror("calloc error.9");    

    if(!((*data)->mutations = (long int *)calloc((unsigned)20,sizeof(long int)))) perror("calloc error.10");
    if(!((*data)->ssize_pop = (int **)calloc((unsigned)20,sizeof(int *)))) perror("calloc error.11");
    for(x=0;x<20;x++) if(!((*data)->ssize_pop[x] = (int *)calloc((unsigned)20,sizeof(int)))) perror("calloc error.11b");
    if(!((*data)->factor_pop = (double *)calloc((unsigned)20,sizeof(double)))) perror("calloc error.34");
    if(!((*data)->ratio_sv = (double *)calloc((unsigned)20,sizeof(double)))) perror("calloc error.34");
    if(!((*data)->npop_sampled = (int *)calloc((unsigned)20,sizeof(int)))) perror("calloc error.34");
    if(!((*data)->ifselection = (int *)calloc((unsigned)20,sizeof(int)))) perror("calloc error.34");
    if(!((*data)->pop_sel = (double *)calloc((unsigned)20,sizeof(double)))) perror("calloc error.34");
    if(!((*data)->sel_nt = (long int *)calloc((unsigned)20,sizeof(long int)))) perror("calloc error.34");
    if(!((*data)->sinit = (double *)calloc((unsigned)20,sizeof(double)))) perror("calloc error.34");
    if(!((*data)->nrec = (double *)calloc((unsigned)20,sizeof(double)))) perror("calloc error.14");
    if(!((*data)->npast = (double *)calloc((unsigned)20,sizeof(double)))) perror("calloc error.15");
    if(!((*data)->tpast = (double *)calloc((unsigned)20,sizeof(double)))) perror("calloc error.16");
    if(!((*data)->freq = (double *)calloc((unsigned)20,sizeof(double)))) perror("calloc error.16");
    
	if(!((*data)->p_gamma = (double *)calloc((unsigned)20,sizeof(double)))) perror("calloc error.56");
    if(!((*data)->alpha_gamma = (double *)calloc((unsigned)20,sizeof(double)))) perror("calloc 57");
    if(!((*data)->correct_gamma = (double *)calloc((unsigned)20,sizeof(double)))) perror("calloc 57");
    if(!((*data)->factor_chrn = (double *)calloc((unsigned)20,sizeof(double)))) perror("calloc 58");
    
	if(!((*data)->p_gammar = (double *)calloc((unsigned)20,sizeof(double)))) perror("calloc error.56");
	if(!((*data)->alpha_gammar = (double *)calloc((unsigned)20,sizeof(double)))) perror("calloc error.56");
    if(!((*data)->correct_gammar = (double *)calloc((unsigned)20,sizeof(double)))) perror("calloc 57");
	if(!((*data)->Rm = (int *)calloc((unsigned)20,sizeof(int)))) perror("calloc error.56");
	if(!((*data)->nhapl = (int *)calloc((unsigned)20,sizeof(int)))) perror("calloc error.56");

	if(!((*data)->linked_segsites = (int *)calloc((unsigned)20,sizeof(int)))) perror("calloc error.89");
	if(!((*data)->linked_rm = (int *)calloc((unsigned)20,sizeof(int)))) perror("calloc error.90");
	if(!((*data)->linked_nhapl = (int *)calloc((unsigned)20,sizeof(int)))) perror("calloc error.90");

	if(!((*data)->heter_theta_alphag = (double *)calloc((unsigned)20,sizeof(double)))) perror("calloc error.91");
	if(!((*data)->heter_rm_alphag = (double *)calloc((unsigned)20,sizeof(double)))) perror("calloc error.88");
	if(!((*data)->invariable_mut_sites = (long int *)calloc((unsigned)20,sizeof(long int)))) perror("calloc error.9");

	if(!((*data)->thetant_min = (double *)calloc((unsigned)20,sizeof(double)))) perror("calloc error.88");
	if(!((*data)->thetant_max = (double *)calloc((unsigned)20,sizeof(double)))) perror("calloc error.88");

	if(!((*data)->recnt_min = (double *)calloc((unsigned)20,sizeof(double)))) perror("calloc error.88");
	if(!((*data)->recnt_max = (double *)calloc((unsigned)20,sizeof(double)))) perror("calloc error.88");

    (*data)->pr_matrix = 1;
    (*data)->mhits = 0;
    (*data)->n_iter = 200;
    (*data)->n_loci = 1;
    (*data)->linked = 0;
    (*data)->despl = 0;
    (*data)->window = 0;
    (*data)->npop = 1;
    (*data)->mig_rate = 0.0;
    
    (*data)->seed1 = 12345678;
    (*data)->seed2 = 546373;

    (*data)->ran_factorpop = 0;
    (*data)->same_factorpop = 1;
    (*data)->neutral_tests = 0;
    (*data)->print_neuttest = 0;

    (*data)->sfix_allthetas = 0;
    (*data)->mc_jump = 1;
    (*data)->pop_size = 1E+06;
    
    (*data)->range_thetant = 0;
    
    (*data)->T_out = 0.0;
    (*data)->nintn = 0;
       
    (*data)->split_pop = 0;
    (*data)->time_split = 0.;
    (*data)->time_scoal = 0.;
    (*data)->factor_anc = 1.;
    
    (*data)->tlimit = 1000.;
    
	(*data)->iflogistic = 1;
    (*data)->ts = 0.;

    (*data)->ifgamma = 0;

    (*data)->rmfix = 0;
    (*data)->range_rnt = 0;
    (*data)->ifgammar = 0;
    (*data)->no_rec_males = 0;

	/*observed values*/
	for(x=0;x<NOBS_STATISTICS;x++) {
		if(!((*data)->obs_statistics[x] = (double *)calloc((unsigned)20,sizeof(double)))) perror("calloc error.561");
		(*data)->obs_statistics[x][0] = (double)0;/*redundant*/
		(*data)->likelihood_error[x] = (double)1e-6;
	}

	(*data)->likelihood_line = 0;

    if(!(f = (char *)malloc(BUFSIZ))) perror("malloc error.1");
    setbuf(file_input,f);
    
    x = 0;
    c = 1;
    
    if(!(numb_c = (double **)malloc(20*sizeof(double *)))) perror("malloc error.2");
    for(x=0;x<20;x++) if(!(numb_c[x] = (double *)malloc(20*sizeof(double)))) perror("realloc error.1");

    while(c > 0) {
        
        /* look for the name of the variable */
        nam = 0;
        while(!nam) {
            if(c == 34) { /* 34 is '"' */
                c = fgetc(file_input);
                while((c > 0) && (c != 34))
                    c = fgetc(file_input);
            }
            if(c <= 0) break;
            if(c != 34) {/*Modificat el 8.5.03 */
                if((c < 91 && c > 64) || (c < 123 && c > 96) || c == 95) {/* if letters or '_' */
                    x=0;
                    while(c && ((c < 91 && c > 64) || (c < 123 && c > 96) || c == 95 || (c > 47 && c < 58))) {
                        name_var[x] = c;	/* assign c to name_var */
                        x++;
                        c = fgetc(file_input);
                        if(c <= 0) break;
                    }
                    name_var[x] = '\0';
                    for(y=0;y<x;y++) if(name_var[y] < 91 && name_var[y] > 64) name_var[y] += 32;/* do lowercase */
                    nam = 1;		/* we have new name_var */
                    if(c <= 0) break;
                }
                else {
                    if(c < 58 && c > 47) {
                        puts("error in input. Numbers with no variable"); 
                        exit(1);
                    }
                    c = fgetc(file_input);
                }
                if(c <= 0) break;
                
                if(name_var[0] == 101 && x == 1) {	/* in case we have E+- */
                    puts("error in input. Numbers with no variable"); 
                    exit(1);
                }
            }
            else c = fgetc(file_input);/*Modificat el 8.5.03 */
        }
        
        /* look for the numbers of the variable */                
        numb_1 = x = 0;
        numb_2 = 1;
        /*if(!(numb_c[0] = realloc(numb_c[0],2*sizeof(double)))) perror("realloc error.1");*/
        numb_c[numb_1][0] = 0;
    
        while(nam) {
            while((c < 58 && c > 42) || c == 69 || c == 101) {
                if((c < 58 && c > 47) || c == 69 || c == 101 || c == 43 || c == 46) { 	/* numbers, E, dec, + */
                    number[x] = c;
                    x++;
                }
                if(c == '-') {		/* for a minus, but also for first-end */
                    if(x && number[x-1] < 58 && number[x-1] > 47) {
                        number[x] = '\0';
                        numb_c[numb_1][numb_2] = atof(number);
                        numb_c[numb_1][0]++;
                        numb_2++;
                        if(numb_2 >= 20) if(!(numb_c[numb_1] = realloc(numb_c[numb_1],(numb_2+1)*sizeof(double))))
                            perror("realloc error.2");
                        x = 0;
                    }
                    else {
                            number[x] = c;
                            x++;
                    }                                    
                }
				/* intervals *//*
                if(c == '/') {			
                    if(x) {
                        number[x] = '\0';
                        numb_c[numb_1][numb_2] = atof(number);
                        numb_c[numb_1][0]++;
                    }
                    if(numb_2 != 2) {
                        puts("Error in input data. Intervals usage using '/' is: 'first'-'end'/'interval'");
                        exit(1);
                    }
                    numb_2++;
                    if(numb_2 >= 20) 
                        if(!(numb_c[numb_1] = realloc(numb_c[numb_1],(numb_2+1)*sizeof(double)))) 
                            perror("realloc error.3");
                    c=fgetc(file_input);                                     
                    while(!((c < 58 && c > 42) || c == 69 || c == 101)) c=fgetc(file_input);
                    x = 0;
                    while((c < 58 && c > 47) || c == 69 || c == 101 || c == 43 || c == 45 || c == 46) {
                        number[x] = c;
                        x++;
                        c = fgetc(file_input);
                    }
                    number[x] = '\0';
                    end = numb_c[numb_1][2];
                    interval = atof(number);
                    st = numb_c[numb_1][1];
                    numb_2 = 1;
                    val = (int)floor((double)((end-st)/interval)) + 2;
                    if(val >= 20) 
                        if(!(numb_c[numb_1] = realloc(numb_c[numb_1],val*(sizeof(double)))))
                            perror("realloc error.4");
                    numb_c[numb_1][0] = 0;
                    
                    while(st <= end) {
                        numb_c[numb_1][numb_2] = st;
                        numb_c[numb_1][0] ++;
                        st += interval;
                        numb_2++;
                    }
                    x = 0;
                }
				*/
                if(c == ',') {			/* do a new vector */
                    if(x) {
                        number[x] = '\0';
                        numb_c[numb_1][numb_2] = atof(number);
                        numb_c[numb_1][0] ++;
                    }
                    numb_2 = 1;
                    numb_1++;
                    if(numb_1 >= 20) {
                        if(!(numb_c = realloc(numb_c,(numb_1+1)*sizeof(double *)))) perror("realloc error.5");
                        if(!(numb_c[numb_1] = (double *) malloc(20 * sizeof(double)))) perror("realloc error.6");
                    }
                    numb_c[numb_1][0] = 0;
                    x = 0;
                }
                c = fgetc(file_input);
                if(c == 34) {
                    c = fgetc(file_input);
                    while((c > 0) && (c != 34))
                        c = fgetc(file_input);
                }
                if(!((c < 58 && c > 42) || c == 69 || c == 101)) {
                    if((c < 91 && c > 64) || (c < 123 && c > 96) || c == 95) {	/* if letters, stop */
                        nam = 1;
                        break;
                    }
                    else {
                        if(x) {
                            number[x] = '\0';
                            numb_c[numb_1][numb_2] = atof(number);
                            numb_c[numb_1][0]++;
                            numb_2++;
                            if(numb_2 >= 20) 
                                if(!(numb_c[numb_1] = realloc(numb_c[numb_1],(numb_2+1)*sizeof(double)))) perror("realloc error.7");
                            x = 0;
                        }
                    }
                }
            }
            if((c < 91 && c > 64) || (c < 123 && c > 96) || c == 95) {	/* if letters, stop */
                nam = 1;
                break;
            }
            if(c == 34) {
                c = fgetc(file_input);
                while((c > 0) && (c != 34))
                    c = fgetc(file_input);
            }
            if(c <= 0) break;
            c = fgetc(file_input);
        }
		if(maxn < numb_1) maxn = numb_1;
        
        /* do assignments to the struct variables */
        for(x=0;x<n_var;x++)
            if(strcmp(name_var,var_file[x]) == 0)
                break;
        
        switch(x) {
            case 0:
                (*data)->mhits = (int)numb_c[0][1];
                break;
            case 1:
                (*data)->n_iter = (long int)numb_c[0][1];
                break;
            case 2:
                (*data)->seed1 = (long int)numb_c[0][1];
                break;
            case 3:
                (*data)->n_loci = (int)numb_c[0][1];
                break;
            case 4:
                if(numb_c[0][0] >= 20) 
                    if(!((*data)->nsam = (int *)realloc((*data)->nsam,(unsigned)(numb_c[0][0]+1)*sizeof(int))))
                        perror("realloc error.9");
                for(z=0;z<numb_c[0][0]+1;z++) (*data)->nsam[z] = (int)numb_c[0][z];
                break;
            case 5:
                if(numb_c[0][0] >= 20) 
                    if(!((*data)->r = (double *)realloc((*data)->r,(unsigned)(numb_c[0][0]+1)*sizeof(double))))
                        perror("realloc error.9");
                for(z=0;z<numb_c[0][0]+1;z++)
					(*data)->r[z] = numb_c[0][z];
                break;
            case 6:
                if(numb_c[0][0] >= 20) 
                    if(!((*data)->nsites = (long int *)realloc((*data)->nsites,
                                           (unsigned)(numb_c[0][0]+1)*sizeof(long int))))
                    perror("realloc error.13");
                for(z=0;z<numb_c[0][0]+1;z++) (*data)->nsites[z] = (long int)numb_c[0][z];
                break;
            case 7:
                if(numb_c[0][0] >= 20) 
                    if(!((*data)->f = (double *)realloc((*data)->f,(unsigned)(numb_c[0][0]+1)*sizeof(double))))
                        perror("realloc error.9");
                for(z=0;z<numb_c[0][0]+1;z++) (*data)->f[z] = numb_c[0][z];
                break;
            case 8:
                if(numb_c[0][0] >= 20) 
                    if(!((*data)->track_len = (double *)realloc((*data)->track_len,
                                        (unsigned)(numb_c[0][0]+1)*sizeof(double))))
                        perror("realloc error.9");
                for(z=0;z<numb_c[0][0]+1;z++) (*data)->track_len[z] = numb_c[0][z];
                break;
            case 9:
                if(numb_c[0][0] >= 20) 
                    if(!((*data)->theta_1 = (double *)realloc((*data)->theta_1,
                    (unsigned)(numb_c[0][0]+1)*sizeof(double)))) 
                        perror("realloc error.33");
                for(z=0;z<numb_c[0][0]+1;z++) (*data)->theta_1[z] = numb_c[0][z];
                break;
            case 10:
                if(numb_c[0][0] >= 20) if(!((*data)->mutations = (long int *)realloc((*data)->mutations,
                                                (unsigned)(numb_c[0][0]+1)*sizeof(long int)))) 
                    perror("realloc error.21");
                for(z=0;z<numb_c[0][0]+1;z++) (*data)->mutations[z] = (long int)numb_c[0][z];
                break;
            case 11:
                (*data)->npop = (long int)numb_c[0][1];
                break;
            case 12:
                if(numb_1+1 != (*data)->n_loci && numb_1 != 0) {/*numb_1 is the number of dimensions. i.e. loci*/
                    perror("Error: nloci must be first defined or ssize_pop is different from nloci. ");
                    exit(1);
                }
                if(numb_1 >= 20) {
                    if(!((*data)->ssize_pop = (int **)realloc((*data)->ssize_pop,(unsigned)(numb_1+1)*sizeof(int *))))
                        perror("realloc error.22");
                }
                for(z=0;z<numb_1+1;z++) {
                    if(z>=20) {
                        if(!((*data)->ssize_pop[z] = (int *)malloc((unsigned)(numb_c[z][0]+1)*sizeof(int))))
                            perror("realloc error.11");
                    }
                    else {
                        if(numb_c[z][0] >= 20) 
                            if(!((*data)->ssize_pop[z]=(int *)realloc((*data)->ssize_pop[z],
                                                        (unsigned)(numb_c[z][0]+1)*sizeof(int))))
                                perror("realloc error.22b");
                    }
                    for(y=0;y<numb_c[z][0]+1;y++) (*data)->ssize_pop[z][y] = (int)numb_c[z][y];
                }
                break;
            case 13:
                (*data)->mig_rate = numb_c[0][1];
                break;
            case 14:
                (*data)->linked = (int)numb_c[0][1];
                break;
            case 15:
                if(numb_1 >= 20) {
                    if(!((*data)->loci_linked = (long int **)realloc((*data)->loci_linked,
                    (unsigned)(numb_1+1)*sizeof(long int *))))
                        perror("realloc error.38");
                }
                for(z=0;z<numb_1+1;z++) {
                    if(z>=20) {
                        if(!((*data)->loci_linked[z] = 
                        (long int *)malloc((unsigned)(numb_c[z][0]+1)*sizeof(long int))))
                            perror("realloc error.11");
                    }
                    else {
                        if(numb_c[z][0] >= 20) 
                        if(!((*data)->loci_linked[z] = (long int *)realloc((*data)->loci_linked[z],
                        (unsigned)(numb_c[z][0]+1)* sizeof(long int))))
                            perror("realloc error.39");
                    }
                    for(y=0;y<numb_c[z][0]+1;y++) (*data)->loci_linked[z][y] = (long int)numb_c[z][y];
                }
                break;
            case 16:
                (*data)->pr_matrix = (int)numb_c[0][1];
                break;
            case 17:
                if(numb_c[0][0] >= 20) 
                    if(!((*data)->ratio_sv = (double *)realloc((*data)->ratio_sv,
                    (unsigned)(numb_c[0][0]+1)*sizeof(double)))) 
                        perror("realloc error.33");
                for(z=0;z<numb_c[0][0]+1;z++) (*data)->ratio_sv[z] = numb_c[0][z];
                break;
            case 18:
                (*data)->seed2 = (long int)numb_c[0][1];
                break;
            case 19:
                if(numb_c[0][0] >= 20) 
                    if(!((*data)->factor_pop = (double *)realloc((*data)->factor_pop,
                    (unsigned)(numb_c[0][0]+1)*sizeof(double)))) 
                        perror("realloc error.58");
                for(z=0;z<numb_c[0][0]+1;z++) (*data)->factor_pop[z] = (double)numb_c[0][z];                
                break;
            case 20:
                (*data)->ran_factorpop = (int)numb_c[0][1];
                break;
            case 21:
                (*data)->same_factorpop = (int)numb_c[0][1];
                break;
            case 22:
                (*data)->neutral_tests = (int)numb_c[0][1];
                break;
            case 23:
                (*data)->print_neuttest = (int)numb_c[0][1];
                break;            
            case 24:
                if(numb_c[0][0] >= 20) {
                    if(!((*data)->npop_sampled = (int *)realloc((*data)->npop_sampled,
                    (unsigned)(numb_c[0][0]+1)*sizeof(int)))) 
                        perror("realloc error.58");
                }
                for(z=0;z<numb_c[0][0]+1;z++) (*data)->npop_sampled[z] = (int)numb_c[0][z];                
                break;
            case 25: /*Exactly equal than case 9*/
                if(numb_c[0][0] >= 20) 
                    if(!((*data)->theta_1 = (double *)realloc((*data)->theta_1,
                    (unsigned)(numb_c[0][0]+1)*sizeof(double)))) 
                        perror("realloc error.33");
                for(z=0;z<numb_c[0][0]+1;z++) (*data)->theta_1[z] = numb_c[0][z];
                break;
            case 26:
                (*data)->sfix_allthetas = (int)numb_c[0][1];
                break;
            case 27:
                (*data)->mc_jump = (int)numb_c[0][1];
                break;
            case 28:
                if(numb_c[0][0] >= 20) 
                    if(!((*data)->ifselection = (int *)realloc((*data)->ifselection,(unsigned)(numb_c[0][0]+1)*sizeof(int))))
                        perror("realloc error.9");
                for(z=0;z<numb_c[0][0]+1;z++) (*data)->ifselection[z] = (int)numb_c[0][z];
                break;
            case 29:
                (*data)->pop_size = (double)numb_c[0][1];
                break;
            case 30:
                if(numb_c[0][0] >= 20) 
                    if(!((*data)->pop_sel = (double *)realloc((*data)->pop_sel,(unsigned)(numb_c[0][0]+1)*sizeof(double))))
                        perror("realloc error.9");
                for(z=0;z<numb_c[0][0]+1;z++) (*data)->pop_sel[z] = (double)numb_c[0][z];
                break;
            case 31:
                if(numb_c[0][0] >= 20) 
                    if(!((*data)->sel_nt = (long int *)realloc((*data)->sel_nt,(unsigned)(numb_c[0][0]+1)*sizeof(long int))))
                        perror("realloc error.9");
                for(z=0;z<numb_c[0][0]+1;z++) (*data)->sel_nt[z] = (long int)numb_c[0][z];
                break;
            case 32:
                if(numb_c[0][0] >= 20) 
                    if(!((*data)->sinit = (double *)realloc((*data)->sinit,(unsigned)(numb_c[0][0]+1)*sizeof(double))))
                        perror("realloc error.9");
                for(z=0;z<numb_c[0][0]+1;z++) (*data)->sinit[z] = (double)numb_c[0][z];
                break;
            case 33:
                (*data)->range_thetant = (int)numb_c[0][1];
                break;
            case 34:
                if(numb_c[0][0] >= 20) 
                    if(!((*data)->thetant_min = (double *)realloc((*data)->thetant_min,(unsigned)(numb_c[0][0]+1)*sizeof(double))))
                        perror("realloc error.69");
                for(z=0;z<numb_c[0][0]+1;z++) (*data)->thetant_min[z] = (double)numb_c[0][z];
                break;
            case 35:
                if(numb_c[0][0] >= 20) 
                    if(!((*data)->thetant_max = (double *)realloc((*data)->thetant_max,(unsigned)(numb_c[0][0]+1)*sizeof(double))))
                        perror("realloc error.69");
                for(z=0;z<numb_c[0][0]+1;z++) (*data)->thetant_max[z] = (double)numb_c[0][z];
                break;
            case 36:
                (*data)->T_out = numb_c[0][1];
                break;
            case 37:
                (*data)->despl = (long int)numb_c[0][1];
                break;
            case 38:
                (*data)->window = (long int)numb_c[0][1];
                break;
            case 39:
                (*data)->nintn = (int)numb_c[0][1];
                break;
            case 40:
                if(numb_c[0][0] >= 20) 
                    if(!((*data)->nrec = (double *)realloc((*data)->nrec,(unsigned)(numb_c[0][0]+1)*sizeof(double)))) 
                        perror("realloc error.26");
                for(z=0;z<numb_c[0][0]+1;z++) (*data)->nrec[z] = numb_c[0][z];
                break;
            case 41:
                if(numb_c[0][0] >= 20) 
                    if(!((*data)->npast = (double *)realloc((*data)->npast,(unsigned)(numb_c[0][0]+1)*sizeof(double)))) 
                        perror("realloc error.27");
                for(z=0;z<numb_c[0][0]+1;z++) (*data)->npast[z] = numb_c[0][z];
                break;
            case 42:
                if(numb_c[0][0] >= 20)
                    if(!((*data)->tpast = (double *)realloc((*data)->tpast,(unsigned)(numb_c[0][0]+1)*sizeof(double)))) 
                        perror("realloc error.28");
                for(z=0;z<numb_c[0][0]+1;z++) (*data)->tpast[z] = numb_c[0][z];
                break;
            case 43:
                (*data)->split_pop = (int)numb_c[0][1];
                break;
            case 44:
                (*data)->time_split = (double)numb_c[0][1];
                break;
            case 45:
                (*data)->time_scoal = (double)numb_c[0][1];
                break;
            case 46:
                (*data)->factor_anc = (double)numb_c[0][1];
                break;
            case 47:
                if(numb_c[0][0] >= 20) 
                    if(!((*data)->freq = (double *)realloc((*data)->freq,(unsigned)(numb_c[0][0]+1)*sizeof(double)))) 
                        perror("realloc error.26");
                for(z=0;z<numb_c[0][0]+1;z++) (*data)->freq[z] = numb_c[0][z];
                break;
            case 48:
                (*data)->tlimit = numb_c[0][1];
                break;
            case 49:
                (*data)->npoprefugia = (int)numb_c[0][1];
                break;
            case 50:
                (*data)->iflogistic = (int)numb_c[0][1];
                break;
            case 51:
                (*data)->ts = numb_c[0][1];
                break;
            case 52:
                (*data)->ifgamma = (int)numb_c[0][1];
                break;
            case 53:
                if(numb_c[0][0] >= 20) 
                    if(!((*data)->p_gamma = (double *)realloc((*data)->p_gamma,(unsigned)(numb_c[0][0]+1)*sizeof(double)))) 
                        perror("realloc error.53");
                for(z=0;z<numb_c[0][0]+1;z++) (*data)->p_gamma[z] = (double)numb_c[0][z];
                break;
            case 54:
                if(numb_c[0][0] >= 20) 
                    if(!((*data)->alpha_gamma = (double *)realloc((*data)->alpha_gamma,(unsigned)(numb_c[0][0]+1)*sizeof(double)))) 
                        perror("realloc error.54");
                for(z=0;z<numb_c[0][0]+1;z++) (*data)->alpha_gamma[z] = (double)numb_c[0][z];
                break;
            case 55:
                if(numb_c[0][0] >= 20) 
                    if(!((*data)->factor_chrn = (double *)realloc((*data)->factor_chrn,(unsigned)(numb_c[0][0]+1)*sizeof(double)))) 
                        perror("realloc error.55");
                for(z=0;z<numb_c[0][0]+1;z++) (*data)->factor_chrn[z] = (double)numb_c[0][z];
                break;
            case 56:
                (*data)->rmfix = (int)numb_c[0][1];
                break;
            case 57:
                if(numb_c[0][0] >= 20) 
                    if(!((*data)->Rm = (int *)realloc((*data)->Rm,(unsigned)(numb_c[0][0]+1)*sizeof(int)))) 
                        perror("realloc error.57");
                for(z=0;z<numb_c[0][0]+1;z++) (*data)->Rm[z] = (int)numb_c[0][z];
                break;
            case 58:
                (*data)->method_samp = (int)numb_c[0][1];
                break;
            case 59:
                (*data)->range_rnt = (int)numb_c[0][1];
                break;
            case 60:
                if(numb_c[0][0] >= 20) 
                    if(!((*data)->recnt_min = (double *)realloc((*data)->recnt_min,(unsigned)(numb_c[0][0]+1)*sizeof(double))))
                        perror("realloc error.69");
                for(z=0;z<numb_c[0][0]+1;z++) (*data)->recnt_min[z] = (double)numb_c[0][z];
                break;
            case 61:
                if(numb_c[0][0] >= 20) 
                    if(!((*data)->recnt_max = (double *)realloc((*data)->recnt_max,(unsigned)(numb_c[0][0]+1)*sizeof(double))))
                        perror("realloc error.69");
                for(z=0;z<numb_c[0][0]+1;z++) (*data)->recnt_max[z] = (double)numb_c[0][z];
                break;
            case 62:
                (*data)->ifgammar = (int)numb_c[0][1];
                break;
            case 63:
                if(numb_c[0][0] >= 20) 
                    if(!((*data)->alpha_gammar = (double *)realloc((*data)->alpha_gammar,(unsigned)(numb_c[0][0]+1)*sizeof(double)))) 
                        perror("realloc error.62");
                for(z=0;z<numb_c[0][0]+1;z++) (*data)->alpha_gammar[z] = (double)numb_c[0][z];
                break;
            case 64:
                if(numb_c[0][0] >= 20) 
                    if(!((*data)->p_gammar = (double *)realloc((*data)->p_gammar,(unsigned)(numb_c[0][0]+1)*sizeof(double)))) 
                        perror("realloc error.63");
                for(z=0;z<numb_c[0][0]+1;z++) (*data)->p_gammar[z] = (double)numb_c[0][z];
                break;
            case 65:
                (*data)->no_rec_males = (int)numb_c[0][1];
                break;
            case 66:
                if(numb_c[0][0] >= 20) 
                    if(!((*data)->nhapl = (int *)realloc((*data)->nhapl,(unsigned)(numb_c[0][0]+1)*sizeof(int)))) 
                        perror("realloc error.57");
                for(z=0;z<numb_c[0][0]+1;z++) (*data)->nhapl[z] = (int)numb_c[0][z];
                break;
			/*observed values*/
            case 67:
                if(numb_c[0][0] >= 20) 
                    if(!((*data)->obs_statistics[0] = (double *)realloc((*data)->obs_statistics[0],(unsigned)(numb_c[0][0]+1)*sizeof(double))))
                        perror("realloc error.9");
                for(z=0;z<numb_c[0][0]+1;z++) (*data)->obs_statistics[0][z] = (double)numb_c[0][z];
                break;
            case 68:
                if(numb_c[0][0] >= 20) 
                    if(!((*data)->obs_statistics[1] = (double *)realloc((*data)->obs_statistics[1],(unsigned)(numb_c[0][0]+1)*sizeof(double))))
                        perror("realloc error.9");
                for(z=0;z<numb_c[0][0]+1;z++) (*data)->obs_statistics[1][z] = (double)numb_c[0][z];
                break;
            case 69:
                if(numb_c[0][0] >= 20) 
                    if(!((*data)->obs_statistics[2] = (double *)realloc((*data)->obs_statistics[2],(unsigned)(numb_c[0][0]+1)*sizeof(double))))
                        perror("realloc error.9");
                for(z=0;z<numb_c[0][0]+1;z++) (*data)->obs_statistics[2][z] = (double)numb_c[0][z];
                break;
            case 70:
                if(numb_c[0][0] >= 20) 
                    if(!((*data)->obs_statistics[3] = (double *)realloc((*data)->obs_statistics[3],(unsigned)(numb_c[0][0]+1)*sizeof(double))))
                        perror("realloc error.9");
                for(z=0;z<numb_c[0][0]+1;z++) (*data)->obs_statistics[3][z] = (double)numb_c[0][z];
                break;
            case 71:
                if(numb_c[0][0] >= 20) 
                    if(!((*data)->obs_statistics[4] = (double *)realloc((*data)->obs_statistics[4],(unsigned)(numb_c[0][0]+1)*sizeof(double))))
                        perror("realloc error.9");
                for(z=0;z<numb_c[0][0]+1;z++) (*data)->obs_statistics[4][z] = (double)numb_c[0][z];
                break;
            case 72:
                if(numb_c[0][0] >= 20) 
                    if(!((*data)->obs_statistics[5] = (double *)realloc((*data)->obs_statistics[5],(unsigned)(numb_c[0][0]+1)*sizeof(double))))
                        perror("realloc error.9");
                for(z=0;z<numb_c[0][0]+1;z++) (*data)->obs_statistics[5][z] = (double)numb_c[0][z];
                break;
            case 73:
                if(numb_c[0][0] >= 20) 
                    if(!((*data)->obs_statistics[6] = (double *)realloc((*data)->obs_statistics[6],(unsigned)(numb_c[0][0]+1)*sizeof(double))))
                        perror("realloc error.9");
                for(z=0;z<numb_c[0][0]+1;z++) (*data)->obs_statistics[6][z] = (double)numb_c[0][z];
                break;
            case 74:
                if(numb_c[0][0] >= 20) 
                    if(!((*data)->obs_statistics[7] = (double *)realloc((*data)->obs_statistics[7],(unsigned)(numb_c[0][0]+1)*sizeof(double))))
                        perror("realloc error.9");
                for(z=0;z<numb_c[0][0]+1;z++) (*data)->obs_statistics[7][z] = (double)numb_c[0][z];
                break;
            case 75:
                if(numb_c[0][0] >= 20) 
                    if(!((*data)->obs_statistics[8] = (double *)realloc((*data)->obs_statistics[8],(unsigned)(numb_c[0][0]+1)*sizeof(double))))
                        perror("realloc error.9");
                for(z=0;z<numb_c[0][0]+1;z++) (*data)->obs_statistics[8][z] = (double)numb_c[0][z];
                break;
            case 76:
                if(numb_c[0][0] >= 20) 
                    if(!((*data)->obs_statistics[9] = (double *)realloc((*data)->obs_statistics[9],(unsigned)(numb_c[0][0]+1)*sizeof(double))))
                        perror("realloc error.9");
                for(z=0;z<numb_c[0][0]+1;z++) (*data)->obs_statistics[9][z] = (double)numb_c[0][z];
                break;
            case 77:
                if(numb_c[0][0] >= 20) 
                    if(!((*data)->obs_statistics[10] = (double *)realloc((*data)->obs_statistics[10],(unsigned)(numb_c[0][0]+1)*sizeof(double))))
                        perror("realloc error.9");
                for(z=0;z<numb_c[0][0]+1;z++) (*data)->obs_statistics[10][z] = (double)numb_c[0][z];
                break;
            case 78:
                if(numb_c[0][0] >= 20) 
                    if(!((*data)->obs_statistics[11] = (double *)realloc((*data)->obs_statistics[11],(unsigned)(numb_c[0][0]+1)*sizeof(double))))
                        perror("realloc error.9");
                for(z=0;z<numb_c[0][0]+1;z++) (*data)->obs_statistics[11][z] = (double)numb_c[0][z];
                break;
            case 79:
                if(numb_c[0][0] >= 20) 
                    if(!((*data)->obs_statistics[12] = (double *)realloc((*data)->obs_statistics[12],(unsigned)(numb_c[0][0]+1)*sizeof(double))))
                        perror("realloc error.9");
                for(z=0;z<numb_c[0][0]+1;z++) (*data)->obs_statistics[12][z] = (double)numb_c[0][z];
                break;
            case 80:
                if(numb_c[0][0] >= 20) 
                    if(!((*data)->obs_statistics[13] = (double *)realloc((*data)->obs_statistics[13],(unsigned)(numb_c[0][0]+1)*sizeof(double))))
                        perror("realloc error.9");
                for(z=0;z<numb_c[0][0]+1;z++) (*data)->obs_statistics[13][z] = (double)numb_c[0][z];
                break;
            case 81:
                if(numb_c[0][0] >= 20) 
                    if(!((*data)->obs_statistics[14] = (double *)realloc((*data)->obs_statistics[14],(unsigned)(numb_c[0][0]+1)*sizeof(double))))
                        perror("realloc error.9");
                for(z=0;z<numb_c[0][0]+1;z++) (*data)->obs_statistics[14][z] = (double)numb_c[0][z];
                break;
            case 82:
                if(numb_c[0][0] >= 20) 
                    if(!((*data)->obs_statistics[15] = (double *)realloc((*data)->obs_statistics[15],(unsigned)(numb_c[0][0]+1)*sizeof(double))))
                        perror("realloc error.9");
                for(z=0;z<numb_c[0][0]+1;z++) (*data)->obs_statistics[15][z] = (double)numb_c[0][z];
                break;
            case 83:
                if(numb_c[0][0] >= 20) 
                    if(!((*data)->obs_statistics[16] = (double *)realloc((*data)->obs_statistics[16],(unsigned)(numb_c[0][0]+1)*sizeof(double))))
                        perror("realloc error.9");
                for(z=0;z<numb_c[0][0]+1;z++) (*data)->obs_statistics[16][z] = (double)numb_c[0][z];
                break;
            case 84:
                if(numb_c[0][0] >= 20) 
                    if(!((*data)->obs_statistics[17] = (double *)realloc((*data)->obs_statistics[17],(unsigned)(numb_c[0][0]+1)*sizeof(double))))
                        perror("realloc error.9");
                for(z=0;z<numb_c[0][0]+1;z++) (*data)->obs_statistics[17][z] = (double)numb_c[0][z];
                break;
            case 85:
                if(numb_c[0][0] >= 20) 
                    if(!((*data)->obs_statistics[18] = (double *)realloc((*data)->obs_statistics[18],(unsigned)(numb_c[0][0]+1)*sizeof(double))))
                        perror("realloc error.9");
                for(z=0;z<numb_c[0][0]+1;z++) (*data)->obs_statistics[18][z] = (double)numb_c[0][z];
                break;
            case 86:
                if(numb_c[0][0] >= 20) 
                    if(!((*data)->obs_statistics[19] = (double *)realloc((*data)->obs_statistics[19],(unsigned)(numb_c[0][0]+1)*sizeof(double))))
                        perror("realloc error.9");
                for(z=0;z<numb_c[0][0]+1;z++) (*data)->obs_statistics[19][z] = (double)numb_c[0][z];
                break;
            case 87:
                if(numb_c[0][0] >= 20) 
                    if(!((*data)->obs_statistics[20] = (double *)realloc((*data)->obs_statistics[20],(unsigned)(numb_c[0][0]+1)*sizeof(double))))
                        perror("realloc error.9");
                for(z=0;z<numb_c[0][0]+1;z++) (*data)->obs_statistics[20][z] = (double)numb_c[0][z];
                break;
            case 88:
                if(numb_c[0][0] >= 20) 
                    if(!((*data)->obs_statistics[21] = (double *)realloc((*data)->obs_statistics[21],(unsigned)(numb_c[0][0]+1)*sizeof(double))))
                        perror("realloc error.9");
                for(z=0;z<numb_c[0][0]+1;z++) (*data)->obs_statistics[21][z] = (double)numb_c[0][z];
                break;
            case 89:
                if(numb_c[0][0] >= 20) 
                    if(!((*data)->obs_statistics[22] = (double *)realloc((*data)->obs_statistics[22],(unsigned)(numb_c[0][0]+1)*sizeof(double))))
                        perror("realloc error.9");
                for(z=0;z<numb_c[0][0]+1;z++) (*data)->obs_statistics[22][z] = (double)numb_c[0][z];
                break;
            case 90:
                if(numb_c[0][0] >= 20) 
                    if(!((*data)->obs_statistics[23] = (double *)realloc((*data)->obs_statistics[23],(unsigned)(numb_c[0][0]+1)*sizeof(double))))
                        perror("realloc error.9");
                for(z=0;z<numb_c[0][0]+1;z++) (*data)->obs_statistics[23][z] = (double)numb_c[0][z];
                break;
            case 91:
                if(numb_c[0][0] >= 20) 
                    if(!((*data)->obs_statistics[24] = (double *)realloc((*data)->obs_statistics[24],(unsigned)(numb_c[0][0]+1)*sizeof(double))))
                        perror("realloc error.9");
                for(z=0;z<numb_c[0][0]+1;z++) (*data)->obs_statistics[24][z] = (double)numb_c[0][z];
                break;
            case 92:
                if(numb_c[0][0] >= 20) 
                    if(!((*data)->linked_segsites = (int *)realloc((*data)->linked_segsites,(unsigned)(numb_c[0][0]+1)*sizeof(int))))
                        perror("realloc error.93");
                for(z=0;z<numb_c[0][0]+1;z++) (*data)->linked_segsites[z] = (int)numb_c[0][z];
                break;
            case 93:
                if(numb_c[0][0] >= 20) 
                    if(!((*data)->linked_rm = (int *)realloc((*data)->linked_rm,(unsigned)(numb_c[0][0]+1)*sizeof(int))))
                        perror("realloc error.93");
                for(z=0;z<numb_c[0][0]+1;z++) (*data)->linked_rm[z] = (int)numb_c[0][z];
                break;
            case 94:
                if(numb_c[0][0] >= 20) 
                    if(!((*data)->heter_theta_alphag = (double *)realloc((*data)->heter_theta_alphag,(unsigned)(numb_c[0][0]+1)*sizeof(double))))
                        perror("realloc error.93");
                for(z=0;z<numb_c[0][0]+1;z++) (*data)->heter_theta_alphag[z] = (double)numb_c[0][z];
                break;
            case 95:
                if(numb_c[0][0] >= 20) 
                    if(!((*data)->heter_rm_alphag = (double *)realloc((*data)->heter_rm_alphag,(unsigned)(numb_c[0][0]+1)*sizeof(double))))
                        perror("realloc error.93");
                for(z=0;z<numb_c[0][0]+1;z++) (*data)->heter_rm_alphag[z] = (double)numb_c[0][z];
                break;
            case 96:
                if(numb_c[0][0] >= 20) 
                    if(!((*data)->linked_nhapl = (int *)realloc((*data)->linked_nhapl,(unsigned)(numb_c[0][0]+1)*sizeof(int))))
                        perror("realloc error.93");
                for(z=0;z<numb_c[0][0]+1;z++) (*data)->linked_nhapl[z] = (int)numb_c[0][z];
                break;
            case 97:
                if(numb_c[0][0] >= 20) 
                    if(!((*data)->correct_gamma = (double *)realloc((*data)->correct_gamma,(unsigned)(numb_c[0][0]+1)*sizeof(double))))
                        perror("realloc error.97");
                for(z=0;z<numb_c[0][0]+1;z++) (*data)->correct_gamma[z] = (double)numb_c[0][z];
                break;
            case 98:
                if(numb_c[0][0] >= 20) 
                    if(!((*data)->correct_gammar = (double *)realloc((*data)->correct_gammar,(unsigned)(numb_c[0][0]+1)*sizeof(double))))
                        perror("realloc error.98");
                for(z=0;z<numb_c[0][0]+1;z++) (*data)->correct_gammar[z] = (double)numb_c[0][z];
                break;
            case 99:
                if(numb_c[0][0] >= 20) 
                    if(!((*data)->invariable_mut_sites = (long int *)realloc((*data)->invariable_mut_sites,(unsigned)(numb_c[0][0]+1)*sizeof(long int))))
                        perror("realloc error.98");
                for(z=0;z<numb_c[0][0]+1;z++) (*data)->invariable_mut_sites[z] = (long int)numb_c[0][z];
                break;
            case 100:
                (*data)->likelihood_line = (int)numb_c[0][1];
                break;
            case 101:
                (*data)->likelihood_error[0] = (double)numb_c[0][1];
                break;
            case 102:
                (*data)->likelihood_error[1] = (double)numb_c[0][1];
                break;
            case 103:
                (*data)->likelihood_error[2] = (double)numb_c[0][1];
                break;
            case 104:
                (*data)->likelihood_error[3] = (double)numb_c[0][1];
                break;
            case 105:
                (*data)->likelihood_error[4] = (double)numb_c[0][1];
                break;
            case 106:
                (*data)->likelihood_error[5] = (double)numb_c[0][1];
                break;
            case 107:
                (*data)->likelihood_error[6] = (double)numb_c[0][1];
                break;
            case 108:
                (*data)->likelihood_error[7] = (double)numb_c[0][1];
                break;
            case 109:
                (*data)->likelihood_error[8] = (double)numb_c[0][1];
                break;
            case 110:
                (*data)->likelihood_error[9] = (double)numb_c[0][1];
                break;
            case 111:
                (*data)->likelihood_error[10] = (double)numb_c[0][1];
                break;
            case 112:
                (*data)->likelihood_error[11] = (double)numb_c[0][1];
                break;
            case 113:
                (*data)->likelihood_error[12] = (double)numb_c[0][1];
                break;
            case 114:
                (*data)->likelihood_error[13] = (double)numb_c[0][1];
                break;
            case 115:
                (*data)->likelihood_error[14] = (double)numb_c[0][1];
                break;
            case 116:
                (*data)->likelihood_error[15] = (double)numb_c[0][1];
                break;
            case 117:
                (*data)->likelihood_error[16] = (double)numb_c[0][1];
                break;
            case 118:
                (*data)->likelihood_error[17] = (double)numb_c[0][1];
                break;
            case 119:
                (*data)->likelihood_error[18] = (double)numb_c[0][1];
                break;
            case 120:
                (*data)->likelihood_error[19] = (double)numb_c[0][1];
                break;
            case 121:
                (*data)->likelihood_error[20] = (double)numb_c[0][1];
                break;
            case 122:
                (*data)->likelihood_error[21] = (double)numb_c[0][1];
                break;
            case 123:
                (*data)->likelihood_error[22] = (double)numb_c[0][1];
                break;
            case 124:
                (*data)->likelihood_error[23] = (double)numb_c[0][1];
                break;
            case 125:
                (*data)->likelihood_error[24] = (double)numb_c[0][1];
                break;
            case 126:
                if(numb_c[0][0] >= 20) 
                    if(!((*data)->obs_statistics[25] = (double *)realloc((*data)->obs_statistics[25],(unsigned)(numb_c[0][0]+1)*sizeof(double))))
                        perror("realloc error.9");
                for(z=0;z<numb_c[0][0]+1;z++) (*data)->obs_statistics[25][z] = (double)numb_c[0][z];
                break;
            case 127:
                (*data)->likelihood_error[25] = (double)numb_c[0][1];
                break;
            case 128:
                if(numb_c[0][0] >= 20) 
                    if(!((*data)->obs_statistics[26] = (double *)realloc((*data)->obs_statistics[26],(unsigned)(numb_c[0][0]+1)*sizeof(double))))
                        perror("realloc error.9");
                for(z=0;z<numb_c[0][0]+1;z++) (*data)->obs_statistics[26][z] = (double)numb_c[0][z];
                break;
            case 129:
                (*data)->likelihood_error[26] = (double)numb_c[0][1];
                break;
            case 130:
                if(numb_c[0][0] >= 20) 
                    if(!((*data)->obs_statistics[27] = (double *)realloc((*data)->obs_statistics[27],(unsigned)(numb_c[0][0]+1)*sizeof(double))))
                        perror("realloc error.9");
                for(z=0;z<numb_c[0][0]+1;z++) (*data)->obs_statistics[27][z] = (double)numb_c[0][z];
                break;
            case 131:
                (*data)->likelihood_error[27] = (double)numb_c[0][1];
                break;
            case 132:
                if(numb_c[0][0] >= 20) 
                    if(!((*data)->obs_statistics[28] = (double *)realloc((*data)->obs_statistics[28],(unsigned)(numb_c[0][0]+1)*sizeof(double))))
                        perror("realloc error.9");
                for(z=0;z<numb_c[0][0]+1;z++) (*data)->obs_statistics[28][z] = (double)numb_c[0][z];
                break;
            case 133:
                (*data)->likelihood_error[28] = (double)numb_c[0][1];
                break;
            case 134:
                if(numb_c[0][0] >= 20) 
                    if(!((*data)->obs_statistics[29] = (double *)realloc((*data)->obs_statistics[29],(unsigned)(numb_c[0][0]+1)*sizeof(double))))
                        perror("realloc error.9");
                for(z=0;z<numb_c[0][0]+1;z++) (*data)->obs_statistics[29][z] = (double)numb_c[0][z];
                break;
            case 135:
                (*data)->likelihood_error[29] = (double)numb_c[0][1];
                break;
            case 136:
                if(numb_c[0][0] >= 20) 
                    if(!((*data)->obs_statistics[30] = (double *)realloc((*data)->obs_statistics[30],(unsigned)(numb_c[0][0]+1)*sizeof(double))))
                        perror("realloc error.9");
                for(z=0;z<numb_c[0][0]+1;z++) (*data)->obs_statistics[30][z] = (double)numb_c[0][z];
                break;
            case 137:
                (*data)->likelihood_error[30] = (double)numb_c[0][1];
                break;
            default:
                break;

        }
    }
    free(f);
    for(x=0;x<=maxn;x++) free(numb_c[x]);
	free(numb_c);
    
    /* INIT SEED FOR RANDOM VALUES */
    if((*data)->seed1 < 0) {
        printf("Error: seed1 must be between 0 and 2147483562");
        exit(1);
    }
    init_seed1((*data)->seed1);
    srand((*data)->seed1);	/*emprem rands() per outgroup i mhits en especiacio*/
          
    /* FILTERS, AJ! VERY BAD !*/
    if((*data)->print_neuttest > 0) (*data)->neutral_tests = 1;

    if((*data)->mhits > 1) {
        printf("Error: mhits must be 0 or 1\n");
        exit(1);
    }
    if((*data)->T_out < 0.0) {
        printf("Error: dist_out must be >= 0\n");
        exit(1);
    }
    if((*data)->seed2 < 0) {
        printf("Error: seed2 must be between 0 and 2147483398\n");
        exit(1);
    }
    if((*data)->n_loci <  1) {
        printf("Error: n_loci must be at least 1\n");
        exit(1);
    }
    if((*data)->n_loci >  1 && (*data)->linked > 0) {
        printf("Error: When loci are linked, is not allowed more than 1 independent loci. \n");
        exit(1);
    }
    if((*data)->n_loci >  1 && (*data)->despl > 0. && (*data)->window > 0.) {
        printf("Error: When more than one loci are linked is not allowed sliding windows.\n ");
        exit(1);
    }
    if((*data)->linked !=  1 && (*data)->despl > 0. && (*data)->window > 0.) {
        printf("Error: Sliding windows is only allowed with linked value equal to 1. \n");
        exit(1);
    }
    if((*data)->linked ==  1 && ((*data)->despl == 0. || (*data)->window == 0.)) {
        printf("Error: When pos_linked value is equal to 1, sliding windows have to be defined. \n");
        exit(1);
    }
	if((*data)->linked > 0) {
		if((*data)->loci_linked[(*data)->linked-1][2] >= (*data)->nsites[1]) {
			printf("Error: size of the linked loci exceeds the total size of the region (0-%ld)\n",(*data)->nsites[1]-1);
			exit(1);
		}
		for(x=0;x<(*data)->linked;x++) {
			if((*data)->loci_linked[x][2] < (*data)->loci_linked[x][1]) {
				printf("Error: size of the linked loci %d is less than zero\n",x);
				exit(1);
			}
			if(x) {
				if((*data)->loci_linked[x][1] < (*data)->loci_linked[x-1][2]) {
					printf("Error: linked locus not sorted\n");
					exit(1);
				}
			}
		}
		if((*data)->linked > 1) {
			if((*data)->linked_segsites[0] > 0) {
				if((*data)->linked_segsites[0] != (*data)->linked) {
						printf("Error: not the same linked loci than defined in linked_segsites\n");
						exit(1);
				}
				if((*data)->mutations[0] == (*data)->n_loci && (*data)->mutations[1] != -1) {
						printf("Error: Only allowed linked_segsites or mutations (for the entire region). Not together\n");
						exit(1);
				}
			}
			else {
				if((*data)->linked_segsites[0] < (*data)->linked) {
					if(!((*data)->linked_segsites = (int *)realloc((*data)->linked_segsites,(unsigned)((*data)->linked+1)*sizeof(int)))) 
						perror("realloc error.93");
					(*data)->linked_segsites[0] = (*data)->linked;
					for(x=1;x<(*data)->linked+1;x++) (*data)->linked_segsites[x] = -1;
				}
			}
			if((/*(*data)->linked_segsites[0] > 0 || */(*data)->linked_segsites[1] != (double)-1) && 
			   ((*data)->sfix_allthetas == 0 /*&& ((*data)->theta_1[0] != (*data)->n_loci) || (*data)->theta_1[1] <= (double)0*/)) {
				printf("Error: Not allowed to fix segsites in linked fragments wihthout the option 'Sfix_allthetas'.\n");
				exit(1);
			}
			
			if((*data)->linked_rm[0] > 0) {
				if((*data)->linked_rm[0] != (*data)->linked) {
						printf("Error: not the same linked loci than defined in linked_rm\n");
						exit(1);
				}
				if((*data)->Rm[0] == (*data)->n_loci && (*data)->Rm[0] != -1) {
						printf("Error: Only allowed linked_rm or rm (for the entire region). Not together\n");
						exit(1);
				}
			}
			else {
				if((*data)->linked_rm[0] < (*data)->linked) {
					if(!((*data)->linked_rm = (int *)realloc((*data)->linked_rm,(unsigned)((*data)->linked+1)*sizeof(int)))) 
						perror("realloc error.93");
					(*data)->linked_rm[0] = (*data)->linked;
					for(x=1;x<(*data)->linked+1;x++) (*data)->linked_rm[x] = -1;
				}
			}
			if((*data)->linked_nhapl[0] > 0) {
				if((*data)->linked_nhapl[0] != (*data)->linked) {
						printf("Error: not the same linked loci than defined in linked_rm\n");
						exit(1);
				}
			}
			else {
				if((*data)->linked_nhapl[1] < (*data)->linked) {
					if(!((*data)->linked_nhapl = (int *)realloc((*data)->linked_nhapl,(unsigned)((*data)->linked+1)*sizeof(int)))) 
						perror("realloc error.93");
					(*data)->linked_nhapl[0] = (*data)->linked;
					for(x=1;x<(*data)->linked+1;x++) (*data)->linked_nhapl[x] = 0;
				}
			}
			if((*data)->linked_rm[1] != (double)-1 && (*data)->rmfix == 0) {
				printf("Error: Not allowed to fix Rm in linked fragments wihthout the option 'rmfix'.\n");
				exit(1);
			}
			if((*data)->linked_nhapl[1] != (double)0 && (*data)->rmfix == 0) {
				printf("Error: Not allowed to fix the number of haplotypes in linked fragments wihthout the option 'rmfix'.\n");
				exit(1);
			}
		}
    }		
    if((*data)->nsam[0] != (*data)->n_loci) {
        printf("Error: n_samples needs n_loci inputs\n");
        exit(1);
    }
    /*
	if((*data)->ratio_sv[0] != 0) {
        if(!((*data)->mhits && (*data)->theta_1[1] > 0.)) {
            printf("Warning: mhits transitions/transversions is only used with mhits and theta options. Not used.\n\n");
        }
    }
	*/
    if((*data)->mhits) {
		/*
        for(x=0;x<(*data)->n_loci;x++) {
            if((*data)->nsam[x+1] == 2) {
                printf("Error: sample=2, mhits is not enabled");
                exit(1);
            }
        }
		*/
        if(((*data)->theta_1[0] > 0 && (*data)->theta_1[1] > 0.) || (*data)->sfix_allthetas || (*data)->range_thetant || (*data)->ifgamma) {
            if((*data)->ratio_sv[0] < (*data)->n_loci) {
                if(!((*data)->ratio_sv = (double *)realloc((*data)->ratio_sv,(unsigned)((*data)->n_loci+1)*sizeof(double)))) 
                    perror("realloc error.33");
                (*data)->ratio_sv[0] = (*data)->n_loci;
                if((*data)->ratio_sv[1] == 0. || (*data)->ratio_sv[1] == 0.5) (*data)->ratio_sv[1] = 0.5;
                for(x=2;x<(*data)->n_loci+1;x++) (*data)->ratio_sv[x] = (*data)->ratio_sv[1];
            }
            for(x=1;x<(*data)->n_loci+1;x++) {
                if((*data)->ratio_sv[x] < 0.) {
                    printf("Error: mhratio_sv must be higher than 0.\n");
                    exit(1);
                }
            }
        }
    }
    if((*data)->r[0]) {
		if((*data)->r[0] != (*data)->n_loci) {
			printf("Error: recombination input error\n");
			exit(1);
		}
    }
	
    if((*data)->heter_theta_alphag[0] > 0) {
		if((*data)->heter_theta_alphag[0] != (*data)->n_loci) {
				printf("Error: not the same loci than defined in heter_theta_alphag\n");
				exit(1);
		}
	}
	else {
		if(!((*data)->heter_theta_alphag = (double *)realloc((*data)->heter_theta_alphag,(unsigned)((*data)->n_loci+1)*sizeof(double)))) 
			perror("realloc error.93");
		(*data)->heter_theta_alphag[0] = (double)(*data)->n_loci;
		for(x=1;x<(*data)->n_loci+1;x++) (*data)->heter_theta_alphag[x] = (double)-1;
	}
		
    if((*data)->invariable_mut_sites[0] > 0) {
		if((*data)->invariable_mut_sites[0] != (*data)->n_loci) {
				printf("Error: not the same loci than defined in invariable_mut_sites\n");
				exit(1);
		}
		else {
			for(x=1;x<(*data)->n_loci+1;x++) {
				if((*data)->invariable_mut_sites[x] >= (*data)->nsites[x]) {
					printf("Error: invariable_mut_sites must be smaller than n_sites.\n");
					exit(1);
				}
			}
		}
	}
	else {
		if(!((*data)->invariable_mut_sites = (long int *)realloc((*data)->invariable_mut_sites,(long int)((*data)->n_loci+1)*sizeof(long int)))) 
			perror("realloc error.99");
		(*data)->invariable_mut_sites[0] = (long int)(*data)->n_loci;
		for(x=1;x<(*data)->n_loci+1;x++) (*data)->invariable_mut_sites[x] = (long int)-1;
	}

    if((*data)->heter_rm_alphag[0] > 0) {
		if((*data)->heter_rm_alphag[0] != (*data)->n_loci) {
				printf("Error: not the same loci than defined in heter_rm_alphag\n");
				exit(1);
		}
	}
	else {
		if(!((*data)->heter_rm_alphag = (double *)realloc((*data)->heter_rm_alphag,(unsigned)((*data)->n_loci+1)*sizeof(double)))) 
			perror("realloc error.93");
		(*data)->heter_rm_alphag[0] = (double)(*data)->n_loci;
		for(x=1;x<(*data)->n_loci+1;x++) (*data)->heter_rm_alphag[x] = (double)-1;
	}
	/*
    if((*data)->r[0]==0) {
		if(!((*data)->r = (double *)realloc((*data)->r,(unsigned)(1+1)*sizeof(double)))) {
			perror("realloc error.33");
			exit(1);
		}
		(*data)->r[0] = (double)1;
		(*data)->r[1] = (double)0;
    }
	*/
    if((int)(*data)->nsites[0] != (*data)->n_loci) {
        printf("Error: nsites needs n_loci inputs\n");
        exit(1);
    }
    if((*data)->f[0] >0 && (*data)->f[1] >0) {
            if((*data)->f[0] == 0) {
                printf("Error: recombination input error\n");
                exit(1);
            }
    }
    if((*data)->track_len[0] && (*data)->track_len[1] > 0) {
            if((*data)->track_len[0] == 0) {
                printf("Error: recombination input error\n");
                exit(1);
            }
    }
    if((*data)->theta_1[0] > 0 && !((*data)->theta_1[0] == 1 && (*data)->theta_1[1] <= 0.)) {
            if((*data)->theta_1[0] == 0) {
                printf("Error: theta_1 needs n_loci inputs\n");
                exit(1);
            }
    }
    if((*data)->mutations[0] > 0 && !((*data)->mutations[0] == 1 && (*data)->mutations[1] == -1)) {
        if((int)(*data)->mutations[0] != (*data)->n_loci) {
            printf("Error: mutations needs n_loci inputs\n");
            exit(1);
        }
    }
    if((*data)->mutations[0] == 0 || ((*data)->mutations[0] == 1 && (*data)->mutations[1] == -1)) {
       if(!((*data)->mutations = (long int *)realloc((*data)->mutations,((*data)->n_loci+1)*sizeof(long int)))) 
            perror("realloc error.38bdb");
		(*data)->mutations[0] = (*data)->n_loci;
		for(z=1;z<(*data)->n_loci+1;z++) (*data)->mutations[z] = -1;
    }
    if((*data)->split_pop == 1) {
        if((*data)->n_loci > 20)
            if(!((*data)->ssize_pop = (int **)realloc((*data)->ssize_pop,(*data)->n_loci*sizeof(int *)))) 
                perror("realloc error.38bb");
        if(!((*data)->npop_sampled = (int *)realloc((*data)->npop_sampled,((*data)->n_loci+1)*sizeof(int)))) 
            perror("realloc error.38bb");/*afegit21.4.2003*/
        (*data)->npop_sampled[0] = (*data)->n_loci;/*afegit21.4.2003*/
        for(z=0;z<(*data)->n_loci;z++) {
            if(z<20) {
                if(!((*data)->ssize_pop[z]=(int *)realloc((*data)->ssize_pop[z],2*sizeof(int)))) 
                    perror("realloc error.22bb");
            }
            else {
                if(!((*data)->ssize_pop[z]=(int *)malloc(2*sizeof(int)))) 
                    perror("realloc error.22bb");
            }
        }
        for(x=0;x<(*data)->n_loci;x++) {
            (*data)->ssize_pop[x][0] = 1;
            (*data)->ssize_pop[x][1] = (*data)->nsam[x+1];
            (*data)->npop_sampled[x+1] = 1;/*afegit 21.4.2003*/
        }
    }
    
    if((*data)->npop == 1 && (*data)->split_pop == 0) {
        if((*data)->n_loci > 20)
            if(!((*data)->ssize_pop = (int **)realloc((*data)->ssize_pop,((*data)->n_loci)*sizeof(int *)))) 
                perror("realloc error.38bb");
        if(!((*data)->npop_sampled = (int *)realloc((*data)->npop_sampled,((*data)->n_loci+1)*sizeof(int)))) 
            perror("realloc error.38bb");/*afegit21.4.2003*/
        (*data)->npop_sampled[0] = (*data)->n_loci;/*afegit21.4.2003*/
        for(z=0;z<(*data)->n_loci;z++) {
            if(z<20) {
                if(!((*data)->ssize_pop[z]=(int *)realloc((*data)->ssize_pop[z],2*sizeof(int)))) 
                    perror("realloc error.22bb");
            }
            else {
                if(!((*data)->ssize_pop[z]=(int *)malloc(2*sizeof(int)))) 
                    perror("realloc error.22bb");
            }
        }
        for(x=0;x<(*data)->n_loci;x++) {
            (*data)->ssize_pop[x][0] = 1;
            (*data)->ssize_pop[x][1] = (*data)->nsam[x+1];
            (*data)->npop_sampled[x+1] = 1;/*afegit 21.4.2003*/
        }
        (*data)->ran_factorpop = 0;
        (*data)->same_factorpop = 1;
        (*data)->factor_pop[0] = 1;
        (*data)->factor_pop[1] = 1;
    }
    if((*data)->npop > 0  && (*data)->split_pop == 0) {
        if((*data)->npop > 1) {
            if((*data)->npop_sampled[0] == 1 && (*data)->n_loci > 1) { 
                /*Si nomŽs hi ha 1 valor Žs que no s«han definit els altres loci*/
                if(!((*data)->npop_sampled = (int *)realloc((*data)->npop_sampled,((*data)->n_loci+1)*sizeof(int))))
                    perror("realloc error.38bb");
                (*data)->npop_sampled[0] = (*data)->n_loci; /*definim tots els loci*/
                for(z=1;z<(*data)->n_loci;z++) (*data)->npop_sampled[z+1] = (*data)->npop_sampled[1];
                /*tots els loci tenen la mateixa mostra que el primer loci*/
            }
			if((*data)->ssize_pop[0][0] == 0) {
				/*if ssize is not included in input file*/
				(*data)->ssize_pop[0][0] = (*data)->npop_sampled[0];
				if(!((*data)->ssize_pop[0] = (int *)realloc((*data)->ssize_pop[0],(unsigned)((*data)->npop_sampled[0]+1)*sizeof(int)))) perror("calloc error.11b");
				(*data)->ssize_pop[0][1] = (*data)->nsam[0+1];
				for(x=2;x<=(*data)->ssize_pop[0][0];x++) {
					(*data)->ssize_pop[0][x] = 0;
				}
			}
			if(((*data)->ssize_pop[1][0] == 0 && (*data)->n_loci > 1) || ((*data)->ssize_pop[0][0] < (*data)->npop_sampled[0])) {
				/*if only the first loci is defined or if ssize_pop is not entirey defined*/
				if(!((*data)->ssize_pop = (int **)realloc((*data)->ssize_pop,(unsigned)(*data)->n_loci*sizeof(int *)))) perror("calloc error.11");
				/*define size for ssize_pop per loci*/
				for(x=0;x<(*data)->n_loci;x++) {
					if((*data)->n_loci <= 20) {
						if(!((*data)->ssize_pop[x] = (int *)realloc((*data)->ssize_pop[x],(unsigned)((*data)->npop_sampled[0]+1)*sizeof(int)))) perror("calloc error.11b");
					}
					else {
						if(!((*data)->ssize_pop[x] = (int *)calloc((unsigned)(*data)->npop_sampled[0]+1,sizeof(int)))) perror("calloc error.11b");
					}	
					(*data)->ssize_pop[x][0] = (*data)->npop_sampled[0];
					(*data)->ssize_pop[x][1] = (*data)->nsam[x+1];
					for(z=2;z<=(*data)->ssize_pop[x][0];z++) {
						(*data)->ssize_pop[x][z] = 0;
					}
				}
			}			
        }
        else {
            if((*data)->npop == 1) {
                if(!((*data)->npop_sampled = (int *)realloc((*data)->npop_sampled,((*data)->n_loci+1)*sizeof(int))))
                    perror("realloc error.38bb");
                for(x=0;x<(*data)->n_loci;x++) { /*definir valors*/
                    if(!((*data)->ssize_pop[x]=(int *)realloc((*data)->ssize_pop[x],2*sizeof(int)))) 
                        perror("realloc error.22bb");
                    (*data)->ssize_pop[x][0] = 1;
                    (*data)->ssize_pop[x][1] = (*data)->nsam[x+1];
                    (*data)->npop_sampled[x+1] = 1;
                }
            }
        }
    }

    if((*data)->nintn > 0) {
        if((*data)->nrec[0] != (*data)->nintn) {
            printf("Error: nrec has not the same intervals indicated in nintn\n");
            exit(1);
        }
        if((*data)->nrec[1] != 1.0) {
            printf("Error: nrec must be 1.0 in the first value\n");
            exit(1);
        }
        if((*data)->npast[0] != (*data)->nintn) {
			if((*data)->npast[0] == 0) {
				if(!((*data)->npast=(double *)realloc((*data)->npast,((*data)->nintn+1)*sizeof(double)))) {
					printf("Error: malloc error in npast\n");
					exit(1);
				}
				for(x=0;x<(*data)->nintn;x++) {
					(*data)->npast[x+1] = (*data)->nrec[x+1];
				}
				(*data)->npast[0] = (*data)->nintn;
			}
			else {
				printf("Error: npast has not the same intervals indicated in nintn\n");
				exit(1);
			}
        }
        if((*data)->tpast[0] != (*data)->nintn) {
			printf("Error: tpast has not the same intervals indicated in nintn\n");
			exit(1);
        }
    }

    if((*data)->ran_factorpop < 0 || (*data)->ran_factorpop > 1) {
        printf("Error ran_factorpop: it should be  0 or 1\n");
        exit(1);
    }
    if((*data)->same_factorpop < 0 || (*data)->same_factorpop > 1) {
        printf("Error ran_factorpop: it should be 0 or 1\n");
        exit(1);
    }
    if((*data)->ran_factorpop != 0 && (*data)->same_factorpop != 0) {
        printf("Error ran_factorpop and same_factorpop can not be activated at the same time\n");
        exit(1);
    }
    if((*data)->neutral_tests < 0 || (*data)->neutral_tests > 1) {
        printf("Error neutral_tests: it should be 0 or 1\n");
        exit(1);
    }
    for(x=0;x<(*data)->n_loci;x++) {
        y = 0;
        for(z=1;z<(*data)->ssize_pop[x][0]+1;z++) y += (*data)->ssize_pop[x][z];
        if((*data)->nsam[x+1] != y) {
            printf("Error: npop_sampled is not coincident with the total number of samples\n");
            exit(1);
        }
    }
    for(x=0;x<(*data)->n_loci;x++) {
        if((*data)->npop_sampled[x+1] > 0 && (*data)->npop_sampled[x+1] != (*data)->ssize_pop[x][0]) {
            printf("Error: npop_sampled is not coincident with the #pops in ssize_pop\n");
            exit(1);
        }
    }
    /*vull tenir un nombre de poblacions mostrejades per cada locus i que tingui un nombre de mostres per poblaci—*/
    /*tots els loci han de tenir el mateix nombre de poblacions, pero el nombre de mostres poden ser diferents*/
    for(x=0;x<(*data)->n_loci;x++)/*aqu’ es posen els valors de npop_sampled en cas no estigui definit...*/
        if((*data)->npop_sampled[x+1] == 0)
            (*data)->npop_sampled[x+1] = (*data)->ssize_pop[x][0];
    (*data)->npop_sampled[0] = (*data)->n_loci;
    /*De fet, el que vull es que en les mostres, en alguns loci es fan servir poblacion diferent i hi han zeros pel mig */
    /*Aix˜ no ho se per qu ho he definit, miro mŽs endavant... ƒS IMPORTANT!! */
    if((*data)->ssize_pop[0][0] < (*data)->npop) {   
        for(x=0;x<(*data)->n_loci;x++) {
            if(!((*data)->ssize_pop[x] = (int *)realloc((*data)->ssize_pop[x],((*data)->npop+1)*sizeof(int))))
                perror("realloc error.78");
            for(z=(*data)->ssize_pop[x][0]+1;z<(*data)->npop+1;z++) (*data)->ssize_pop[x][z] = 0;
            (*data)->ssize_pop[x][0] = (*data)->npop;
        }
    }
    /**/
    if((*data)->factor_pop[0] > 0 && (*data)->factor_pop[1] > 0.0) {
        for(x=1;x<(*data)->factor_pop[0]+1;x++) {
            if((*data)->factor_pop[x] <= 0.0 ) {
                printf("Error factor_pop values: they should be higher than 0\n");
                exit(1);
            }
        }
        if((*data)->factor_pop[1] != 1.0 && (*data)->split_pop == 0) {
            printf("Error factor_pop vector: First value must be 1.0\n");
            exit(1);
        }
        /* popsizes: 1 or between 1 and 10 relative to the first pop if all are not indicated*/
        /**/
        if((*data)->factor_pop[0] < (*data)->npop) {
            if(!((*data)->factor_pop = (double *)realloc((*data)->factor_pop,((*data)->npop+1)*sizeof(double))))
                perror("realloc error.1");
            if((*data)->ran_factorpop == 0 && (*data)->same_factorpop == 0) z = (int)(*data)->factor_pop[0]+1;
            else z = 1;
            (*data)->factor_pop[0] = (*data)->npop;
            for(x=z;x<(*data)->npop+1;x++) {
                if((*data)->same_factorpop == 1 || ((*data)->ran_factorpop == 0 && (*data)->same_factorpop == 0)) (*data)->factor_pop[x] = 1.0;
                else {
                    if(x==1) (*data)->factor_pop[x] = 1.0;
                    else {
                        (*data)->factor_pop[x] = (double)ran1()*9.0 + 1.0; 
                        if((double)ran1() < 0.5)  (*data)->factor_pop[x] = 1./(*data)->factor_pop[x];
                    }
                }
            }
        }
        /**/
    }
	
    if((*data)->theta_1[0] == 0 || ((*data)->theta_1[0] == 1 && (*data)->theta_1[1] <= (double)0)) {
		if((*data)->mutations[0] == 0 || ((*data)->mutations[0] == 1 && (*data)->mutations[1] == -1 && (int)(*data)->mutations[0] != (*data)->n_loci)) {
			if((*data)->range_thetant == 0 && (*data)->ifgamma == 0) {
                printf("Error: thetaw, mutations, range_thetant and/or ifgamma must be defined\n");
                exit(1);
			}
		}
	}
	
    if((*data)->r[0] == 0 || ((*data)->r[0] == 1 && (*data)->r[1] <= 0.  && (*data)->r[0] != (*data)->n_loci)) {
		if((*data)->range_rnt == 0 && (*data)->ifgammar == 0) {
			printf("Error: recombination, range_rnt or ifgammar must be defined\n");
			exit(1);
		}
	}
	
    if((*data)->neutral_tests == 1) {
        if((*data)->pr_matrix != 0) {
            printf("Error: netral_tests and pr_matrix are incompatibles\n");
            exit(1);
        }
        if((*data)->print_neuttest < 0 || (*data)->print_neuttest > 2) {
            printf("Error print_neuttest: it should be between 0 and 2\n");
            exit(1);
        }
    }
    if((*data)->sfix_allthetas > 1) {
		printf("Error: sfix_allthetas can only be 0 or 1.\n");
		exit(1);
    }
    if((*data)->method_samp > 2) {
		printf("Error: method_samp can only be 0, 1 or 2.\n");
		exit(1);
    }

    if((*data)->ifgamma == 1 && ((*data)->theta_1[0] > 0 && !((*data)->theta_1[0] > 0 && (*data)->theta_1[1] <= 0))) {
        printf("Error: 'theta' and 'ifgamma' can not be defined at the same time.\n");
        exit(1);
    }
    if((*data)->range_thetant && ((*data)->theta_1[0] > 0 && !(*data)->theta_1[0] == 1 && (*data)->theta_1[1] <= 0)) {
        printf("Error: 'theta' and 'range_thetant' can not be defined at the same time.\n");
        exit(1);
    }
    if((*data)->ifgammar == 1 && ((*data)->r[0] > 0 && !((*data)->r[0] > 0 && (*data)->r[1] == 0))) {
        printf("Error: 'recombination' and 'ifgammar' can not be defined at the same time.\n");
        exit(1);
    }
    if((*data)->range_rnt && ((*data)->r[0] > 0 && !(*data)->r[0] == 1 && (*data)->r[1] == 0)) {
        printf("Error: 'recombination' and 'range_rnt' can not be defined at the same time.\n");
        exit(1);
    }

    if((*data)->sfix_allthetas && ((*data)->theta_1[0] > 0 && !((*data)->theta_1[0] == 1 && (*data)->theta_1[1] <= 0.))) {
        printf("Error: 'Thetaw' and 'sfix_allthetas' can not be defined at the same time.\n");
        exit(1);
    }
    if(((*data)->sfix_allthetas == 1 || (*data)->rmfix == 1) && (*data)->method_samp == 2) {
        if((*data)->mc_jump < 1){
            printf("Error: mc_jump must be a positive integer\n");
            exit(1);
        }
    }
    if((*data)->sfix_allthetas) {
        if((*data)->range_thetant != 0 && (*data)->range_thetant != 1 && (*data)->range_thetant != 2) {
            printf("Error: range_recnt must be 0, 1 or 2\n");
            exit(1);
        }
        else {
            if((*data)->range_thetant == 1 || (*data)->range_thetant == 2) {
                if((*data)->thetant_min[0] == (*data)->thetant_max[0] && (*data)->thetant_min[0] > 0) {
					for(x=1;x<=(*data)->thetant_min[0];x++) {
						if((*data)->thetant_min[x] < 0) {
							printf("Error: thetant_min must be positive or a zero value\n");
							exit(1);
						}
						if((*data)->thetant_max[x] < (*data)->thetant_min[x]) {
							printf("Error: thetant_max must be positive or a zero value, and lower or equal than thetant_min\n\n");
							exit(1);
						}
					}
				}
				else {
					printf("Error: thetant_min and thetant_max must be defined\n");
					exit(1);
				}
            }
			if((*data)->range_thetant == 0) {
				if((*data)->ifgamma == 0) {
                    printf("Error in input: options sfix_allthetas, range_thetant and ifgamma. \n");
                    exit(1);
				}
			}
        }
    }
    if((*data)->rmfix && ((*data)->r[0] > 0 && !((*data)->r[0] == 1 && (*data)->r[1] == (double)0))) {
        printf("Error: 'recombination' and 'rmfix' can not be defined at the same time.\n");
        exit(1);
    }
    
    if((*data)->ifselection[0] < (*data)->n_loci) { 
        /*Si nomŽs hi ha 1 valor Žs que no s«han definit els altres loci*/
        if(!((*data)->ifselection = (int *)realloc((*data)->ifselection,((*data)->n_loci+1)*sizeof(int))))
            perror("realloc error.38bb");
        for(z=(*data)->ifselection[0];z<(*data)->n_loci;z++) (*data)->ifselection[z+1] = 0;/*no selection*/
        (*data)->ifselection[0] = (*data)->n_loci; /*definim tots els loci*/
    }
    if((*data)->r[0] < (*data)->n_loci) { 
        /*Si nomŽs hi ha 1 valor Žs que no s«han definit els altres loci*/
        if(!((*data)->r = (double *)realloc((*data)->r,((*data)->n_loci+1)*sizeof(double))))
            perror("realloc error.38bb");
        for(z=0;z<(*data)->n_loci;z++)
			(*data)->r[z+1] = (double)0;/*no recombination*/
        (*data)->r[0] = (*data)->n_loci; /*definim tots els loci*/
    }

    if((*data)->split_pop ==1 && (*data)->time_scoal <= (*data)->time_split) {
        printf("Error: time_scoal must be larger than time_split\n");
        exit(1);
    }
	
    if((*data)->split_pop == 1 && (*data)->freq[0] != (*data)->npoprefugia) {
        printf("Error: freq_refugia must be defined for all refugia\n");
        exit(1);
    }

    if((*data)->split_pop == 1 && (*data)->factor_pop[0] != (*data)->npoprefugia) {
        printf("Error: factor_pop must be defined for all refugia\n");
        exit(1);
    }

    if((*data)->split_pop == 0 && (*data)->npop > 1 && (*data)->factor_pop[0] != (*data)->npop) {
        printf("Error: factor_pop must be defined for all populations\n");
        exit(1);
    }
	
    if((*data)->split_pop == 0 && (*data)->npop > 1 && (*data)->mig_rate == 0.) {
        printf("Error: migration must be defined.\n");
        exit(1);
    }
	if((*data)->ifgamma == 1) {
		if((*data)->p_gamma[0] != (double)(*data)->n_loci) {
			printf("Error: p gamma parameter must be defined for ALL loci and be higher than 0.\n");
			exit(1);
		}
		if((*data)->alpha_gamma[0] != (double)(*data)->n_loci) {
			printf("Error: alpha gamma parameter must be defined for ALL loci and be higher than 0.\n");
			exit(1);
		}
		if((*data)->correct_gamma[0] == (double)0) {
			if(!((*data)->correct_gamma = (double *)realloc((*data)->correct_gamma,(unsigned)((*data)->n_loci+1)*sizeof(double)))) 
				perror("realloc error.155");
			(*data)->correct_gamma[0] = (double)(*data)->n_loci;
			for(x=1;x<=(*data)->n_loci;x++) {
				(*data)->correct_gamma[x] = (double)1.0;
			}
		}
		if((*data)->correct_gamma[0] != (double)(*data)->n_loci) {
			printf("Error: correction gamma parameter must be defined for ALL loci and be higher than 0.\n");
			exit(1);
		}
	}
	for(x=1;x<=(*data)->n_loci;x++) {
		if(((*data)->ifgamma == 1 && (*data)->sfix_allthetas > 0) && ((*data)->p_gamma[x] <= (double)0.0 || (*data)->alpha_gamma[x] <= (double)0. || (*data)->correct_gamma[x] <= (double)0.)) {
			printf("Error: gamma parameters must be defined for ALL loci and be higher than 0.\n");
			exit(1);
		}
	}
	
	if((*data)->factor_chrn[0] == (double)0) {
		if(!((*data)->factor_chrn = (double *)realloc((*data)->factor_chrn,(unsigned)((*data)->n_loci+1)*sizeof(double)))) 
			perror("realloc error.155");
		(*data)->factor_chrn[0] = (double)(*data)->n_loci;
		for(x=1;x<=(*data)->n_loci;x++) {
			(*data)->factor_chrn[x] = (double)1.0;
		}
	}
	if((*data)->factor_chrn[0] < (*data)->n_loci) {
		if(!((*data)->factor_chrn = (double *)realloc((*data)->factor_chrn,(unsigned)((*data)->n_loci+1)*sizeof(double)))) 
			perror("realloc error.155");
		for(x=(int)(*data)->factor_chrn[0];x<=(int)(*data)->n_loci;x++) {
			(*data)->factor_chrn[x] = (double)1.0;
		}
		(*data)->factor_chrn[0] = (double)(*data)->n_loci;
	}

	if((*data)->ratio_sv[0] == (double)0) {
		if(!((*data)->ratio_sv = (double *)realloc((*data)->ratio_sv,(unsigned)((*data)->n_loci+1)*sizeof(double)))) 
			perror("realloc error.155");
		(*data)->ratio_sv[0] = (double)(*data)->n_loci;
		for(x=1;x<=(*data)->n_loci;x++) {
			(*data)->ratio_sv[x] = (double)0.50;
		}
	}
	if((*data)->ratio_sv[0] < (*data)->n_loci) {
		if(!((*data)->ratio_sv = (double *)realloc((*data)->ratio_sv,(unsigned)((*data)->n_loci+1)*sizeof(double)))) 
			perror("realloc error.155");
		for(x=(int)(*data)->ratio_sv[0];x<=(int)(*data)->n_loci;x++) {
			(*data)->ratio_sv[x] = (double)0.50;
		}
		(*data)->ratio_sv[0] = (double)(*data)->n_loci;
	}

	if((*data)->p_gamma[0] == (double)0) {
		if(!((*data)->p_gamma = (double *)realloc((*data)->p_gamma,(unsigned)((*data)->n_loci+1)*sizeof(double)))) 
			perror("realloc error.155");
		(*data)->p_gamma[0] = (double)(*data)->n_loci;
		for(x=1;x<=(*data)->n_loci;x++) {
			(*data)->p_gamma[x] = (double)1.0;
		}
	}
	if((*data)->p_gamma[0] < (*data)->n_loci) {
		if(!((*data)->p_gamma = (double *)realloc((*data)->p_gamma,(unsigned)((*data)->n_loci+1)*sizeof(double)))) 
			perror("realloc error.155");
		for(x=(int)(*data)->p_gamma[0];x<=(int)(*data)->n_loci;x++) {
			(*data)->p_gamma[x] = (double)1.0;
		}
		(*data)->p_gamma[0] = (double)(*data)->n_loci;
	}

	if((*data)->alpha_gamma[0] == (double)0) {
		if(!((*data)->alpha_gamma = (double *)realloc((*data)->alpha_gamma,(unsigned)((*data)->n_loci+1)*sizeof(double)))) 
			perror("realloc error.155");
		(*data)->alpha_gamma[0] = (double)(*data)->n_loci;
		for(x=1;x<=(*data)->n_loci;x++) {
			(*data)->alpha_gamma[x] = (double)1.0;
		}
	}
	if((*data)->alpha_gamma[0] < (*data)->n_loci) {
		if(!((*data)->alpha_gamma = (double *)realloc((*data)->alpha_gamma,(unsigned)((*data)->n_loci+1)*sizeof(double)))) 
			perror("realloc error.155");
		for(x=(int)(*data)->alpha_gamma[0];x<=(int)(*data)->n_loci;x++) {
			(*data)->alpha_gamma[x] = (double)1.0;
		}
		(*data)->alpha_gamma[0] = (double)(*data)->n_loci;
	}

    if((*data)->mutations[0] > 0 && !((*data)->mutations[1] == -1)) {
		for(x=1;x<(int)(*data)->n_loci+1;x++) if((*data)->mutations[x] != 0) break;
		if(x < (*data)->n_loci) {
			if((*data)->sfix_allthetas == 0 && (*data)->mhits != 0) {
				printf("Error. mhits option is not allowed when theta is not defined (define thetaw or a distribution of theta)).\n");
				exit(1);
			}
		}
	}

    if((*data)->rmfix > 1) {
		printf("Error: rmfix can only be 0 or 1.\n");
		exit(1);
    }
    if((*data)->rmfix && (*data)->linked < 2) {
		if((*data)->Rm[0] != (*data)->n_loci) { 
			printf("Error: Once 'rmfix' is defined, 'Rm' must be defined for each locus.\n");
			exit(1);
		}
	}
	if((*data)->nhapl[0] == (int)0 && (*data)->rmfix) {
		if(!((*data)->nhapl = (int *)realloc((*data)->nhapl,(unsigned)((*data)->n_loci+1)*sizeof(int)))) 
			perror("realloc error.1558");
		(*data)->nhapl[0] = (int)(*data)->n_loci;
		for(x=1;x<=(*data)->n_loci;x++) {
			(*data)->nhapl[x] = (int)0;
		}
	}
    if((*data)->rmfix && (*data)->linked < 2) {
		if((int)(*data)->nhapl[0] != (int)(*data)->n_loci) { 
			printf("Error: Once 'rmfix' is defined, 'nhapl' must be defined for each locus.\n");
			exit(1);
		}
	}
    if((*data)->rmfix) {
        if((*data)->range_rnt != 0 && (*data)->range_rnt != 1 && (*data)->range_rnt != 2) {
            printf("Error: range_recnt must be 0, 1 or 2\n");
            exit(1);
        }
        else {
            if((*data)->range_rnt == 1 || (*data)->range_rnt == 2) {
                if((*data)->recnt_min[0] == (*data)->recnt_max[0] && (*data)->recnt_min[0] > 0) {
					for(x=1;x<=(*data)->recnt_min[0];x++) {
						if((*data)->recnt_min[x] < 0) {
							printf("Error: recnt_min must be positive or a zero value\n");
							exit(1);
						}
						if((*data)->recnt_max[x] < (*data)->recnt_min[x]) {
							printf("Error: recnt_max must be positive or a zero value, and lower or equal than recnt_min\n\n");
							exit(1);
						}
					}
				}
				else {
					printf("Error: recnt_max and recnt_max must be defined\n");
					exit(1);
				}
			}
			if((*data)->range_rnt == 0) {
				if((*data)->ifgammar == 0) {
                    printf("Error in input: options rmfix, range_rant and ifgammar. \n");
                    exit(1);
				}
			}
        }
    }

	if((*data)->ifgammar == 1) {
		if((*data)->p_gammar[0] != (double)(*data)->n_loci) {
			printf("Error: p gammar parameter must be defined for ALL loci and be higher than 0.\n");
			exit(1);
		}
		if((*data)->alpha_gammar[0] != (double)(*data)->n_loci) {
			printf("Error: alpha gammar parameter must be defined for ALL loci and be higher than 0.\n");
			exit(1);
		}
		if((*data)->correct_gammar[0] == (double)0) {
			if(!((*data)->correct_gammar = (double *)realloc((*data)->correct_gammar,(unsigned)((*data)->n_loci+1)*sizeof(double)))) 
				perror("realloc error.155");
			(*data)->correct_gammar[0] = (double)(*data)->n_loci;
			for(x=1;x<=(*data)->n_loci;x++) {
				(*data)->correct_gammar[x] = (double)1.0;
			}
		}
		if((*data)->correct_gammar[0] != (double)(*data)->n_loci) {
			printf("Error: correction gammar parameter must be defined for ALL loci and be higher than 0.\n");
			exit(1);
		}
	}
	for(x=1;x<=(*data)->n_loci;x++) {
		if(((*data)->ifgammar == 1 && (*data)->rmfix > 0) && ((*data)->p_gammar[x] <= (double)0.0 || (*data)->alpha_gammar[x] <= (double)0. || (*data)->correct_gammar[x] <= (double)0.)) {
			printf("Error: gammar parameters must be defined for ALL loci and be higher than 0.\n");
			exit(1);
		}
	}

	if((*data)->p_gammar[0] == (double)0) {
		if(!((*data)->p_gammar = (double *)realloc((*data)->p_gammar,(unsigned)((*data)->n_loci+1)*sizeof(double)))) 
			perror("realloc error.158");
		(*data)->p_gammar[0] = (double)(*data)->n_loci;
		for(x=1;x<=(*data)->n_loci;x++) {
			(*data)->p_gammar[x] = (double)1.0;
		}
	}
	if((*data)->p_gammar[0] < (*data)->n_loci) {
		if(!((*data)->p_gammar = (double *)realloc((*data)->p_gammar,(unsigned)((*data)->n_loci+1)*sizeof(double)))) 
			perror("realloc error.159");
		for(x=(int)(*data)->p_gammar[0];x<=(int)(*data)->n_loci;x++) {
			(*data)->p_gammar[x] = (double)1.0;
		}
		(*data)->p_gammar[0] = (double)(*data)->n_loci;
	}

	if((*data)->alpha_gammar[0] == (double)0) {
		if(!((*data)->alpha_gammar = (double *)realloc((*data)->alpha_gammar,(unsigned)((*data)->n_loci+1)*sizeof(double)))) 
			perror("realloc error.159n");
		(*data)->alpha_gammar[0] = (double)(*data)->n_loci;
		for(x=1;x<=(*data)->n_loci;x++) {
			(*data)->alpha_gammar[x] = (double)1.0;
		}
	}
	if((*data)->alpha_gammar[0] < (*data)->n_loci) {
		if(!((*data)->alpha_gammar = (double *)realloc((*data)->alpha_gammar,(unsigned)((*data)->n_loci+1)*sizeof(double)))) 
			perror("realloc error.159p");
		for(x=(int)(*data)->alpha_gammar[0];x<=(int)(*data)->n_loci;x++) {
			(*data)->alpha_gammar[x] = (double)1.0;
		}
		(*data)->alpha_gammar[0] = (double)(*data)->n_loci;
	}
	
	if(((*data)->sfix_allthetas == 1 || (*data)->rmfix == 1) && (*data)->method_samp == 0) {
		printf("Error: when sfix_allthetas or rmfix is defined, method_samp must be 1 or 2.\n");
		exit(1);
	}

	if(!((*data)->no_rec_males == 1 || (*data)->no_rec_males == 0)) {
		printf("Error: no_rec_males must be 1 or 0.\n");
		exit(1);
	}

	/*observed data*/
    if((*data)->despl > 0 && (*data)->window > 0) {
        totalnloci = (int)ceil(((double)(*data)->nsites[1] - (double)(*data)->window) / ((double)(*data)->despl)) + (int)1;
	}	
    else {
        if((*data)->linked > 1) totalnloci = (*data)->linked;
        else totalnloci = (*data)->n_loci;
	}

	if((*data)->likelihood_line) {
		y = 0;
		for(x=0;x<NOBS_STATISTICS;x++) {
			if((*data)->obs_statistics[x][1] == (double)1) {
				y = 1;
				break;
			}
		}
		if(y==0) {
			printf("Error: At least one observed value must be defined when likelihood_line is active.\n");
			exit(1);
		}
	}
	
	for(x=0;x<NOBS_STATISTICS;x++) {
		if((*data)->obs_statistics[x][0] == (double)0) {
			if(!((*data)->obs_statistics[x] = (double *)realloc((*data)->obs_statistics[x],(unsigned)(totalnloci+2)*sizeof(double)))) 
				perror("realloc error.obs_statistics");
			(*data)->obs_statistics[x][0] = (double)totalnloci+(double)1;
			(*data)->obs_statistics[x][1] = (double)0;
		}
		else {
			if((*data)->obs_statistics[x][1] == (double)1 && (*data)->obs_statistics[x][0] != (double)totalnloci+(double)1) {
				switch(x) {
					case 0:
						printf("Error: TD observed values must have nloci values when is defined\n");
						exit(1);
						break;
					case 1:
						printf("Error: Fs observed values must have nloci values when is defined\n");
						exit(1);
						break;
					case 2:
						printf("Error: FDn observed values must have nloci values when is defined\n");
						exit(1);
						break;
					case 3:
						printf("Error: FFn observed values must have nloci values when is defined\n");
						exit(1);
						break;
					case 4:
						printf("Error: FD observed values must have nloci values when is defined\n");
						exit(1);
						break;
					case 5:
						printf("Error: FF observed values must have nloci values when is defined\n");
						exit(1);
						break;
					case 6:
						printf("Error: H observed values must have nloci values when is defined\n");
						exit(1);
						break;
					case 7:
						printf("Error: B observed values must have nloci values when is defined\n");
						exit(1);
						break;
					case 8:
						printf("Error: Q observed values must have nloci values when is defined\n");
						exit(1);
						break;
					case 9:
						printf("Error: ZA observed values must have nloci values when is defined\n");
						exit(1);
						break;
					case 10:
						printf("Error: Fst observed values must have nloci values when is defined\n");
						exit(1);
						break;
					case 11:
						printf("Error: Kw observed values must have nloci values when is defined\n");
						exit(1);
						break;
					case 12:
						printf("Error: Hw observed values must have nloci values when is defined\n");
						exit(1);
						break;
					case 13:
						printf("Error: R2 observed values must have nloci values when is defined\n");
						exit(1);
						break;
					case 14:
						printf("Error: Ssites observed values must have nloci values when is defined\n");
						exit(1);
						break;
					case 15:
						printf("Error: piw observed values must have nloci values when is defined\n");
						exit(1);
						break;
					case 16:
						printf("Error: pib observed values must have nloci values when is defined\n");
						exit(1);
						break;
					case 17:
						printf("Error: thetaWatt observed values must have nloci values when is defined\n");
						exit(1);
						break;
					case 18:
						printf("Error: thetaTaj observed values must have nloci values when is defined\n");
						exit(1);
						break;
					case 19:
						printf("Error: thetaFW observed values must have nloci values when is defined\n");
						exit(1);
						break;
					case 20:
						printf("Error: D_Dmin observed values must have nloci values when is defined\n");
						exit(1);
						break;
					case 21:
						printf("Error: H_norm observed values must have nloci values when is defined\n");
						exit(1);
						break;
					case 22:
						printf("Error: maxhap observed values must have nloci values when is defined\n");
						exit(1);
						break;
					case 23:
						printf("Error: maxhap1 observed values must have nloci values when is defined\n");
						exit(1);
						break;
					case 24:
						printf("Error: Rm observed values must have nloci values when is defined\n");
						exit(1);
						break;
					case 25:
						printf("Error: thetafl observed values must have nloci values when is defined\n");
						exit(1);
						break;
					case 26:
						printf("Error: thetal observed values must have nloci values when is defined\n");
						exit(1);
						break;
					case 27:
						printf("Error: zengE observed values must have nloci values when is defined\n");
						exit(1);
						break;
					case 28:
						printf("Error: EW observed values must have nloci values when is defined\n");
						exit(1);
						break;
					case 29:
						printf("Error: Fstw observed values must have nloci values when is defined\n");
						exit(1);
						break;
					case 30:
						printf("Error: Pwh observed values must have nloci values when is defined\n");
						exit(1);
						break;
					default:
						break;
				}
			}
			else {
				if((*data)->obs_statistics[x][1] == (double)0) {
					if(!((*data)->obs_statistics[x] = (double *)realloc((*data)->obs_statistics[x],(unsigned)(totalnloci+2)*sizeof(double)))) 
						perror("realloc error.obs_statistics");
					(*data)->obs_statistics[x][0] = (double)totalnloci+(double)1;
					(*data)->obs_statistics[x][1] = (double)0;
				}
				else {
					if((*data)->obs_statistics[x][1] != (double)0 && (*data)->obs_statistics[x][1] != (double)1) {
						switch(x) {
							case 0:
								printf("Error: TD observed values must be 0 in the first value when is not defined and 1 when defined\n");
								exit(1);
								break;
							case 1:
								printf("Error: Fs observed values must have nloci values when is defined\n");
								exit(1);
								break;
							case 2:
								printf("Error: FDn observed values must have nloci values when is defined\n");
								exit(1);
								break;
							case 3:
								printf("Error: FFn observed values must have nloci values when is defined\n");
								exit(1);
								break;
							case 4:
								printf("Error: FD observed values must have nloci values when is defined\n");
								exit(1);
								break;
							case 5:
								printf("Error: FF observed values must have nloci values when is defined\n");
								exit(1);
								break;
							case 6:
								printf("Error: H observed values must have nloci values when is defined\n");
								exit(1);
								break;
							case 7:
								printf("Error: B observed values must have nloci values when is defined\n");
								exit(1);
								break;
							case 8:
								printf("Error: Q observed values must have nloci values when is defined\n");
								exit(1);
								break;
							case 9:
								printf("Error: ZA observed values must have nloci values when is defined\n");
								exit(1);
								break;
							case 10:
								printf("Error: Fst observed values must have nloci values when is defined\n");
								exit(1);
								break;
							case 11:
								printf("Error: Kw observed values must have nloci values when is defined\n");
								exit(1);
								break;
							case 12:
								printf("Error: Hw observed values must have nloci values when is defined\n");
								exit(1);
								break;
							case 13:
								printf("Error: R2 observed values must have nloci values when is defined\n");
								exit(1);
								break;
							case 14:
								printf("Error: Ssites observed values must have nloci values when is defined\n");
								exit(1);
								break;
							case 15:
								printf("Error: piw observed values must have nloci values when is defined\n");
								exit(1);
								break;
							case 16:
								printf("Error: pib observed values must have nloci values when is defined\n");
								exit(1);
								break;
							case 17:
								printf("Error: thetaWatt observed values must have nloci values when is defined\n");
								exit(1);
								break;
							case 18:
								printf("Error: thetaTaj observed values must have nloci values when is defined\n");
								exit(1);
								break;
							case 19:
								printf("Error: thetaFW observed values must have nloci values when is defined\n");
								exit(1);
								break;
							case 20:
								printf("Error: D_Dmin observed values must have nloci values when is defined\n");
								exit(1);
								break;
							case 21:
								printf("Error: H_norm observed values must have nloci values when is defined\n");
								exit(1);
								break;
							case 22:
								printf("Error: maxhap observed values must have nloci values when is defined\n");
								exit(1);
								break;
							case 23:
								printf("Error: maxhap1 observed values must have nloci values when is define\nd");
								exit(1);
								break;
							case 24:
								printf("Error: Rm observed values must have nloci values when is defined\n");
								exit(1);
								break;
							case 25:
								printf("Error: thetafl observed values must have nloci values when is defined\n");
								exit(1);
								break;
							case 26:
								printf("Error: thetal observed values must have nloci values when is defined\n");
								exit(1);
								break;
							case 27:
								printf("Error: zengE observed values must have nloci values when is defined\n");
								exit(1);
								break;
							case 28:
								printf("Error: EW observed values must have nloci values when is defined\n");
								exit(1);
								break;
							case 29:
								printf("Error: Fstw observed values must have nloci values when is defined\n");
								exit(1);
								break;
							case 30:
								printf("Error: Pwh observed values must have nloci values when is defined\n");
								exit(1);
								break;
							default:
								break;
						}
					}
				}
			}
		}
	}
}

void output_data(FILE *out, char *in, struct var **data)
{
    void print_var(int,FILE *);

    int x,y,z,zz;
    double j;
    time_t now;
    struct tm *date;
    char s[80];
    
    time(&now);
    date = localtime(&now);
    strftime(s,80,"%c",date);
    
    fprintf(out,"\nOUTPUT FILE: date %s \n\nInput data from the file: %s\n",s,in);
    /*fputs("\nFinite island model. Population parameters in function of 4No\n",out);*/
    /*fputs("except for the parameters that indicate the size of the populations (in No).\n",out);*/
    print_var(23,out);
    fprintf(out," %d",(*data)->print_neuttest);
	print_var(16,out);
    fprintf(out," %d",(*data)->pr_matrix);

	/*
	print_var(100,out);
    fprintf(out,"%d",(*data)->likelihood_line);
	*/
    if((*data)->mhits && ((*data)->theta_1[1] > 0. || (*data)->sfix_allthetas || (*data)->range_thetant || (*data)->ifgamma)) {
        print_var(0,out);
        fprintf(out," %d",(*data)->mhits);
        print_var(36,out);
        fprintf(out," %f",(*data)->T_out);

        print_var(17,out);
        j = (*data)->ratio_sv[1];
        for(y=2;y<(*data)->n_loci+1;y++) if((*data)->ratio_sv[y] != j) break;
        if(y < (*data)->n_loci + 1) for(z=1;z<(*data)->n_loci+1;z++) fprintf(out," %.3G",(*data)->ratio_sv[z]);
        else fprintf(out," %.3G",j);
        print_var(18,out);
		fprintf(out," %ld",(*data)->seed2);
    }
    print_var(1,out);
    fprintf(out," %ld",(*data)->n_iter);                
    print_var(2,out);        
    fprintf(out," %ld",(*data)->seed1);
    fputs("\n",out);
    print_var(3,out);
    fprintf(out," %d",(*data)->n_loci);
    
    if((*data)->linked) {
        print_var(14,out);
        fprintf(out," %d",(*data)->linked);
        if((*data)->linked > 1) {
            print_var(15,out);
            for(x=0;x<(*data)->linked;x++) {
                if(x) fputs(",",out);
                for(z=1;z<3;z++) fprintf(out," %ld",(*data)->loci_linked[x][z]); 
            }
			/**/
			if((*data)->linked_segsites[1] != -1) {
				print_var(92,out);
				for(x=1;x<=(*data)->linked;x++) {
					fprintf(out," %d",(*data)->linked_segsites[x]);
				}
			}
			if((*data)->linked_rm[1] != -1) {
				print_var(93,out);
				for(x=1;x<=(*data)->linked;x++) {
					fprintf(out," %d",(*data)->linked_rm[x]);
				}
			}
			if((*data)->linked_nhapl[1] != 0) {
				print_var(96,out);
				for(x=1;x<=(*data)->linked;x++) {
					fprintf(out," %d",(*data)->linked_nhapl[x]);
				}
			}
			/**/
        }
        else {
            print_var(37,out);
            fprintf(out," %ld",(*data)->despl);
            print_var(38,out);
            fprintf(out," %ld",(*data)->window);
        }
    }
                    
    print_var(55,out);
	for(y=1;y<(*data)->factor_chrn[0]+1;y++) fprintf(out," %f",(*data)->factor_chrn[y]);
    print_var(6,out);        
    for(y=1;y<(int)(*data)->nsites[0]+1;y++) fprintf(out," %ld",(*data)->nsites[y]);
    print_var(4,out);        
    for(y=1;y<(*data)->nsam[0]+1;y++) fprintf(out," %d",(*data)->nsam[y]);
    
	if((*data)->range_thetant) {
		print_var(33,out);
		fprintf(out," %d",(*data)->range_thetant);
		if((*data)->ifgamma == 0) {
			print_var(34,out);
			for(y=1;y<(*data)->thetant_min[0]+1;y++) fprintf(out," %f",(*data)->thetant_min[y]);
			print_var(35,out);
			for(y=1;y<(*data)->thetant_max[0]+1;y++) fprintf(out," %f",(*data)->thetant_max[y]);
		}
	}
	if((*data)->ifgamma) {
		print_var(52,out);
		fprintf(out," %d",(*data)->ifgamma);
		print_var(54,out);        
		for(y=1;y<(*data)->alpha_gamma[0]+1;y++) fprintf(out," %G",(*data)->alpha_gamma[y]);
		print_var(53,out);        
		for(y=1;y<(*data)->p_gamma[0]+1;y++) fprintf(out," %g",(*data)->p_gamma[y]);
		print_var(97,out);        
		for(y=1;y<(*data)->correct_gamma[0]+1;y++) fprintf(out," %g",(*data)->correct_gamma[y]);
	}
	if(!((*data)->mutations[0] == 0 || ((*data)->mutations[0] == 1 && (*data)->mutations[1] == -1))) {
        for(y=1;y<(int)(*data)->mutations[0]+1;y++) if((*data)->mutations[y] > -1) break;
		if(y<(int)(*data)->mutations[0]+1) {
			print_var(10,out);
			for(y=1;y<(int)(*data)->mutations[0]+1;y++) fprintf(out," %ld",(*data)->mutations[y]);
			if((*data)->sfix_allthetas && (*data)->method_samp) {
				print_var(26,out);
				fprintf(out," %d",(*data)->sfix_allthetas);
				print_var(58,out);
				fprintf(out," %d",(*data)->method_samp);
				if((*data)->method_samp > 1) {
					print_var(27,out);
					fprintf(out," %d",(*data)->mc_jump);
				}
			}
		}
    }	
    if(!((*data)->theta_1[0] == (double)0 || ((*data)->theta_1[0] == (double)1 && (*data)->theta_1[1] == (double)0))) {
        print_var(9,out);
        for(y=1;y<(*data)->theta_1[0]+1;y++) fprintf(out," %G",(*data)->theta_1[y]);
    }
	if((*data)->heter_theta_alphag[1] > (double)0) {
		print_var(94,out);
		for(y=1;y<(*data)->n_loci+1;y++) fprintf(out," %f",(*data)->heter_theta_alphag[y]);
	}
    print_var(65,out);
	fprintf(out," %d",(*data)->no_rec_males);
	if((*data)->invariable_mut_sites[1] > (double)0) {
		print_var(99,out);        
		for(y=1;y<(*data)->invariable_mut_sites[0]+1;y++) fprintf(out," %ld",(*data)->invariable_mut_sites[y]);
	}
	if((*data)->range_rnt == 1 || (*data)->range_rnt == 2) {
		print_var(59,out);
		fprintf(out," %d",(*data)->range_rnt);
		print_var(60,out);
		for(y=1;y<(*data)->recnt_min[0]+1;y++) fprintf(out," %f",(*data)->recnt_min[y]);
		print_var(61,out);
		for(y=1;y<(*data)->recnt_max[0]+1;y++) fprintf(out," %f",(*data)->recnt_max[y]);
	}
	if((*data)->ifgammar == 1) {
		print_var(62,out);
		fprintf(out," %d",(*data)->ifgammar);
		print_var(63,out);        
		for(y=1;y<(*data)->alpha_gammar[0]+1;y++) fprintf(out," %G",(*data)->alpha_gammar[y]);
		print_var(64,out);        
		for(y=1;y<(*data)->p_gammar[0]+1;y++) fprintf(out," %g",(*data)->p_gammar[y]);
		print_var(98,out);        
		for(y=1;y<(*data)->correct_gammar[0]+1;y++) fprintf(out," %g",(*data)->correct_gammar[y]);
	}
	if((*data)->ifgammar == 0 && (*data)->range_rnt == 0) {
		print_var(5,out);        
		for(y=1;y<(*data)->r[0]+1;y++) fprintf(out," %G",(*data)->r[y]); 
    }
	if((*data)->rmfix && (*data)->method_samp) {
		print_var(57,out);
		for(y=1;y<(int)(*data)->Rm[0]+1;y++) fprintf(out," %d",(*data)->Rm[y]);
		print_var(66,out);
		for(y=1;y<(int)(*data)->nhapl[0]+1;y++) {
			if((*data)->nhapl[y] > 0)
				fprintf(out," %d",(*data)->nhapl[y]);
			else 
				fprintf(out," 0");
		}
		print_var(56,out);
		fprintf(out," %d",(*data)->rmfix);
		print_var(58,out);
		fprintf(out," %d",(*data)->method_samp);
		if((*data)->method_samp > 1) {
			print_var(27,out);
			fprintf(out," %d",(*data)->mc_jump);
		}
	}	
	if((*data)->heter_rm_alphag[1] > (double)0) {
		print_var(95,out);
		for(y=1;y<(*data)->n_loci+1;y++) fprintf(out," %f",(*data)->heter_rm_alphag[y]);
	}
	
    z = 0;
    for(y=1;y<(*data)->f[0];y++) {/*only print f and track_len if defined*/
        if((*data)->f[y]>0.) {
            z=1;
            break;
        }
    }
    if(z==1) {
        fputs("\n",out);
        if(!((*data)->f[0] == 0 || (*data)->f[1] == 0.)) {
            print_var(7,out);
            for(y=0;y<(*data)->n_loci;y++) { 
                fprintf(out," %G",(*data)->f[y]); 
            }
        }
        if(!((*data)->track_len[0] == 0 || (*data)->track_len[1] == 0.)) {
            print_var(8,out);        
            for(y=0;y<(*data)->n_loci;y++) { 
                fprintf(out," %G",(*data)->track_len[y]); 
            }
        }
    }
	
    z = 0;
	if((*data)->ifselection[0]>0) {
		for(y=1;y<=(*data)->ifselection[0];y++) {
            if((*data)->ifselection[y]==1) {
				z=1;
				break;
			}
        }
    }
    if(z==0) {
        if((*data)->mig_rate != 0. && (*data)->split_pop == 0) {
            fputs("\n",out);
            print_var(11,out);
            fprintf(out," %ld",(*data)->npop);
            print_var(13,out);
            fprintf(out," %G",(*data)->mig_rate);
            fputs("\n",out);
            print_var(24,out);
            /*for(y=0;y<(*data)->npop-1;y++)*/ fprintf(out," %d",(*data)->npop_sampled[/*y+*/1]);
            
            if((*data)->ran_factorpop == 0 && (*data)->same_factorpop == 0) {
                zz= 1;
                if((*data)->split_pop == 0) {
                    for(y=2;y<(*data)->npop_sampled[0]+1;y++) 
                        if((*data)->npop_sampled[y] > (*data)->npop_sampled[zz]) zz = y;
                    if(!((*data)->factor_pop[0] == 0 || ((*data)->factor_pop[0] == 1 && (*data)->factor_pop[1] == 0.0))){
                        print_var(19,out);
                        for(y=1;y<(*data)->factor_pop[0]+1;y++) fprintf(out," %.3G",(*data)->factor_pop[y]); 
                    }
                }
                else for(y=1;y<(*data)->npop_sampled[zz]+1;y++) fprintf(out," %.3G",(*data)->factor_pop[y]);
            }
            print_var(20,out);
            fprintf(out," %d",(*data)->ran_factorpop);
            print_var(21,out);
            fprintf(out," %d",(*data)->same_factorpop);
            if(!((*data)->ssize_pop[0][0] == 0 || ((*data)->ssize_pop[0][0] == 1 && (*data)->ssize_pop[0][1] == 0))){
                print_var(12,out);
                for(y=0;y<(*data)->n_loci;y++) { 
                    /*per˜ sembla ho tinc restringit a nomes la mateixa mostra per tots locus...*/
                    if(y) fputs(",",out);
                    for(zz=1;zz<(*data)->ssize_pop[y][0]+1;zz++)
                        fprintf(out," %d",(*data)->ssize_pop[y][zz]); 
                }
            }
            fputs("\n",out);
        }
        if(!((*data)->nintn == 0)) {
            print_var(39,out);
            fprintf(out," %d",(*data)->nintn);
            print_var(50,out);
            fprintf(out," %d",(*data)->iflogistic);
            if((*data)->iflogistic) {
				print_var(51,out);
				fprintf(out," %G",(*data)->ts);
			}
            print_var(40,out);
            for(y=1;y<(*data)->nintn+1;y++) fprintf(out," %G",(*data)->nrec[y]);
            print_var(41,out);        
            for(y=1;y<(*data)->nintn+1;y++) fprintf(out," %G",(*data)->npast[y]);
            print_var(42,out);        
            for(y=1;y<(*data)->nintn+1;y++) fprintf(out," %G",(*data)->tpast[y]);
        }
        fputs("\n",out);
        if((*data)->split_pop == 1) {
            print_var(11,out);
            fprintf(out," %d",(*data)->npoprefugia);
            print_var(13,out);
            fprintf(out," %G",(*data)->mig_rate);
            print_var(43,out);
            fprintf(out," %d",(*data)->split_pop);
            print_var(44,out);
            fprintf(out," %f",(*data)->time_split);
            print_var(45,out);
            fprintf(out," %f",(*data)->time_scoal);
            print_var(46,out);
            fprintf(out," %f",(*data)->factor_anc);
            print_var(47,out);        
            for(y=1;y<(*data)->npoprefugia+1;y++) fprintf(out," %G",(*data)->freq[y]);
            print_var(19,out);        
            for(y=1;y<(*data)->npoprefugia+1;y++) fprintf(out," %G",(*data)->factor_pop[y]);
        }
    }
    else{ /*in case selection*/
        fputs("\n",out);
		print_var(29,out);
        fprintf(out," %.2E",(*data)->pop_size);
        print_var(28,out);
        for(y=1;y<(*data)->ifselection[0]+1;y++) fprintf(out," %d",(*data)->ifselection[y]);
        print_var(30,out);
        for(y=1;y<(*data)->pop_sel[0]+1;y++) fprintf(out," %G",(*data)->pop_sel[y]);
        print_var(31,out);
        for(y=1;y<(*data)->sel_nt[0]+1;y++) fprintf(out," %ld",(*data)->sel_nt[y]);
        print_var(32,out);
        for(y=1;y<(*data)->sinit[0]+1;y++) fprintf(out," %G",(*data)->sinit[y]);
        fputs("\n\nRecombination parameter value for the studied region and all the region until the selective position.\n",out);
    }
    if((*data)->sfix_allthetas == 1 && (*data)->method_samp == 1) 
        fputs("\nREJECTION ALGORITHM OF TAVARE et al. 1997. Fix S and screening all probable thetas.  Rejection algorithm n. 2.",out);
    if((*data)->sfix_allthetas == 1 && (*data)->method_samp == 2)
        fputs("\nMCMC ALGORITHM. Fix S and screening all probable thetas.",out);
    if((*data)->rmfix == 1 && (*data)->method_samp == 1) 
        fputs("\nREJECTION ALGORITHM. Fix Rm (and optionally nhapl) and screening all probable R.",out);
    if((*data)->rmfix == 1 && (*data)->method_samp == 2)
        fputs("\nMCMC ALGORITHM. Fix Rm (and optionally nhapl) and screening all probable R.",out);
    fputs("\n\n",out);
    
    if(fflush(out) != 0) {
        puts("\nError. Buffer print error\n");
        exit(1);
    }
}

void print_var(int x,FILE *out)
{
    int y,z;
    char til20[20];

    til20[0] = var_file[x][0] - 32;
    for(y=1;y<19;y++) {
        til20[y] = var_file[x][y];
        if(til20[y] == '\0') break;
    }
    til20[y] = ':';
    for(z=y+1;z<19;z++) til20[z] = 32;
    til20[19] = '\0';
    fprintf(out,"\n%s   ",til20);
    
    return;
}

/*************** FREE *********************/
void free_inputdata(struct var **data)
{
    int x,y;
    int loc;
    
    free((*data)->nsam);
    free((*data)->r);
    free((*data)->nsites);
    if((*data)->linked > 20) loc = (*data)->linked;
    else loc = 20;   
    for(x=0;x<loc;x++) free((*data)->loci_linked[x]);
    free((*data)->loci_linked);
    free((*data)->f);
    free((*data)->track_len);
    free((*data)->theta_1);    
    free((*data)->mutations);
    if((*data)->n_loci > 20) loc = (*data)->n_loci;
    else loc = 20;   
    for(y=0;y<loc;y++) free((*data)->ssize_pop[y]);
    free((*data)->ssize_pop);
    free((*data)->factor_pop);
    free((*data)->ratio_sv);
    free((*data)->npop_sampled);
    free((*data)->ifselection);
    free((*data)->pop_sel);
    free((*data)->sel_nt);
    free((*data)->sinit);
    free((*data)->nrec);
    free((*data)->npast);
    free((*data)->tpast);
    free((*data)->freq);
	
    free((*data)->p_gamma);
    free((*data)->alpha_gamma);
    free((*data)->correct_gamma);
    free((*data)->factor_chrn);

    free((*data)->p_gammar);
    free((*data)->alpha_gammar);
    free((*data)->correct_gammar);
    free((*data)->Rm);
    free((*data)->nhapl);
	
    free((*data)->linked_segsites);
    free((*data)->linked_rm);
    free((*data)->linked_nhapl);
    
	free((*data)->heter_theta_alphag);
    free((*data)->invariable_mut_sites);
    free((*data)->heter_rm_alphag);

    free((*data)->thetant_min);
    free((*data)->thetant_max);

    free((*data)->recnt_min);
    free((*data)->recnt_max);

	for(x=0;x<NOBS_STATISTICS;x++) free((*data)->obs_statistics[x]);

    return;
}
