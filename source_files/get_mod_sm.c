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

#include <stdlib.h>
#include "mlsp_sm.h"

/* fer una nova estructura que contingui totes les variables en proces */
void getpars_fix(struct var **data, struct var2 **inputp)
{    
    int x;
    double ran1(void);

    if((*data)->sfix_allthetas == 1 && (*data)->method_samp == 1) 
        (*inputp)->Sfix_alltheta = 1; /*Rejection algorithm*/
    if((*data)->sfix_allthetas == 1 && (*data)->method_samp == 2) {
        (*inputp)->Sfix_alltheta = 2; /*MHMCMC algorithm*/
        (*inputp)->mc_jump = (*data)->mc_jump;
    }
    
    if((*data)->rmfix == 1 && (*data)->method_samp == 1) 
        (*inputp)->rmfix = 1; /*Rejection algorithm*/
    if((*data)->rmfix == 1 && (*data)->method_samp == 2) {
        (*inputp)->rmfix = 2; /*MHMCMC algorithm*/
        (*inputp)->mc_jump = (*data)->mc_jump;
    }

    (*inputp)->howmany = (*data)->n_iter;			
    (*inputp)->pr_matrix = (*data)->pr_matrix;						
    if((*data)->split_pop == 1) (*inputp)->npop = (*data)->npoprefugia;
	else (*inputp)->npop = (*data)->npop;
    (*inputp)->config = (int *)malloc((*inputp)->npop*sizeof(int));
    (*inputp)->mhits = (*data)->mhits;					
    (*inputp)->T_out = (*data)->T_out;					
    (*inputp)->f = (double)0;
    (*inputp)->track_len = (double)0;
    (*inputp)->tloci = (*data)->n_loci;

    (*inputp)->linked = (*data)->linked;
    if((*inputp)->linked > 1) {
        (*inputp)->loci_linked = (long int **)malloc((*inputp)->linked*sizeof(long int *));		
        (*inputp)->linked_segsites = (int *)malloc((*inputp)->linked*sizeof(int));
        (*inputp)->linked_rm = (int *) malloc((*inputp)->linked*sizeof(int));		
        (*inputp)->linked_nhapl = (int *) malloc((*inputp)->linked*sizeof(int));		
        for(x=0;x<(*inputp)->linked;x++) {
            (*inputp)->loci_linked[x] = (long int *)malloc(2*sizeof(long int));
            (*inputp)->loci_linked[x][0] = (*data)->loci_linked[x][1];
            (*inputp)->loci_linked[x][1] = (*data)->loci_linked[x][2];
			(*inputp)->linked_segsites[x] = (*data)->linked_segsites[x+1];
			(*inputp)->linked_rm[x] = (*data)->linked_rm[x+1];
			(*inputp)->linked_nhapl[x] = (*data)->linked_nhapl[x+1];
        }
    }

	if((*inputp)->linked == 1) {
		(*inputp)->despl = (*data)->despl;
		(*inputp)->window = (*data)->window;
	}
	else {
		(*inputp)->despl = 0;
		(*inputp)->window = 0;
	}

    if((*inputp)->npop > 1 || (*data)->split_pop == 1 ) {
        (*inputp)->factor_pop = (double *) malloc(((*inputp)->npop)*sizeof(double));
        for(x=0;x<(*data)->factor_pop[0];x++) (*inputp)->factor_pop[x] = (*data)->factor_pop[x+1];
        (*inputp)->ran_factorpop = (*data)->ran_factorpop;
        (*inputp)->same_factorpop = (*data)->same_factorpop;
 
        if((*inputp)->ran_factorpop == 1) {
            (*inputp)->factor_pop[0] = (double)1;
            for(x=1;x<(*inputp)->npop;x++) {
                (*inputp)->factor_pop[x] = (double)ran1()*(double)9 + (double)1;
                if((double)ran1() < 0.5) (*inputp)->factor_pop[x] = (double)1/(*inputp)->factor_pop[x];
            }
        }
        if((*inputp)->same_factorpop == 1) {
            for(x=0;x<(*inputp)->npop;x++) {
                (*inputp)->factor_pop[x] = (double)1;
            }
        }        
    }
    else {
        (*inputp)->factor_pop = (double *)malloc(1*sizeof(double));
        (*inputp)->factor_pop[0] = 1.;
        (*inputp)->ran_factorpop = 0;
        (*inputp)->same_factorpop = 1;
    }
    
    (*inputp)->neutral_tests = (*data)->neutral_tests;
    (*inputp)->print_neuttest = (*data)->print_neuttest;

    (*inputp)->pop_size = (*data)->pop_size;/*selection*/
    /*
    (*inputp)->range_thetant = (*data)->range_thetant;
    (*inputp)->thetant_min   = (*data)->thetant_min;
    (*inputp)->thetant_max   = (*data)->thetant_max;
    
    (*inputp)->range_rnt = (*data)->range_rnt;
    (*inputp)->recnt_min   = (*data)->recnt_min;
    (*inputp)->recnt_max   = (*data)->recnt_max;
    */
    (*inputp)->iflogistic = (*data)->iflogistic;
    (*inputp)->ts = (*data)->ts;
    (*inputp)->nintn = (*data)->nintn;
    (*inputp)->nrec = (double *)malloc((2+(*inputp)->nintn)*sizeof(double));
    (*inputp)->npast = (double *)malloc((2+(*inputp)->nintn)*sizeof(double));
    (*inputp)->tpast = (double *)malloc((2+(*inputp)->nintn)*sizeof(double));
    
    (*inputp)->nrec[1] = (double)1;
    (*inputp)->npast[1] = (double)1;
    (*inputp)->tpast[0] = (double)0;
    (*inputp)->tpast[1] = (double)1;

    for(x=1;x<=(*inputp)->nintn;x++) {
        (*inputp)->nrec[x] = (*data)->nrec[x];    
        (*inputp)->npast[x] = (*data)->npast[x];    
        (*inputp)->tpast[x] = (*data)->tpast[x] + (*inputp)->tpast[x-1];    
    }
    
    (*inputp)->split_pop = (*data)->split_pop;
    (*inputp)->time_split = (*data)->time_split;
    (*inputp)->time_scoal = (*data)->time_scoal;
    (*inputp)->factor_anc = (*data)->factor_anc;
    if((*inputp)->split_pop) {
        (*inputp)->freq = (double *)malloc((*inputp)->npop*sizeof(double));
        for(x=0;x<(*inputp)->npop;x++) {
            (*inputp)->freq[x] = (*data)->freq[x+1];
        }
    }
    (*inputp)->tlimit = (*data)->tlimit;    
    (*inputp)->no_rec_males = (*data)->no_rec_males;

    return;
}    

void getpars_mod(struct var **data, struct var2 **inputp,int p1)
{
    int x;
    
    (*inputp)->nloci = p1;
    (*inputp)->nsam = (*data)->nsam[p1+1];
    /*selection*/
    (*inputp)->ifselection = (*data)->ifselection[p1+1];
    if((*data)->ifselection[p1+1] == 1) {
        (*inputp)->npop = 2;
        (*inputp)->config = (int *)realloc((*inputp)->config,2*sizeof(int));
        (*inputp)->pop_sel = (*data)->pop_sel[p1+1];
        (*inputp)->sel_nt = (*data)->sel_nt[p1+1];
        (*inputp)->sinit = (*data)->sinit[p1+1];
    }
    if((*data)->ifselection[p1+1] == 0 && (p1 > 0 && (*data)->ifselection[p1] == 1)) {
        (*inputp)->npop = 1;
	}

    for(x=0;x<(*data)->ssize_pop[p1][0];x++)
        (*inputp)->config[x] = (*data)->ssize_pop[p1][x+1];
    for(x=(*data)->ssize_pop[p1][0];x<(*inputp)->npop;x++)
        (*inputp)->config[x] = 0;
            
    (*inputp)->nsites = (*data)->nsites[p1+1];    
    /*theta vs S*/
	if((*data)->theta_1[0]) {
		if(!((*data)->theta_1[0]==1 && (*data)->theta_1[1] == (double)0)) {
			(*inputp)->theta = (*data)->theta_1[p1+1];
			(*inputp)->segsitesin = -1;
		}
		else {
			(*inputp)->segsitesin = (*data)->mutations[p1+1];
			(*inputp)->theta = (double)0;
		}
	}
    else {
		(*inputp)->segsitesin = (*data)->mutations[p1+1];
		(*inputp)->theta = (double)0;
	}
	
	if((*data)->sfix_allthetas) {
		(*inputp)->segsitesin = (*data)->mutations[p1+1];
		(*inputp)->theta = (double)0;
	}
	
	if((*data)->rmfix) {
		(*inputp)->Rm = (*data)->Rm[p1+1];
		(*inputp)->nhapl = (*data)->nhapl[p1+1];
		(*inputp)->r = (double)-1;
	}
	else {
		(*inputp)->Rm = -1;
		(*inputp)->nhapl = -1;
		if((*data)->ifgammar==0 && (*data)->range_rnt==0) 
			(*inputp)->r = (*data)->r[p1+1];
		else 
			(*inputp)->r = (double)-1;
	}

	(*inputp)->ratio_sv = (*data)->ratio_sv[p1+1];
    (*inputp)->npop_sampled = (*data)->npop_sampled[p1+1];    
    
    (*inputp)->migrate = (*data)->mig_rate;    
    
    (*inputp)->range_thetant = (*data)->range_thetant;
    if((*data)->thetant_min[0] >= p1+1) x = p1+1;
	else x = 1;
	(*inputp)->thetant_min   = (*data)->thetant_min[x] * (double)(*data)->nsites[p1+1];
    (*inputp)->thetant_max   = (*data)->thetant_max[x] * (double)(*data)->nsites[p1+1];
	
    (*inputp)->range_rnt = (*data)->range_rnt;
    if((*data)->recnt_min[0] >= p1+1) x = p1+1;
	else x = 1;
    (*inputp)->recnt_min   = (*data)->recnt_min[x] * (double)(*data)->nsites[p1+1];
    (*inputp)->recnt_max   = (*data)->recnt_max[x] * (double)(*data)->nsites[p1+1];
	
	(*inputp)->ifgamma = (*data)->ifgamma;
	(*inputp)->p_gamma = (*data)->p_gamma[p1+1];
	(*inputp)->alpha_gamma = (*data)->alpha_gamma[p1+1];
	(*inputp)->correct_gamma = (*data)->correct_gamma[p1+1];
	
	(*inputp)->ifgammar = (*data)->ifgammar;
	(*inputp)->p_gammar = (*data)->p_gammar[p1+1];
	(*inputp)->alpha_gammar = (*data)->alpha_gammar[p1+1];
	(*inputp)->correct_gammar = (*data)->correct_gammar[p1+1];
	
	(*inputp)->factor_chrn = (double)1/(double)(*data)->factor_chrn[p1+1];
	
    (*inputp)->heter_theta_alphag = (*data)->heter_theta_alphag[p1+1];
    (*inputp)->invariable_mut_sites = (*data)->invariable_mut_sites[p1+1];
    (*inputp)->heter_rm_alphag = (*data)->heter_rm_alphag[p1+1];

    return;
}
/*************************************** FREE **************************************/
void free_getpars_fix(struct var **data, struct var2 **inputp) 
{
    int x;
    
    free((*inputp)->config);
    free((*inputp)->nrec);
    free((*inputp)->npast);
    free((*inputp)->tpast);
    
    if((*data)->linked > 1) {
        for(x=0;x<(*data)->linked;x++) 
            free((*inputp)->loci_linked[x]);
    }
    if((*inputp)->linked > 1) {
        free((*inputp)->loci_linked);
        free((*inputp)->linked_segsites);
        free((*inputp)->linked_rm);
        free((*inputp)->linked_nhapl);
	}
    free((*inputp)->factor_pop);
    
    if((*inputp)->split_pop)
        free((*inputp)->freq);

    return;
}

