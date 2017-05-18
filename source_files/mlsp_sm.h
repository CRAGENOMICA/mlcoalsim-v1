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

#define MLCOALSIM "mlcoalsim version 1.42 (20080318)"

#define SHOWPROGRESS 1
/*set to 0 when the user do not want to see the progress of simulations in stdout */
#define NOBS_STATISTICS 31
#define ADDITIONAL_PRINTS 12
/*In case sfixall_thetas == 1 or rmfix == 1*/
/*0: No additional prints. 1: print theta/rec observed. 2: print treelength. 3:print thetae estimated (for Sfix)*/
/* 12,13,23 and 123: combinations*/
#define MAXREJECT 1000000
/*maximum rejection if accepted are less than 1/MAXREJECT, abort program*/
#define TIMEDURATION 0

struct var
{
    int pr_matrix;
    int mhits;
    double T_out;

    long int  n_iter;
    long int seed1;
    long int seed2;
    int n_loci;
    double *ratio_sv;
    int linked;
    long int **loci_linked;
    long int despl;
    long int window;
    
    int *nsam;
    double *r;
    long int *nsites;        
    double *f;
    double *track_len;
    double *theta_1;
    long int *mutations;
    
    long int npop;
    int **ssize_pop;
    double mig_rate;
    
    int *ifselection;
    double *pop_sel;
    double *sinit;
    long int *sel_nt;
    double pop_size;

    double *factor_pop;
    int ran_factorpop;
    int same_factorpop;
    int *npop_sampled;

    int neutral_tests;
    int print_neuttest;
    int sfix_allthetas;
    int mc_jump;
    
    int range_thetant;
    double *thetant_min;
    double *thetant_max;
	
	int ifgamma;
	double *p_gamma;
	double *alpha_gamma;
	double *correct_gamma;

	int iflogistic;
	double ts;

    int nintn;
    double *nrec;
    double *npast;
    double *tpast;        

    int split_pop;
    double time_split;
    double time_scoal;
    double factor_anc;
    double *freq;
    
    double tlimit;
	int npoprefugia;
	
	double *factor_chrn;

    int rmfix;
	int method_samp;
	int *Rm;
	int *nhapl;

    int range_rnt;
    double *recnt_min;
    double *recnt_max;
	
	int ifgammar;
	double *p_gammar;
	double *alpha_gammar;
	double *correct_gammar;
	
	int no_rec_males;

	int *linked_segsites;
	int *linked_rm;
	int *linked_nhapl;
	
	double *heter_theta_alphag;
	long int *invariable_mut_sites;
	double *heter_rm_alphag;

	/*observed statistics*/
	double *obs_statistics[NOBS_STATISTICS];

	int likelihood_line;
	double likelihood_error[NOBS_STATISTICS];
};

struct var2
{
    int rmfix;
	int Rm;
	int nhapl;

    int range_rnt;
    double recnt_min;
    double recnt_max;
	
	int ifgammar;
	double p_gammar;
	double alpha_gammar;
	double correct_gammar;
    
    int range_thetant;
    double thetant_min;
    double thetant_max;

	int ifgamma;
	double p_gamma;
	double alpha_gamma;
	double correct_gamma;

    int Sfix_alltheta;
    int mc_jump;
    
    long int howmany;
    int nsam;
    double r;
    double f;
    double track_len;
    long int nsites;
    double theta;
    long int segsitesin;
    long int npop;
    int *config;
    double migrate;
   
    int linked;
    long int **loci_linked;
    long int despl;
    long int window;
    
    int nloci;
    int tloci;
    int pr_matrix;
    int mhits;
    double ratio_sv;
    double T_out;
    int Sout;
     
    double *factor_pop;
    int ran_factorpop;
    int same_factorpop;
    int npop_sampled;

    int neutral_tests;
    int print_neuttest;
    
    int ifselection;
    double pop_size;
    double pop_sel;
    long int sel_nt;
    double sinit;

	int iflogistic;
	double ts;
	
    int nintn;
    double *nrec;
    double *npast;
    double *tpast;        

    int split_pop;
    double time_split;
    double time_scoal;
    double factor_anc;
    double *freq;
    
    double tlimit;

	double factor_chrn;
	int no_rec_males;

	int *linked_segsites;
	int *linked_rm;
	int *linked_nhapl;
	
	double heter_theta_alphag;
	long int invariable_mut_sites;
	double heter_rm_alphag;

};

