/************************************************************************************/
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

/*HUDSON ms.c file modified*/
/*ALSO INCLUDED A FUNCTION FOR RM FROM J. Wall.*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "mlsp_sm.h"

#define PRIOR_THETA 0
/*1: Prior from the avg of empirical theta. 0: We use a uniform or a gamma distribution*/
#define NIT_PRIOR (int)500
#define FACTOR (double)1.5
#define INTERVALS 2*(int)sqrt((double)NIT_PRIOR)
/*please, NIT_PRIOR no more than 32000*/
#define SITESINC 100
#define PRINTTHETAS 0
#define NEUTVALUES2 31

#define DEBUGSOUT 0

long int maxsites = 1000;	/* la llargada de la regió */
struct node {
    int abv;
    int ndes;
    double time;
};
struct segl {
    long int beg;
    struct node *ptree;
    long int next;
};

struct dnapar {
    double k;
    long int S;
    int nhapl;
    int *fhapl;
    int B1;
    int Q1;
    int *freq;
    double piw;
    double pib;
    long int *unic;
    int maxhapl;
	int maxhapl1;
	int Rm;
	double thetaL;
	
	double withinw;
	double *pid;
};

struct prob_par {/*posterior probabilities for theta and rec*/
	double thetap;
	double recp;
	double Ttotp;
};
static struct prob_par **postp = NULL; /*posterior probabilities for theta and rec*/
static double **thetafromS = NULL; /*if the theta range is inside the simulated values*/
/*static double **rfromrm = NULL; *//*if the rec range is inside the simulated values*/

static double **matrix_test = NULL;
double **coef = NULL;
struct dnapar *neutpar = NULL;    
long int *posit;	/*important! externes perque al reallocar memoria no localitza la nova posicio. */
                        /*Tambe es fa posant la posicio de memoria*/
char **list;

int ms(struct var2 **inputp,FILE *output)
{
    long int count;
    long int segsites;
    int i;

    char **cmatrix(int,long int);
    long int gensam(long int,int,int *,long int,double,long int,double,double,double,double,int,
        long int,double *,double *,int,double,double,double,long int,double,int *,int,double *,double *,double *,
        int, double, double, double, double *, double,int,double,double,double,double *,double *);
    void print_matrix(long int, struct var2 **,FILE *,int,long int);
    void mod_mhits(long int, struct var2 **);
    double logp(double);
    double ran1(void);
    
    void init_coef(double *,int);
    void calc_neutpar(int,long int,struct var2 **,struct dnapar *,double);
    double tajima_d(double,int, double *);
    double Fs(int,double,int);
    double fl_d(int,int,int,double *);
    double fl_f(int,int,int,double,double *);
    double fl_d2(int,int,int,double *);
    double fl_f2(int,int,int,double,double *);
    double fay_wu(int,int *,double);
    double fay_wu_normalized2(int,double,double,double,double *,double);
    double Fst(double,double,int);
    double Zns(int,long int,struct var2 **);
    double Zns_window(struct var2 **,long int,long int);
    double ZnA_(int,long int,struct var2 **);
    double ZnA_window(struct var2 **,long int,long int);
    double testHap(int, int *);
    double koutg(int,int,int *,long int);
    double R2(long int *,double,int,int);
    double Gxi(int ,int *, double);
    double Gximod(int ,int *, double);
    double frabs(int ,int *, double);
    double tajima_dvsdmin(double,int, double *,int);
	double E_zeng(int,double,double,double,double *);
    double Fstw(double *, int *,double,int);
    double EWtest(int, int *);
	double pwh(int,double *);
    /*canviat per debug*/
    double estnm(int,int,int *,long int,char **);
    /*double gst1,gst2;*/
    double logPPoisson2(long int, double);
    
    double thetae,thetaemin,thetaemax,lengtht;
    int j,k,neutv,sep=0;
    long int jcount,kcount,mcount,*listnumbers;
    double u,logv;
	long int ui;
    double theta_canddt,theta_mc,logprob_mc,logprob_canddt,logPoissonkk;
    double div,divt,divr;
    int burn_in=0;
    double weight_thcanddt=1.,weight_thmc=1.;
    int compare_(const void *,const void *);
    double thprior[NIT_PRIOR+2];
    int windows,nwindow,aa;
    long int s0,s1,kk,ll;
	long int inputS;
	
	int Min_rec(int,int,int,int);
	int Rmi;
	/*double ZAi;*/int nhi;
	double rece,recemin,recemax,rec_mc,rec_canddt;
    double correction_rec(double,double,int);
	double recombinationv;
    void calc_neutpar_window(struct var2 **,struct dnapar *,long int,long int,double);

    double tdgamma[NIT_PRIOR];
	double gammadist(double);

    double lengtht_mc;    

	/*counting simulations in stdout*/
    static double counterp,restp;
    static long int p;
	int x;
	
	/*heterogeneity in mutation and recombination*/
	double *weightmut=0;
	double *weightrec=0;
	int do_heter_gamma_sites(double *,double,double,long int,long int);
	int do_heter_gamma_sitesrec(double *,double,double,long int);
	double thetaSv;
	
	long int nrej; /*avoid too much rejections*/
    
    #if PRIOR_THETA
    double sumslope;
    double thrange[INTERVALS+2];
    double thslope[INTERVALS+2];
    double prob0;
    int l;
    double ranprior,ranprior0;    
    #endif
	
	#if DEBUGSOUT
    FILE *file_debug;
    if (!(file_debug = fopen("debug.out","w"))) perror("Error in input/output");
	#endif
    #if PRINTTHETAS
    FILE *file_debug_sel;
    if (!(file_debug_sel = fopen("debug_sel.out","w"))) perror("Error in input/output");
    fputs("piTaj\tpiWatt\tpiFW\n",file_debug_sel);
    #endif
    
	/*heterogeneity in mutation and recombination*/
	/*init*/
	if(!(weightmut = (double *)calloc((*inputp)->nsites+1,sizeof(double)))) {
		perror("calloc error ms.89");
		exit(1);
	}
	if(!(weightrec = (double *)calloc((*inputp)->nsites,sizeof(double)))) {
		perror("calloc error ms.89");
		exit(1);
	}

	/*counting simulations in stdout*/
	if((*inputp)->Sfix_alltheta == 2 || (*inputp)->rmfix == 2) {
		if((*inputp)->nloci == 0) {
			counterp = (double)(((*inputp)->howmany+(*inputp)->mc_jump) * (*inputp)->tloci)/(double)50;
			p = 1;
			restp = (double)0;
		}
	}
	else {
		if((*inputp)->nloci == 0) {
			counterp = (double)((*inputp)->howmany * (*inputp)->tloci)/(double)50;
			p = 1;
			restp = (double)0;
		}
	}

    /*two important parameters for Sfixallthetas = 2 !!*/
    if((*inputp)->Sfix_alltheta == 2 || (*inputp)->rmfix == 2) {
        sep = (*inputp)->mc_jump; 	/*separation of the values we pick in the markov chain: 1 fastests...*/
        if((*inputp)->Sfix_alltheta == 2 || (*inputp)->rmfix == 2)
            div = 1.;
        else div = 1.;
        /*Fixed to 1 works well. Range of the chosen value in the uniform distribution: (0,1]*/
    }
    /*Specially in case MHMCMC/MCMCRA, we need to mix values to avoid high correlations. */
    /*If not, all values are consecutives*/
    if(!(listnumbers = (long int *)malloc((unsigned)((*inputp)->howmany+sep)*sizeof(long int)))) {
        perror("calloc error ms.theta_acc");
        exit(1);
    }
    for(jcount=0;jcount<(*inputp)->howmany+sep;jcount++) listnumbers[jcount] = jcount;
    /*end case MCMC/MCMCRA*/
	/*Define the matrix for posterior probabilities for theta and other values*/
    if((*inputp)->ifgamma == 1 || (*inputp)->range_thetant || (*inputp)->ifgammar == 1 || (*inputp)->range_rnt) {
		if(postp == NULL) {
			if(!(postp = (struct prob_par **)calloc((*inputp)->tloci,sizeof(struct prob_par *)))) {
				perror("calloc error ms.09");
				exit(1);
			}
			if((*inputp)->ifgamma == 1 || (*inputp)->range_thetant) {
				if(!(thetafromS = (double **)calloc((*inputp)->tloci,sizeof(double *)))) {
					perror("calloc error ms.19b");
					exit(1);
				}
			}
			/*
			if((*inputp)->ifgammar == 1 || (*inputp)->range_rnt) {
				if(!(rfromrm = (double **)calloc((*inputp)->tloci,sizeof(double *)))) {
					perror("calloc error ms.19b");
					exit(1);
				}
			}
			*/
			for(i=0;i<(*inputp)->tloci;i++) {
				if(!(postp[i] = (struct prob_par *)calloc(((*inputp)->howmany + sep + 1),sizeof(struct prob_par)))) {
					perror("calloc error ms.09");
					exit(1);
				}
				if((*inputp)->ifgamma == 1 || (*inputp)->range_thetant) {
					if(!(thetafromS[i] = (double *)calloc(NIT_PRIOR+2,sizeof(double)))) {
						perror("calloc error ms.19b");
						exit(1);
					}
				}
				/*
				if((*inputp)->ifgammar == 1 || (*inputp)->range_rnt == 1) {
					if(!(rfromrm[i] = (double *)calloc(NIT_PRIOR+2,sizeof(double)))) {
						perror("calloc error ms.19b");
						exit(1);
					}
				}
				*/
			}
		}
	}
    
    if((*inputp)->theta > 0.0 || (((*inputp)->ifgamma == 1 || (*inputp)->range_thetant) && (*inputp)->Sfix_alltheta == 0) || ((*inputp)->linked > 1 && (*inputp)->Sfix_alltheta) || ((*inputp)->segsitesin > -1 && (*inputp)->mhits == 1)) {        
        list = cmatrix((*inputp)->nsam,maxsites+1);
        if(!(posit = (long int *)malloc((long int)(maxsites*sizeof(long int)))))
            perror("ms error.1");
    }
    else {
        list = cmatrix((*inputp)->nsam,(long int)(*inputp)->segsitesin+1);
        if(!(posit = (long int *)malloc((long int)(((long int)(*inputp)->segsitesin+1)*sizeof(long int)))))
            perror("ms error.2");
    }
    /*if((*inputp)->neutral_tests) {*/
        /* matriu double test taj,fs,fd,ff,h,B,Q,ZnS,Fst,#hap,divhap... per tots loci */
        if(matrix_test == NULL) {
            if((*inputp)->linked == 0) {/*independent regions*/
                if((matrix_test = (double **)calloc((*inputp)->tloci*NEUTVALUES2,sizeof(double *))) == NULL) {
                    perror("calloc error ms.1b");
                    exit(1);
                }
                for(i=0;i<(*inputp)->tloci*NEUTVALUES2;i++) {
                    if((matrix_test[i] = (double *)calloc(((*inputp)->howmany + sep),sizeof(double))) == NULL) {
                        perror("calloc error ms.1c");
                        exit(1);
                    }
                }
            }
            else {
                if((*inputp)->linked == 1) /*sliding windows*/
                    windows = (int)ceil(((double)(*inputp)->nsites - (double)(*inputp)->window) / ((double)(*inputp)->despl)) + 1;
                else /*separated but linked regions*/
                    windows = (*inputp)->linked;
                /*define matrix for linked*/
                if((matrix_test = (double **)calloc(windows*NEUTVALUES2,sizeof(double *))) == NULL) {
                    perror("calloc error ms.1b");
                    exit(1);
                }
                for(i=0;i<windows*NEUTVALUES2;i++) {
                    if((matrix_test[i] = (double *)calloc(((*inputp)->howmany + sep),sizeof(double))) == NULL) {
                        perror("calloc error ms.1c");
                        exit(1);
                    }
                }
            }
            if((neutpar = (struct dnapar *)calloc(1,sizeof(struct dnapar))) == NULL) {
                perror("calloc error ms.1d1");
                exit(1);
            }
            if((coef = (double **)calloc(1,sizeof(double *))) == NULL) {
                perror("calloc error ms.1d1");
                exit(1);
            }
            for(i=0;i<1;i++) {
                if((coef[i] = (double *)calloc(12,sizeof(double))) == NULL) {
                    perror("calloc error ms.1d1");
                    exit(1);
                }
            }
        }
        init_coef(coef[0],(*inputp)->nsam);

        if((neutpar[0].freq = (int *)calloc((*inputp)->nsam,sizeof(int))) == NULL) {
            perror("calloc error ms.1d3");
            exit(1);
        }
        if((neutpar[0].fhapl = (int *)calloc((*inputp)->nsam,sizeof(int))) == NULL) {
            perror("calloc error ms.1d3");
            exit(1);
        }
        if((neutpar[0].unic = (long int *)calloc((*inputp)->nsam,sizeof(long int))) == NULL) {
            perror("calloc error ms.1d3");
            exit(1);
        }
        if((neutpar[0].pid = (double *)calloc(((*inputp)->nsam * ((*inputp)->nsam-1))/2,sizeof(double))) == NULL) {
            perror("calloc error ms.1d6");
            exit(1);
        }
	/*}*/

	if((*inputp)->range_thetant) {
		thetaemax = (*inputp)->thetant_max;
		thetaemin = (*inputp)->thetant_min;
	}
	if((*inputp)->range_rnt) {
		recemax = (*inputp)->recnt_max;
		recemin = (*inputp)->recnt_min;
	}

    if((*inputp)->Sfix_alltheta && (*inputp)->rmfix == 0 /**/&& (*inputp)->linked < 2/**/ && (*inputp)->mhits == 0) { /*Calculate approx. min and max for theta. Do NIT_PRIOR iterations and estimate theta...*/
        count = 0;
        for(j=0;j<NIT_PRIOR;j++) {
            /*do a tree*/
			if((*inputp)->ifgammar == 1) rece = gammadist((*inputp)->alpha_gammar)/(*inputp)->p_gammar * (*inputp)->correct_gammar;
			else {
				if((*inputp)->range_rnt == 1) rece = recemin + (recemax - recemin) * ran1();
				else {
					if((*inputp)->range_rnt == 2) rece = (double)exp((double)log((double)recemin) + ((double)log((double)recemax) - (double)log((double)recemin)) * ran1());
					else rece = (double)(*inputp)->r;
				}
			}			
			recombinationv = correction_rec(rece,(*inputp)->factor_chrn,(*inputp)->no_rec_males);
			thetaSv = (double)1.0;
			do_heter_gamma_sites(weightmut,(*inputp)->heter_theta_alphag,thetaSv+thetaSv/(double)(*inputp)->nsites,(*inputp)->nsites+1,(*inputp)->invariable_mut_sites);
			do_heter_gamma_sitesrec(weightrec,(*inputp)->heter_rm_alphag,(double)recombinationv/*+recombinationv/(double)(*inputp)->nsites*/,(*inputp)->nsites);
            segsites = gensam((*inputp)->npop,(*inputp)->nsam, (*inputp)->config, (*inputp)->nsites,
            (*inputp)->theta, (*inputp)->segsitesin, recombinationv, (*inputp)->f, (*inputp)->track_len,
            (*inputp)->migrate,(*inputp)->mhits,count,(*inputp)->factor_pop, 
            &lengtht,(*inputp)->ifselection, 
            (*inputp)->pop_sel,(*inputp)->sinit,(*inputp)->pop_size,(*inputp)->sel_nt,(*inputp)->T_out,&(*inputp)->Sout,
            (*inputp)->nintn, (*inputp)->nrec, (*inputp)->npast, (*inputp)->tpast,
            (*inputp)->split_pop,(*inputp)->time_split,(*inputp)->time_scoal,(*inputp)->factor_anc,(*inputp)->freq,
            (*inputp)->tlimit,(*inputp)->iflogistic,(*inputp)->ts,(*inputp)->factor_chrn,(*inputp)->ratio_sv,
			weightmut,weightrec);
			
			/*WE INCLUDED LENGTHT!!!*/
            /*estimate theta given lengtht. Approach by its average....*/
			if(segsites==0) thetae = (double)0.1/lengtht;
			else thetae = (double)segsites/lengtht;
			/*keep thetae*/
            thetafromS[(*inputp)->nloci][j] = thetae;

            thprior[j+1] = thetae;
            /*we do two matrix (thprior and ttprior) in case we change the distribution of thetas for other dist.*/
            /*ttprior[j+1] = lengtht; */
            /*keep the theta max and the theta min*/
            if(j==0 || thetaemin > thetae) thetaemin = thetae; 
            if(j==0 || thetaemax < thetae) thetaemax = thetae;
        }
		
		if((*inputp)->ifgamma==1) {
			for(j=0;j<NIT_PRIOR;j++) {
				tdgamma[j] = gammadist((*inputp)->alpha_gamma)/(*inputp)->p_gamma * (*inputp)->correct_gamma;
			}
			qsort(tdgamma,(int)NIT_PRIOR,sizeof(double),compare_);/*all values are sorted*/
			thetaemax = tdgamma[(int)floor((double)((NIT_PRIOR+2) - (int)((NIT_PRIOR+2)/(double)40)))];
			thetaemin = tdgamma[(int)floor((double)((int)((NIT_PRIOR+2)/(double)40)))];
		}
		else {
			/*we'll use arbitrarially FACTOR times lower and bigger than the observed theta values*/
			thetaemin /= FACTOR;
			thetaemax *= FACTOR;
			/*WE CHECK THE 2.5% C.INTERVAL FOR EACH TAIL. OUT OF THESE C.I. THE THETA INTERVAL IS NOT COMPATIBLE*/
			/*for(j=0;j<NIT_PRIOR;j++) thprior[j+1] = thetaemin + ran1()*(thetaemax-thetaemin);*//*debug prior: uniform*/
			thprior[0] = thetaemin;
			thprior[NIT_PRIOR+1] = thetaemax;
			qsort(thprior,(int)NIT_PRIOR+2,sizeof(double),compare_);/*all values are sorted*/
			/*IN CASE UPPER AND LOWER BOUNDS OF THETA: WE USE THE NARROWER LIMITS. */
			/*THE EXPLANATION: LIMITS FROM THE SIMULATED DISTRIBUTION ARE LIKE INFINITE. */
			/*IF WE USE BIGGER VALUES, THE SIMULATIONS ARE VERY VERY SLOW*/
			/*WE CHECK THE 2.5% C.INTERVAL FOR EACH TAIL. OUT OF THESE C.I. THE THETA INTERVAL IS NOT COMPATIBLE*/
			if((*inputp)->range_thetant) {
				/*thetaemin = thprior[(int)((NIT_PRIOR+2)/40.0)];*//*use 95% of the distribution of theta*/
				/*thetaemax = thprior[(NIT_PRIOR+2) - (int)((NIT_PRIOR+2)/40.0)];*//*use 95% of the distribution of theta*/
				if(thetaemin < (*inputp)->thetant_min) thetaemin = (*inputp)->thetant_min;
				if(thetaemax > (*inputp)->thetant_max) thetaemax = (*inputp)->thetant_max;
			}
		}
		if(thetaemin > thetaemax) {
			puts("\nSorry. Theta values are not compatible with the number of segregating sites indicated.\n");
			fprintf(output,"\nSorry. Theta values are not compatible with the number of segregating sites indicated in locus[%d].\n",(*inputp)->nloci);
			return(1);
		}
		if(thetaemin > thprior[NIT_PRIOR] || thetaemax < thprior[1]) {
			printf("\nSorry. Theta values are not compatible with the number of segregating sites indicated in locus[%d].\n",(*inputp)->nloci);
			fprintf(output,"\nSorry. Theta values are not compatible with the number of segregating sites indicated in locus[%d].\n",(*inputp)->nloci);
			return(1);
		}
	}
    if(((*inputp)->Sfix_alltheta || (*inputp)->rmfix) && (*inputp)->linked < 2) { /*Calculate approx. min and max for theta. Do NIT_PRIOR iterations and estimate theta...*/
		if(!((*inputp)->Sfix_alltheta && (*inputp)->rmfix == 0)) {
			thetaemax = (*inputp)->thetant_max;
			thetaemin = (*inputp)->thetant_min;
		}
		if((*inputp)->rmfix && (*inputp)->linked < 2) {
			recemax = (*inputp)->recnt_max;
			recemin = (*inputp)->recnt_min;
		}

		if((*inputp)->Sfix_alltheta==2 || (*inputp)->rmfix==2) {
			divt = div*(thetaemax-thetaemin); /*range of the chosen value in the uniform distribution for theta*/
			divr = div*(recemax-recemin); /*range of the chosen value in the uniform distribution for rec*/
			if((*inputp)->range_thetant == 2) divt = div*((double)log(thetaemax)-(double)log(thetaemin)); /*range of the chosen value in the uniform distribution for theta*/
			if((*inputp)->range_rnt == 2) divr = div*((double)log(recemax)-(double)log(recemin)); /*range of the chosen value in the uniform distribution for rec*/
			burn_in = 0;

			/*MHMCMC we need to mix values to avoid high correlations*/
			for(jcount=0;jcount<(*inputp)->howmany-1;jcount++) {
				kcount = (long int)(ran1()*(double)((*inputp)->howmany-jcount)) + jcount;
				mcount = listnumbers[kcount];
				listnumbers[kcount] = listnumbers[jcount];
				listnumbers[jcount] = mcount;
			}   
		}
        /*other parameters*/        
        if((*inputp)->Sfix_alltheta) {
			if((*inputp)->segsitesin > 0) logPoissonkk = (double)logPPoisson2((long int)(*inputp)->segsitesin,(double)(*inputp)->segsitesin);
			else logPoissonkk = (double)logPPoisson2((long int)(*inputp)->segsitesin,(double)0.1);
		}
	}
	if((*inputp)->Sfix_alltheta || (*inputp)->rmfix) {
        j = 0; 
        k = 0; 
        jcount = 0; 
        mcount = 0; 
        kcount = 0; 
		neutv = NEUTVALUES2;
    }
	if((*inputp)->Sfix_alltheta == 0 && (*inputp)->rmfix == 0  && (*inputp)->linked < 2) {
		if((*inputp)->range_thetant) {
			thetaemax = (*inputp)->thetant_max;
			thetaemin = (*inputp)->thetant_min;
		}
		if((*inputp)->range_rnt) {
			recemax = (*inputp)->recnt_max;
			recemin = (*inputp)->recnt_min;
		}
	}
    /************************************************  ITERATIONS: Routine from Hudson and modified ****************************************************/
    count = 0; 
    while(((*inputp)->howmany + sep) - count++) {    
        /* ranfactor: between 0.1 and 1, or between 1 and 10. Equally divided */
        if((*inputp)->ran_factorpop == 2) {
            for(i=1;i<(*inputp)->npop;i++) {
                (*inputp)->factor_pop[i] = (double)ran1()* (double)9 + (double)1;
                if((double)ran1() < 0.5) (*inputp)->factor_pop[i] = 1./(*inputp)->factor_pop[i];
            }
        }
        if((*inputp)->Sfix_alltheta == 0 && (*inputp)->rmfix == 0) {/*"normal" simulations: Fix theta (or fix S) and fix R*/
			if((*inputp)->ifgamma == 1) thetae = gammadist((*inputp)->alpha_gamma)/(*inputp)->p_gamma * (*inputp)->correct_gamma;
			else {
				if((*inputp)->range_thetant == 1) thetae = thetaemin + (thetaemax - thetaemin) * ran1();
				else {
					if((*inputp)->range_thetant == 2) thetae = (double)exp((double)log((double)thetaemin) + ((double)log((double)thetaemax) - (double)log((double)thetaemin)) * ran1());
					else thetae = (double)(*inputp)->theta;
				}
			}
			if((*inputp)->ifgammar == 1) rece = gammadist((*inputp)->alpha_gammar)/(*inputp)->p_gammar * (*inputp)->correct_gammar;
			else {
				if((*inputp)->range_rnt == 1) rece = recemin + (recemax - recemin) * ran1();
				else {
					if((*inputp)->range_rnt == 2) rece = (double)exp((double)log((double)recemin) + ((double)log((double)recemax) - (double)log((double)recemin)) * ran1());
					else rece = (double)(*inputp)->r;
				}
			}			
            recombinationv = correction_rec(rece,(*inputp)->factor_chrn,(*inputp)->no_rec_males);
			if((*inputp)->segsitesin == -1) thetaSv = thetae;
			else thetaSv = (double)1.0;
			do_heter_gamma_sites(weightmut,(*inputp)->heter_theta_alphag,thetaSv+thetaSv/(double)(*inputp)->nsites,(*inputp)->nsites+1,(*inputp)->invariable_mut_sites);
			do_heter_gamma_sitesrec(weightrec,(*inputp)->heter_rm_alphag,(double)recombinationv/*+recombinationv/(double)(*inputp)->nsites*/,(*inputp)->nsites);
			segsites = gensam((*inputp)->npop,(*inputp)->nsam, (*inputp)->config, (*inputp)->nsites,
            thetae, (*inputp)->segsitesin, recombinationv, (*inputp)->f, (*inputp)->track_len,
            (*inputp)->migrate,(*inputp)->mhits,count,(*inputp)->factor_pop, 
            &lengtht,(*inputp)->ifselection, 
            (*inputp)->pop_sel,(*inputp)->sinit,(*inputp)->pop_size,(*inputp)->sel_nt,(*inputp)->T_out,&(*inputp)->Sout,
            (*inputp)->nintn, (*inputp)->nrec, (*inputp)->npast, (*inputp)->tpast,
            (*inputp)->split_pop,(*inputp)->time_split,(*inputp)->time_scoal,(*inputp)->factor_anc,(*inputp)->freq,
            (*inputp)->tlimit,(*inputp)->iflogistic,(*inputp)->ts,(*inputp)->factor_chrn,(*inputp)->ratio_sv,
			weightmut,weightrec);
            if((*inputp)->mhits) mod_mhits(segsites,inputp); /******** mhits ******************/
			
			/*printf("%ld\n",segsites);*/
			if((*inputp)->ifgamma == 1 || (*inputp)->range_thetant) 
				postp[(*inputp)->nloci][listnumbers[count-1]].thetap = thetae;
			if((*inputp)->ifgammar == 1 || (*inputp)->range_rnt) 
				postp[(*inputp)->nloci][listnumbers[count-1]].recp = rece;
			if((*inputp)->ifgamma == 1 || (*inputp)->range_thetant || (*inputp)->ifgammar == 1 || (*inputp)->range_rnt)
				postp[(*inputp)->nloci][listnumbers[count-1]].Ttotp  = lengtht;

			/*printing a point every 2% of the total iterations*/
			#if SHOWPROGRESS == 1
			if((double)p+restp >= counterp) {
				restp += (double)p - counterp; 
				if((double)restp/(double)counterp > (double)1) {
					for(x=0;x<(int)floor(restp/counterp);x++) printf(".");
					restp -= (double)floor(restp/counterp) * counterp;
				}
				else printf(".");
				fflush(stdout);
				p = 1;
			}
			else p += 1;
			#endif
        }
        /*Do the rejection algorithm (RA) of Tavaré et al. 1997.*/
        else {
			if((*inputp)->Sfix_alltheta == 1 || (*inputp)->rmfix == 1) { 
				nrej = jcount;
				do {
					if((*inputp)->ifgamma == 1) thetae = gammadist((*inputp)->alpha_gamma)/(*inputp)->p_gamma * (*inputp)->correct_gamma;
					else {
						if((*inputp)->range_thetant == 1) thetae = thetaemin + (thetaemax - thetaemin) * ran1();
						else {
							if((*inputp)->range_thetant == 2) thetae = (double)exp((double)log((double)thetaemin) + ((double)log((double)thetaemax) - (double)log((double)thetaemin)) * ran1());
							else thetae = (double)(*inputp)->theta;
						}
					}
					if((*inputp)->ifgammar == 1) rece = gammadist((*inputp)->alpha_gammar)/(*inputp)->p_gammar * (*inputp)->correct_gammar;
					else {
						if((*inputp)->range_rnt == 1) rece = recemin + (recemax - recemin) * ran1();
						else {
							if((*inputp)->range_rnt == 2) rece = (double)exp((double)log((double)recemin) + ((double)log((double)recemax) - (double)log((double)recemin)) * ran1());
							else rece = (double)(*inputp)->r;
						}
					}			
					recombinationv = correction_rec(rece,(*inputp)->factor_chrn,(*inputp)->no_rec_males);
					if((*inputp)->segsitesin == -1 || (*inputp)->mhits == 1) thetaSv = thetae;
					else thetaSv = (double)1.0;
					/*be careful with mhits*/
					inputS = (*inputp)->segsitesin;
					if((*inputp)->mhits == 1) inputS = -1;

					do_heter_gamma_sites(weightmut,(*inputp)->heter_theta_alphag,thetaSv+thetaSv/(double)(*inputp)->nsites,(*inputp)->nsites+1,(*inputp)->invariable_mut_sites);
					do_heter_gamma_sitesrec(weightrec,(*inputp)->heter_rm_alphag,(double)recombinationv/*+recombinationv/(double)(*inputp)->nsites*/,(*inputp)->nsites);
					segsites = gensam((*inputp)->npop,(*inputp)->nsam, (*inputp)->config, (*inputp)->nsites,
					thetae, /*(*inputp)->segsitesin*/inputS, recombinationv, (*inputp)->f, (*inputp)->track_len,
					(*inputp)->migrate,(*inputp)->mhits,count,(*inputp)->factor_pop, 
					&lengtht,(*inputp)->ifselection, 
					(*inputp)->pop_sel,(*inputp)->sinit,(*inputp)->pop_size,(*inputp)->sel_nt,(*inputp)->T_out,&(*inputp)->Sout,
					(*inputp)->nintn, (*inputp)->nrec, (*inputp)->npast, (*inputp)->tpast,
					(*inputp)->split_pop,(*inputp)->time_split,(*inputp)->time_scoal,(*inputp)->factor_anc,(*inputp)->freq,
					(*inputp)->tlimit,(*inputp)->iflogistic,(*inputp)->ts,(*inputp)->factor_chrn,(*inputp)->ratio_sv,
					weightmut,weightrec);
					if((*inputp)->mhits) mod_mhits(segsites,inputp); /******** mhits ******************/
					
					if((*inputp)->linked > 1) {/*linked fragments with subset fixed values*/
						s0=s1=0;
						nwindow=0;
						for(aa=0;aa<(*inputp)->linked;aa++) {
							kk=(*inputp)->loci_linked[aa][0];
							ll=(*inputp)->loci_linked[aa][1]+1;
							while(s0 < segsites && posit[s0] < kk) s0++;
							while(s1 < segsites && posit[s1] < ll) s1++;
							/*calc_estadistics*/
							calc_neutpar_window(inputp,neutpar,s0,s1,(double)recombinationv);
							if((*inputp)->rmfix == 1) {
								Rmi = neutpar[0].Rm;
								nhi = neutpar[0].nhapl;
								if(Rmi != (*inputp)->linked_rm[aa] || ((*inputp)->linked_nhapl[aa] != 0 && nhi != (*inputp)->linked_nhapl[aa])) {
									jcount += 1;
									if(jcount - nrej >= MAXREJECT) {
										printf("\nSorry: Too much rejections using the current parameter. Aborting simulation\n");
										fprintf(output,"\n\nSorry: Too much rejections using the current parameter. Aborting simulation\n");
										exit(2);
									}
									break; /*reject*/
								}
							}
							if((*inputp)->Sfix_alltheta == 1) {
								ui = neutpar[0].S;
								if(ui != (long int)(*inputp)->linked_segsites[aa]) {
									jcount += 1;
									if(jcount - nrej >= MAXREJECT) {
										printf("\nSorry: Too much rejections using the current parameter. Aborting simulation\n");
										fprintf(output,"\n\nSorry: Too much rejections using the current parameter. Aborting simulation\n");
										exit(2);
									}
									break; /*reject*/
								}
							}
						}
						if(aa<(*inputp)->linked) continue; /*reject*/
						else {
							/*accept*/ /*We will be in this loop until we have "howmany" success*/
							if((*inputp)->ifgamma == 1 || (*inputp)->range_thetant) 
								postp[(*inputp)->nloci][listnumbers[count-1]].thetap = thetae;
							if((*inputp)->ifgammar == 1 || (*inputp)->range_rnt) 
								postp[(*inputp)->nloci][listnumbers[count-1]].recp = rece;
							if((*inputp)->ifgamma == 1 || (*inputp)->range_thetant || (*inputp)->ifgammar == 1 || (*inputp)->range_rnt)
								postp[(*inputp)->nloci][listnumbers[count-1]].Ttotp  = lengtht;
							/*
							if((*inputp)->Sfix_alltheta == 1) 
								postp[(*inputp)->nloci][listnumbers[count-1]].thetap = thetae;
							if((*inputp)->rmfix == 1) 
								postp[(*inputp)->nloci][listnumbers[count-1]].recp = rece;
							postp[(*inputp)->nloci][listnumbers[count-1]].Ttotp  = lengtht;
							*/
							jcount++;
							
							/*printing a point every 2% of the total iterations*/
							#if SHOWPROGRESS == 1
							if((double)p+restp >= counterp) {
								restp += (double)p - counterp; 
								if((double)restp/(double)counterp > (double)1) {
									for(x=0;x<(int)floor(restp/counterp);x++) printf(".");
									restp -= (double)floor(restp/counterp) * counterp;
								}
								else printf(".");
								fflush(stdout);
								p = 1;
							}
							else p += 1;
							#endif

							break;
						}
					}
					else {/*like not linked fragments*/
						/**/calc_neutpar(0,segsites,inputp,neutpar+0,(double)recombinationv);/**/
						if((*inputp)->rmfix == 1) {
							/*
							if(recombinationv) Rmi = Min_rec(0,segsites,(*inputp)->nsam,0);
							else Rmi = 0;
							ZAi = (double)ZnA_(0,segsites,inputp);
							if(Rmi != (*inputp)->Rm || (ZAi < (*inputp)->ZA - (double)0.1 || ZAi > (*inputp)->ZA +(double)0.1)) {
							*/
							Rmi = neutpar[0].Rm;
							nhi = neutpar[0].nhapl;
							if(Rmi != (*inputp)->Rm || ((*inputp)->nhapl != 0 && nhi != (*inputp)->nhapl)) {
								jcount += 1;
								if(jcount - nrej >= MAXREJECT) {
									printf("\nSorry: Too much rejections using the current parameter. Aborting simulation\n");
									fprintf(output,"\n\nSorry: Too much rejections using the current parameter. Aborting simulation\n");
									exit(2);
								}
								continue; /*reject*/
							}
						}
						if((*inputp)->Sfix_alltheta == 1 && (*inputp)->mhits == 0) {
							u = (double)logPPoisson2((long int)(*inputp)->segsitesin,thetae*lengtht) - logPoissonkk;
							logv = (double)log((double)ran1());
							if(u < logv) {
								jcount++;
								if(jcount - nrej >= MAXREJECT) {
									printf("\nSorry: Too much rejections using the current parameter. Aborting simulation\n");
									fprintf(output,"\n\nSorry: Too much rejections using the current parameter. Aborting simulation\n");
									exit(2);
								}
								continue; /*reject*/
							}
						}
						if((*inputp)->Sfix_alltheta == 1 && (*inputp)->mhits == 1) {
							ui = neutpar[0].S;
							if(ui != (long int)(*inputp)->segsitesin) {
								jcount += 1;
								if(jcount - nrej >= MAXREJECT) {
									printf("\nSorry: Too much rejections using the current parameter. Aborting simulation\n");
									fprintf(output,"\n\nSorry: Too much rejections using the current parameter. Aborting simulation\n");
									exit(2);
								}
								continue; /*reject*/
							}
						}
						/*accept*/ /*We will be in this loop until we have "howmany" success*/
						if((*inputp)->ifgamma == 1 || (*inputp)->range_thetant) 
							postp[(*inputp)->nloci][listnumbers[count-1]].thetap = thetae;
						if((*inputp)->ifgammar == 1 || (*inputp)->range_rnt) 
							postp[(*inputp)->nloci][listnumbers[count-1]].recp = rece;
						if((*inputp)->ifgamma == 1 || (*inputp)->range_thetant || (*inputp)->ifgammar == 1 || (*inputp)->range_rnt)
							postp[(*inputp)->nloci][listnumbers[count-1]].Ttotp  = lengtht;
						/*
						if((*inputp)->Sfix_alltheta == 1) 
							postp[(*inputp)->nloci][listnumbers[count-1]].thetap = thetae;
						if((*inputp)->rmfix == 1) 
							postp[(*inputp)->nloci][listnumbers[count-1]].recp = rece;
						postp[(*inputp)->nloci][listnumbers[count-1]].Ttotp  = lengtht;
						*/
						jcount++;

						/*printing a point every 2% of the total iterations*/
						#if SHOWPROGRESS == 1
						if((double)p+restp >= counterp) {
							restp += (double)p - counterp; 
							if((double)restp/(double)counterp > (double)1) {
								for(x=0;x<(int)floor(restp/counterp);x++) printf(".");
								restp -= (double)floor(restp/counterp) * counterp;
							}
							else printf(".");
							fflush(stdout);
							p = 1;
						}
						else p += 1;
						#endif

						break;
					}
				}while(1);
			}    
			/*MHMCMC */
			/*We do the comparison u = Po(S,theta_canddt*Ln_canddt)/Po(S,theta_mc*Ln_mc).*/
			/*When values were rejected, we chose the anterior value of the chain.*/
			else { /*if((*inputp)->Sfix_alltheta == 2 || (*inputp)->rmfix == 2)*/
				do {
					if(mcount == 0) {
						if((*inputp)->ifgamma == 1) theta_mc = gammadist((*inputp)->alpha_gamma)/(*inputp)->p_gamma * (*inputp)->correct_gamma;
						else {
							if((*inputp)->range_thetant == 1) theta_mc = thetaemin + (thetaemax - thetaemin) * ran1();
							else {
								if((*inputp)->range_thetant == 2) theta_mc = (double)exp((double)log((double)thetaemin) + ((double)log((double)thetaemax) - (double)log((double)thetaemin)) * ran1());
								else theta_mc = (double)(*inputp)->theta;
							}
						}
						if((*inputp)->ifgammar == 1) rec_mc = gammadist((*inputp)->alpha_gammar)/(*inputp)->p_gammar * (*inputp)->correct_gammar;
						else {
							if((*inputp)->range_rnt == 1) rec_mc = recemin + (recemax - recemin) * ran1();
							else {
								if((*inputp)->range_rnt == 2) rec_mc = (double)exp((double)log((double)recemin) + ((double)log((double)recemax) - (double)log((double)recemin)) * ran1());
								else rec_mc = (double)(*inputp)->r;
							}
						}			
						theta_canddt = theta_mc;
						rec_canddt = rec_mc;
						/*assuming a value betweeen 0 and 1, small values are important...*/
						/*theta_mc = thetaemax / (1. + (thetaemax/thetaemin - 1.) * ran1());*/
					}
					else {
						if((*inputp)->ifgamma==1) {
							theta_canddt = gammadist((*inputp)->alpha_gamma)/(*inputp)->p_gamma * (*inputp)->correct_gamma;
						}
						else {
							if((*inputp)->Sfix_alltheta == 2) {
								if((*inputp)->range_thetant == 1) {
									if(thetaemax <= theta_mc + divt/(double)2) 
										theta_canddt = thetaemax - divt*ran1();
									else if(thetaemin > theta_mc - divt/(double)2)
										theta_canddt = thetaemin + divt*ran1();
									else theta_canddt = theta_mc + divt*(ran1() - (double)0.5); /*uniform prior (but extrems)*/
								}
								else {
									if((*inputp)->range_thetant == 2) {
										if((double)log(thetaemax) <= (double)log(theta_mc) + divt/(double)2) 
											theta_canddt = (double)exp((double)log(thetaemax) - divt*ran1());
										else if((double)log(thetaemin) > (double)log(theta_mc) - divt/(double)2)
											theta_canddt = (double)exp((double)log(thetaemin) + divt*ran1());
										else theta_canddt = (double)exp((double)log(theta_mc) + divt*(ran1() - (double)0.5)); /*uniform prior (but extrems)*/
									}
									else theta_canddt = (double)(*inputp)->theta;
								}
							}
							else {
								if((*inputp)->range_thetant == 1) theta_canddt = thetaemin + (thetaemax - thetaemin) * ran1();
								{
									if((*inputp)->range_thetant == 2) theta_canddt = (double)exp((double)log((double)thetaemin) + ((double)log((double)thetaemax) - (double)log((double)thetaemin)) * ran1());
									else theta_canddt = (double)(*inputp)->theta;
								}
							}
						}
						
						if((*inputp)->ifgammar==1) {
							rec_canddt = gammadist((*inputp)->alpha_gammar)/(*inputp)->p_gammar * (*inputp)->correct_gamma;
						}
						else {
							if((*inputp)->rmfix == 2) {
								if((*inputp)->range_rnt == 1) {
									if(recemax <= rec_mc + divr/(double)2) 
										rec_canddt = recemax - divr*ran1();
									else if(recemin > rec_mc - divr/(double)2)
										rec_canddt = recemin + divr*ran1();
									else rec_canddt = rec_mc + divt*(ran1() - (double)0.5); /*uniform prior (but extrems)*/
								}
								else {
									if((*inputp)->range_rnt == 2) {
										if((double)log(recemax) <= (double)log(rec_mc) + divt/(double)2) 
											rec_canddt = (double)exp((double)log(recemax) - divt*ran1());
										else if((double)log(recemin) > (double)log(rec_mc) - divt/(double)2)
											rec_canddt = (double)exp((double)log(recemin) + divt*ran1());
										else rec_canddt = (double)exp((double)log(rec_mc) + divt*(ran1() - (double)0.5)); /*uniform prior (but extrems)*/
									}
									rec_canddt = (double)(*inputp)->r;
								}
							}
							else {
								if((*inputp)->range_rnt == 1) rec_canddt = recemin + (recemax - recemin) * ran1();
								{
									if((*inputp)->range_rnt == 2) rec_canddt = (double)exp((double)log((double)recemin) + ((double)log((double)recemax) - (double)log((double)recemin)) * ran1());
									else rec_canddt = (double)(*inputp)->r;
								}
							}
						}
					}
					recombinationv = correction_rec(rec_canddt,(*inputp)->factor_chrn,(*inputp)->no_rec_males);
					if((*inputp)->segsitesin == -1 || (*inputp)->mhits == 1) thetaSv = theta_canddt;
					else thetaSv = (double)1.0;
					/*be careful with mhits*/
					inputS = (*inputp)->segsitesin;
					if((*inputp)->mhits == 1) inputS = -1;

					do_heter_gamma_sites(weightmut,(*inputp)->heter_theta_alphag,thetaSv+thetaSv/(double)(*inputp)->nsites,(*inputp)->nsites+1,(*inputp)->invariable_mut_sites);
					do_heter_gamma_sitesrec(weightrec,(*inputp)->heter_rm_alphag,(double)recombinationv/*+recombinationv/(double)(*inputp)->nsites*/,(*inputp)->nsites);
					segsites = gensam((*inputp)->npop,(*inputp)->nsam, (*inputp)->config, (*inputp)->nsites,
					theta_canddt, /*(*inputp)->segsitesin*/inputS, recombinationv, (*inputp)->f, (*inputp)->track_len,
					(*inputp)->migrate,(*inputp)->mhits,count,(*inputp)->factor_pop, 
					&lengtht,(*inputp)->ifselection, 
					(*inputp)->pop_sel,(*inputp)->sinit,(*inputp)->pop_size,(*inputp)->sel_nt,(*inputp)->T_out,
					&(*inputp)->Sout,(*inputp)->nintn, (*inputp)->nrec, (*inputp)->npast, (*inputp)->tpast,
					(*inputp)->split_pop,(*inputp)->time_split,(*inputp)->time_scoal,(*inputp)->factor_anc,(*inputp)->freq,
					(*inputp)->tlimit,(*inputp)->iflogistic,(*inputp)->ts,(*inputp)->factor_chrn,(*inputp)->ratio_sv,
					weightmut,weightrec);
					if((*inputp)->mhits) mod_mhits(segsites,inputp); /******** mhits ******************/
					
					/*to avoid correlations only chose values "sep" positions separated */                        
					if(k == sep) {
						if(kcount > 1) {
							/*move the new value in matrix_test*/
							for(j=0;j<neutv;j++)
								matrix_test[(*inputp)->nloci*neutv+j][listnumbers[count-kcount-1]] 
									= matrix_test[(*inputp)->nloci*neutv+j][listnumbers[count-2]];							
							count -= kcount-1;							
						}
						if(kcount == 0 && count < ((*inputp)->howmany + sep)) {
							/*move the anterior value in matrix_test*/
							for(j=0;j<neutv;j++)
								matrix_test[(*inputp)->nloci*neutv+j][listnumbers[count-1]]
									= matrix_test[(*inputp)->nloci*neutv+j][listnumbers[count-2]];
							count++;
						}
						kcount = 0;
						k = 0;

						if((*inputp)->ifgamma == 1 || (*inputp)->range_thetant) postp[(*inputp)->nloci][listnumbers[count-1]].thetap = theta_mc;
						if((*inputp)->ifgammar == 1 || (*inputp)->range_rnt) postp[(*inputp)->nloci][listnumbers[count-1]].recp = rec_mc;
						postp[(*inputp)->nloci][listnumbers[count-1]].Ttotp  = lengtht_mc;

						/*printing a point every 2% of the total iterations*/
						#if SHOWPROGRESS == 1
						if((double)p+restp >= counterp) {
							restp += (double)p - counterp; 
							if((double)restp/(double)counterp > (double)1) {
								for(x=0;x<(int)floor(restp/counterp);x++) printf(".");
								restp -= (double)floor(restp/counterp) * counterp;
							}
							else printf(".");
							fflush(stdout);
							p = 1;
						}
						else p += 1;
						#endif

						if(count == ((*inputp)->howmany + sep)) {
							break; /*end*/
						}
					}
					/*Acceptation/Rejection*/
					if(mcount == 0) {
						if((*inputp)->linked > 1) {/*linked fragments with fixed subset values*/
							s0=s1=0;
							nwindow=0;
							for(aa=0;aa<(*inputp)->linked;aa++) {
								kk=(*inputp)->loci_linked[aa][0];
								ll=(*inputp)->loci_linked[aa][1]+1;
								while(s0 < segsites && posit[s0] < kk) s0++;
								while(s1 < segsites && posit[s1] < ll) s1++;
								/*calc_estadistics*/
								calc_neutpar_window(inputp,neutpar,s0,s1,(double)recombinationv);
								if((*inputp)->rmfix == 1) {
									Rmi = neutpar[0].Rm;
									nhi = neutpar[0].nhapl;
									if(Rmi != (*inputp)->linked_rm[aa] || ((*inputp)->linked_nhapl[aa] != 0 && nhi != (*inputp)->linked_nhapl[aa])) {
										break; /*reject*/
									}
								}
								if((*inputp)->Sfix_alltheta == 1) {
									ui = neutpar[0].S;
									if(ui != (long int)(*inputp)->linked_segsites[aa]) {
										break; /*reject*/
									}
								}
							}
							if(aa<(*inputp)->linked) continue; /*reject*/
							else {
								/*accept*/
								if((*inputp)->Sfix_alltheta == 2)
									logprob_mc = (double)logPPoisson2((long int)(*inputp)->segsitesin,theta_mc*lengtht);
								/*in case showing the trees sizes*/
								lengtht_mc = lengtht;

								kcount++;
								jcount++;
								mcount++;
								k++;

								if((*inputp)->Sfix_alltheta == 2) postp[(*inputp)->nloci][listnumbers[count-1]].thetap = theta_mc;
								if((*inputp)->rmfix == 2) postp[(*inputp)->nloci][listnumbers[count-1]].recp = rec_mc;
								postp[(*inputp)->nloci][listnumbers[count-1]].Ttotp  = lengtht_mc;

								#if SHOWPROGRESS == 1
								if((double)p+restp >= counterp) {
									restp += (double)p - counterp; 
									if((double)restp/(double)counterp > (double)1) {
										for(x=0;x<(int)floor(restp/counterp);x++) printf(".");
										restp -= (double)floor(restp/counterp) * counterp;
									}
									else printf(".");
									fflush(stdout);
									p = 1;
								}
								else p += 1;
								#endif

								break;
							}
						}
						else {
							/**/calc_neutpar(0,segsites,inputp,neutpar+0,(double)recombinationv);/**/
							if((*inputp)->rmfix == 2) {
								/*
								if(recombinationv) Rmi = Min_rec(0,segsites,(*inputp)->nsam,0);
								else Rmi = 0;
								ZAi = (double)ZnA_(0,segsites,inputp);
								if(Rmi != (*inputp)->Rm || (ZAi < (*inputp)->ZA -(double)0.1 || ZAi > (*inputp)->ZA +(double)0.1)) {
								*/
								Rmi = neutpar[0].Rm;
								nhi = neutpar[0].nhapl;
								if(Rmi != (*inputp)->Rm || ((*inputp)->nhapl != 0 && nhi != (*inputp)->nhapl)) {
									continue;
								}
							}
							if((*inputp)->Sfix_alltheta == 2 && (*inputp)->mhits == 0) {
								u = (double)logPPoisson2((long int)(*inputp)->segsitesin,theta_mc*lengtht) - logPoissonkk;
								if(u < (double)log((double)ran1())) {
									continue;
								}
							}						
							if((*inputp)->Sfix_alltheta == 2)
								logprob_mc = (double)logPPoisson2((long int)(*inputp)->segsitesin,theta_mc*lengtht);

							if((*inputp)->Sfix_alltheta == 2 && (*inputp)->mhits == 1) {
								u = neutpar[0].S;
								if(u != (long int)(*inputp)->segsitesin) {
									continue; /*reject*/
								}
							}
							/*in case showing the trees sizes*/
							lengtht_mc = lengtht;

							kcount++;
							jcount++;
							mcount++;
							k++;

							if((*inputp)->Sfix_alltheta == 2) postp[(*inputp)->nloci][listnumbers[count-1]].thetap = theta_mc;
							if((*inputp)->rmfix == 2) postp[(*inputp)->nloci][listnumbers[count-1]].recp = rec_mc;
							postp[(*inputp)->nloci][listnumbers[count-1]].Ttotp  = lengtht_mc;

							#if SHOWPROGRESS == 1
							if((double)p+restp >= counterp) {
								restp += (double)p - counterp; 
								if((double)restp/(double)counterp > (double)1) {
									for(x=0;x<(int)floor(restp/counterp);x++) printf(".");
									restp -= (double)floor(restp/counterp) * counterp;
								}
								else printf(".");
								fflush(stdout);
								p = 1;
							}
							else p += 1;
							#endif

							break;
						}
					}
					else {
						if((*inputp)->linked > 1) {/*linked fragments with subset fixed values*/
							s0=s1=0;
							nwindow=0;
							for(aa=0;aa<(*inputp)->linked;aa++) {
								kk=(*inputp)->loci_linked[aa][0];
								ll=(*inputp)->loci_linked[aa][1]+1;
								while(s0 < segsites && posit[s0] < kk) s0++;
								while(s1 < segsites && posit[s1] < ll) s1++;
								/*calc_estadistics*/
								calc_neutpar_window(inputp,neutpar,s0,s1,(double)recombinationv);
								if((*inputp)->rmfix == 1) {
									Rmi = neutpar[0].Rm;
									nhi = neutpar[0].nhapl;
									if(Rmi != (*inputp)->linked_rm[aa] || ((*inputp)->linked_nhapl[aa] != 0 && nhi != (*inputp)->linked_nhapl[aa])) {
										mcount++;
										k++;
										break; /*reject*/
									}
								}
								if((*inputp)->Sfix_alltheta == 1) {
									logprob_canddt = (double)logPPoisson2((long int)(*inputp)->segsitesin,theta_canddt*lengtht);
									ui = neutpar[0].S;
									if(ui != (long int)(*inputp)->linked_segsites[aa]) {
										mcount++;
										k++;
										break; /*reject*/
									}
								}
							}
							if(aa<(*inputp)->linked) continue; /*reject*/
							else {
								/*accept*/
								if((*inputp)->rmfix == 2)
									rec_mc = rec_canddt;
								if((*inputp)->Sfix_alltheta == 2) {
									theta_mc = theta_canddt;
									logprob_mc = logprob_canddt;
								}
								/*in case showing the trees sizes*/
								lengtht_mc = lengtht;

								kcount++;
								jcount++;
								mcount++;
								k++;

								break; /*accept*/
							}
						}
						else {
							/**/calc_neutpar(0,segsites,inputp,neutpar+0,(double)recombinationv);/**/
							if((*inputp)->rmfix == 2) {
								/*
								if(recombinationv) Rmi = Min_rec(0,segsites,(*inputp)->nsam,0);
								else Rmi = 0;
								ZAi = (double)ZnA_(0,segsites,inputp);
								if(Rmi != (*inputp)->Rm || (ZAi < (*inputp)->ZA -(double)0.1 || ZAi > (*inputp)->ZA +(double)0.1)) {
								*/
								Rmi = neutpar[0].Rm;
								nhi = neutpar[0].nhapl;
								if(Rmi != (*inputp)->Rm || ((*inputp)->nhapl != 0 && nhi != (*inputp)->nhapl)) {
									mcount++;
									k++;
									continue;
								}
							}
							if((*inputp)->Sfix_alltheta == 2 && (*inputp)->mhits == 0) {
								logprob_canddt = (double)logPPoisson2((long int)(*inputp)->segsitesin,theta_canddt*lengtht);
								u = logprob_canddt - logprob_mc;
								u += (double)log((double)(weight_thcanddt/weight_thmc)); /*logpr_prob - logprior_mc*/
								if(u < (double)log((double)ran1())) {  /*reject*/
									mcount++;
									k++;
									continue; 
								}
							}
							if((*inputp)->Sfix_alltheta == 2 && (*inputp)->mhits == 1) {
								u = neutpar[0].S;
								if(u != (long int)(*inputp)->segsitesin) {
									mcount++;
									k++;
									continue; /*reject*/
								}
							}
							/*accept*/
							if((*inputp)->rmfix == 2)
								rec_mc = rec_canddt;
							if((*inputp)->Sfix_alltheta == 2) {
								theta_mc = theta_canddt;
								logprob_mc = logprob_canddt;
							}
							/*in case showing the trees sizes*/
							lengtht_mc = lengtht;

							kcount++;
							jcount++;
							mcount++;
							k++;

							break; /*accept*/
						}
					}
				} while(1);
			} 
		}
		
		if((*inputp)->pr_matrix) /******** print the whole matrix **********/
			print_matrix(segsites,inputp,output,(*inputp)->pr_matrix,count);
		else { /******************************* STATISTICS *******************************/
			if((*inputp)->linked > 0) {
				/*define the limits for each region (min and max values)*/
				if((*inputp)->linked == 1) {/*sliding windows*/
					s0=s1=0;
					nwindow=0;
					for(kk=0,ll=(*inputp)->window;kk<(*inputp)->nsites;kk += (long int)(*inputp)->despl) {
						while(posit[s0]<kk && s0 < segsites) s0++;
						while(posit[s1]<ll && s1 < segsites) s1++;
						/*calc_estadistics*/
						calc_neutpar_window(inputp,neutpar,s0,s1,(double)recombinationv);
						/*calc_neut_tests*/
						/*we need to define the number of windows and the current window number*/
						matrix_test[nwindow*NEUTVALUES2+0][listnumbers[count-1]] 
							= (double)tajima_d(neutpar[0].k,neutpar[0].S,coef[0]);
						matrix_test[nwindow*NEUTVALUES2+1][listnumbers[count-1]] 
							= (double)Fs((*inputp)->nsam,neutpar[0].k,neutpar[0].nhapl);
						matrix_test[nwindow*NEUTVALUES2+2][listnumbers[count-1]] 
							= (double)fl_d2((*inputp)->config[0],
											neutpar[0].freq[1]+neutpar[0].freq[(*inputp)->config[0]-1],
											neutpar[0].S,coef[0]);
						matrix_test[nwindow*NEUTVALUES2+3][listnumbers[count-1]] 
							= (double)fl_f2((*inputp)->config[0],
											neutpar[0].freq[1]+neutpar[0].freq[(*inputp)->config[0]-1],
											neutpar[0].S, neutpar[0].k,coef[0]);
						matrix_test[nwindow*NEUTVALUES2+4][listnumbers[count-1]] 
							= (double)fl_d((*inputp)->nsam,neutpar[0].freq[1],neutpar[0].S,coef[0]);
						matrix_test[nwindow*NEUTVALUES2+5][listnumbers[count-1]] 
							= (double)fl_f((*inputp)->nsam,neutpar[0].freq[1],neutpar[0].S, neutpar[0].k,coef[0]);
						matrix_test[nwindow*NEUTVALUES2+6][listnumbers[count-1]] 
							= (double)fay_wu((*inputp)->nsam,neutpar[0].freq,neutpar[0].k);
						if(neutpar[0].S > 1) {
							matrix_test[nwindow*NEUTVALUES2+7][listnumbers[count-1]] 
								= (double)neutpar[0].B1/((double)neutpar[0].S - (double)1.0);
							matrix_test[nwindow*NEUTVALUES2+8][listnumbers[count-1]] 
								= (double)neutpar[0].Q1/((double)neutpar[0].S);
							matrix_test[nwindow*NEUTVALUES2+9][listnumbers[count-1]] 
								= /*-10000;*//**/(double)ZnA_window(inputp,s0,s1);/**//*CHECKING*/
						}
						else {
							matrix_test[nwindow*NEUTVALUES2+7][listnumbers[count-1]] = -10000;
							matrix_test[nwindow*NEUTVALUES2+8][listnumbers[count-1]] = -10000;
							matrix_test[nwindow*NEUTVALUES2+9][listnumbers[count-1]] = -10000;
						}
						matrix_test[nwindow*NEUTVALUES2+10][listnumbers[count-1]] 
								= (double)Fst(neutpar[0].piw,neutpar[0].pib,(*inputp)->npop);
						
						matrix_test[nwindow*NEUTVALUES2+11][listnumbers[count-1]] 
							= (double)neutpar[0].nhapl/(double)(*inputp)->config[0];
						matrix_test[nwindow*NEUTVALUES2+12][listnumbers[count-1]] 
							= (double)testHap((*inputp)->config[0],neutpar[0].fhapl);
						if(neutpar[0].S > 0)
							matrix_test[nwindow*NEUTVALUES2+13][listnumbers[count-1]] 
								= (double)R2(neutpar[0].unic,neutpar[0].k,(*inputp)->config[0],neutpar[0].S);
						else 
							matrix_test[nwindow*NEUTVALUES2+13][listnumbers[count-1]] = -10000;
						matrix_test[nwindow*NEUTVALUES2+14][listnumbers[count-1]] 
							= (double)neutpar[0].S;/*Gxi((*inputp)->config[0],neutpar[0].freq,(double)neutpar[0].S/coef[0][0]);*//*MODIFIED*/
						matrix_test[nwindow*NEUTVALUES2+15][listnumbers[count-1]] 
							= (double)neutpar[0].piw;
						matrix_test[nwindow*NEUTVALUES2+16][listnumbers[count-1]] 
							= (double)neutpar[0].pib;
						/*S, pi, thetaFW*/
						matrix_test[nwindow*NEUTVALUES2+17][listnumbers[count-1]] 
							= (double)neutpar[0].S/(double)coef[0][0];
						matrix_test[nwindow*NEUTVALUES2+18][listnumbers[count-1]] 
							= (double)neutpar[0].k;
						matrix_test[nwindow*NEUTVALUES2+19][listnumbers[count-1]] 
							= (double)neutpar[0].k - 
							  (double)(matrix_test[nwindow*NEUTVALUES2+6][listnumbers[count-1]]);
						matrix_test[nwindow*NEUTVALUES2+20][listnumbers[count-1]] 
							= (double)tajima_dvsdmin(neutpar[0].k,neutpar[0].S,coef[0],(*inputp)->nsam);
						matrix_test[nwindow*NEUTVALUES2+21][listnumbers[count-1]] 
							/*= (double)fay_wu_normalized((int)(*inputp)->nsam,neutpar[0].freq,neutpar[0].k,neutpar[0].S);*/
							= (double)fay_wu_normalized2((int)(*inputp)->nsam,(double)neutpar[0].thetaL,(double)neutpar[0].S/(double)coef[0][0],(double)neutpar[0].S,coef[0],neutpar[0].k);
						matrix_test[nwindow*NEUTVALUES2+22][listnumbers[count-1]] 
							= (double)neutpar[0].maxhapl/(double)(*inputp)->nsam;
						matrix_test[nwindow*NEUTVALUES2+23][listnumbers[count-1]] 
							= (double)neutpar[0].maxhapl1/(double)(*inputp)->nsam;
						matrix_test[nwindow*NEUTVALUES2+24][listnumbers[count-1]] 
							= (double)neutpar[0].Rm;
						matrix_test[nwindow*NEUTVALUES2+25][listnumbers[count-1]] 
							= (double)neutpar[0].freq[1];
						matrix_test[nwindow*NEUTVALUES2+26][listnumbers[count-1]] 
							= (double)neutpar[0].thetaL;
						matrix_test[nwindow*NEUTVALUES2+27][listnumbers[count-1]] 
							= (double)E_zeng((int)(*inputp)->nsam,(double)neutpar[0].thetaL,(double)neutpar[0].S/(double)coef[0][0],(double)neutpar[0].S,coef[0]);
						matrix_test[nwindow*NEUTVALUES2+28][listnumbers[count-1]] 
							= (double)EWtest((int)(*inputp)->nsam,neutpar[0].fhapl);
						matrix_test[nwindow*NEUTVALUES2+29][listnumbers[count-1]] 
							= (double)Fst(neutpar[0].withinw,neutpar[0].pib,(*inputp)->npop);
						if(neutpar[0].S > 0) {
						matrix_test[nwindow*NEUTVALUES2+30][listnumbers[count-1]] 
							= (double)pwh((int)(*inputp)->nsam,neutpar[0].pid)/(double)neutpar[0].k;
						}
						else {
							matrix_test[nwindow*NEUTVALUES2+30][listnumbers[count-1]] = (double)-10000;
						}
						/*end calc neutrality test*/
						nwindow += 1;
						if(ll == (*inputp)->nsites) break;
						else ll += (*inputp)->despl;
						if(ll > (*inputp)->nsites) ll = (*inputp)->nsites;
						s1 = s0;
					}
				}
				else {/*linked regions*/
					s0=s1=0;
					nwindow=0;
					for(aa=0;aa<(*inputp)->linked;aa++) {
						kk=(*inputp)->loci_linked[aa][0];
						ll=(*inputp)->loci_linked[aa][1]+1;
						
						while(s0 < segsites && posit[s0] < kk) {
							s0++;
						}
						
						while(s1 < segsites && posit[s1] < ll) {
							s1++;
						}
						
						/*calc_estadistics*/
						calc_neutpar_window(inputp,neutpar,s0,s1,(double)recombinationv);
						/*calc_neut_tests*/
						matrix_test[nwindow*NEUTVALUES2+0][listnumbers[count-1]] 
							= (double)tajima_d(neutpar[0].k,neutpar[0].S,coef[0]);
						matrix_test[nwindow*NEUTVALUES2+1][listnumbers[count-1]] 
							= (double)Fs((*inputp)->nsam,neutpar[0].k,neutpar[0].nhapl);
						matrix_test[nwindow*NEUTVALUES2+2][listnumbers[count-1]] 
							= (double)fl_d2((*inputp)->config[0],
											neutpar[0].freq[1]+neutpar[0].freq[(*inputp)->config[0]-1],
											neutpar[0].S,coef[0]);
						matrix_test[nwindow*NEUTVALUES2+3][listnumbers[count-1]] 
							= (double)fl_f2((*inputp)->config[0],
											neutpar[0].freq[1]+neutpar[0].freq[(*inputp)->config[0]-1],
											neutpar[0].S, neutpar[0].k,coef[0]);
						matrix_test[nwindow*NEUTVALUES2+4][listnumbers[count-1]] 
							= (double)fl_d((*inputp)->nsam,neutpar[0].freq[1],neutpar[0].S,coef[0]);
						matrix_test[nwindow*NEUTVALUES2+5][listnumbers[count-1]] 
							= (double)fl_f((*inputp)->nsam,neutpar[0].freq[1],neutpar[0].S, neutpar[0].k,coef[0]);
						matrix_test[nwindow*NEUTVALUES2+6][listnumbers[count-1]] 
							= (double)fay_wu((*inputp)->nsam,neutpar[0].freq,neutpar[0].k);
						if(neutpar[0].S > 1) {
							matrix_test[nwindow*NEUTVALUES2+7][listnumbers[count-1]] 
								= (double)neutpar[0].B1/((double)neutpar[0].S - (double)1.0);
							matrix_test[nwindow*NEUTVALUES2+8][listnumbers[count-1]] 
								= (double)neutpar[0].Q1/((double)neutpar[0].S);
							matrix_test[nwindow*NEUTVALUES2+9][listnumbers[count-1]] 
								= /*-10000;*//**/(double)ZnA_window(inputp,s0,s1);/**//*CHECKING*/
						}
						else {
							matrix_test[nwindow*NEUTVALUES2+7][listnumbers[count-1]] = -10000;
							matrix_test[nwindow*NEUTVALUES2+8][listnumbers[count-1]] = -10000;
							matrix_test[nwindow*NEUTVALUES2+9][listnumbers[count-1]] = -10000;
						}
						matrix_test[nwindow*NEUTVALUES2+10][listnumbers[count-1]] 
								= (double)Fst(neutpar[0].piw,neutpar[0].pib,(*inputp)->npop);
						
						matrix_test[nwindow*NEUTVALUES2+11][listnumbers[count-1]] 
							= (double)neutpar[0].nhapl/(double)(*inputp)->config[0];
						matrix_test[nwindow*NEUTVALUES2+12][listnumbers[count-1]] 
							= (double)testHap((*inputp)->config[0],neutpar[0].fhapl);
						if(neutpar[0].S > 0)
							matrix_test[nwindow*NEUTVALUES2+13][listnumbers[count-1]] 
								= (double)R2(neutpar[0].unic,neutpar[0].k,(*inputp)->config[0],neutpar[0].S);
						else 
							matrix_test[nwindow*NEUTVALUES2+13][listnumbers[count-1]] = -10000;
						matrix_test[nwindow*NEUTVALUES2+14][listnumbers[count-1]] 
							= (double)neutpar[0].S;/*Gxi((*inputp)->config[0],neutpar[0].freq,(double)neutpar[0].S/coef[0][0]);*//*MODIFIED*/
						matrix_test[nwindow*NEUTVALUES2+15][listnumbers[count-1]] 
							= (double)neutpar[0].piw;
						matrix_test[nwindow*NEUTVALUES2+16][listnumbers[count-1]] 
							= (double)neutpar[0].pib;
						/*S, pi, thetaFW*/
						matrix_test[nwindow*NEUTVALUES2+17][listnumbers[count-1]] 
							= (double)neutpar[0].S/(double)coef[0][0];
						matrix_test[nwindow*NEUTVALUES2+18][listnumbers[count-1]] 
							= (double)neutpar[0].k;
						matrix_test[nwindow*NEUTVALUES2+19][listnumbers[count-1]] 
							= (double)neutpar[0].k - 
							  (double)(matrix_test[nwindow*NEUTVALUES2+6][listnumbers[count-1]]);
						matrix_test[nwindow*NEUTVALUES2+20][listnumbers[count-1]] 
							= (double)tajima_dvsdmin(neutpar[0].k,neutpar[0].S,coef[0],(*inputp)->nsam);
						matrix_test[nwindow*NEUTVALUES2+21][listnumbers[count-1]] 
							/*= (double)fay_wu_normalized((*inputp)->nsam,neutpar[0].freq,neutpar[0].k,neutpar[0].S);*/
							= (double)fay_wu_normalized2((int)(*inputp)->nsam,(double)neutpar[0].thetaL,(double)neutpar[0].S/(double)coef[0][0],(double)neutpar[0].S,coef[0],neutpar[0].k);
						matrix_test[nwindow*NEUTVALUES2+22][listnumbers[count-1]] 
							= (double)neutpar[0].maxhapl/(double)(*inputp)->nsam;
						matrix_test[nwindow*NEUTVALUES2+23][listnumbers[count-1]] 
							= (double)neutpar[0].maxhapl1/(double)(*inputp)->nsam;
						matrix_test[nwindow*NEUTVALUES2+24][listnumbers[count-1]] 
							= (double)neutpar[0].Rm;
						matrix_test[nwindow*NEUTVALUES2+25][listnumbers[count-1]] 
							= (double)neutpar[0].freq[1];
						matrix_test[nwindow*NEUTVALUES2+26][listnumbers[count-1]] 
							= (double)neutpar[0].thetaL;
						matrix_test[nwindow*NEUTVALUES2+27][listnumbers[count-1]] 
							= (double)E_zeng((int)(*inputp)->nsam,(double)neutpar[0].thetaL,(double)neutpar[0].S/(double)coef[0][0],(double)neutpar[0].S,coef[0]);
						matrix_test[nwindow*NEUTVALUES2+28][listnumbers[count-1]] 
							= (double)EWtest((int)(*inputp)->nsam,neutpar[0].fhapl);
						matrix_test[nwindow*NEUTVALUES2+29][listnumbers[count-1]] 
							= (double)Fst(neutpar[0].withinw,neutpar[0].pib,(*inputp)->npop);
						if(neutpar[0].S > 0) {
						matrix_test[nwindow*NEUTVALUES2+30][listnumbers[count-1]] 
							= (double)pwh((int)(*inputp)->nsam,neutpar[0].pid)/(double)neutpar[0].k;
						}
						else {
							matrix_test[nwindow*NEUTVALUES2+30][listnumbers[count-1]] = (double)-10000;
						}
					   /*end calc neutrality test*/
						nwindow += 1;
						s0 = s1;
					}
				}
			}
			else {
				/**/if(!((*inputp)->rmfix)) /**/
					calc_neutpar(0,segsites,inputp,neutpar+0,(double)recombinationv);
					
				matrix_test[(*inputp)->nloci*NEUTVALUES2+0][listnumbers[count-1]] 
					= (double)tajima_d(neutpar[0].k,neutpar[0].S,coef[0]);
				matrix_test[(*inputp)->nloci*NEUTVALUES2+1][listnumbers[count-1]] 
					= (double)Fs((*inputp)->nsam,neutpar[0].k,neutpar[0].nhapl);
				matrix_test[(*inputp)->nloci*NEUTVALUES2+2][listnumbers[count-1]] 
					= (double)fl_d2((*inputp)->config[0],
									neutpar[0].freq[1]+neutpar[0].freq[(*inputp)->config[0]-1],
									neutpar[0].S,coef[0]);
				matrix_test[(*inputp)->nloci*NEUTVALUES2+3][listnumbers[count-1]] 
					= (double)fl_f2((*inputp)->config[0],
									neutpar[0].freq[1]+neutpar[0].freq[(*inputp)->config[0]-1],
									neutpar[0].S, neutpar[0].k,coef[0]);
				matrix_test[(*inputp)->nloci*NEUTVALUES2+4][listnumbers[count-1]] 
					= (double)fl_d((*inputp)->nsam,neutpar[0].freq[1],neutpar[0].S,coef[0]);
				matrix_test[(*inputp)->nloci*NEUTVALUES2+5][listnumbers[count-1]] 
					= (double)fl_f((*inputp)->nsam,neutpar[0].freq[1],neutpar[0].S, neutpar[0].k,coef[0]);
				matrix_test[(*inputp)->nloci*NEUTVALUES2+6][listnumbers[count-1]] 
					= (double)fay_wu((*inputp)->nsam,neutpar[0].freq,neutpar[0].k);
				if(neutpar[0].S > 1) {
					matrix_test[(*inputp)->nloci*NEUTVALUES2+7][listnumbers[count-1]] 
						= (double)neutpar[0].B1/((double)neutpar[0].S - (double)1.0);
					matrix_test[(*inputp)->nloci*NEUTVALUES2+8][listnumbers[count-1]] 
						= (double)neutpar[0].Q1/((double)neutpar[0].S);
					matrix_test[(*inputp)->nloci*NEUTVALUES2+9][listnumbers[count-1]] 
						= /*-10000;*//**/(double)ZnA_(0,segsites,inputp);/**//*INACTIVE*/
				}
				else {
					matrix_test[(*inputp)->nloci*NEUTVALUES2+7][listnumbers[count-1]] = -10000;
					matrix_test[(*inputp)->nloci*NEUTVALUES2+8][listnumbers[count-1]] = -10000;
					matrix_test[(*inputp)->nloci*NEUTVALUES2+9][listnumbers[count-1]] = -10000;
				}
				matrix_test[(*inputp)->nloci*NEUTVALUES2+10][listnumbers[count-1]] 
						= (double)Fst(neutpar[0].piw,neutpar[0].pib,(*inputp)->npop);
				
				matrix_test[(*inputp)->nloci*NEUTVALUES2+11][listnumbers[count-1]] 
					= (double)neutpar[0].nhapl/(double)(*inputp)->config[0];
				matrix_test[(*inputp)->nloci*NEUTVALUES2+12][listnumbers[count-1]] 
					= (double)testHap((*inputp)->config[0],neutpar[0].fhapl);
				if(neutpar[0].S > 0)
					matrix_test[(*inputp)->nloci*NEUTVALUES2+13][listnumbers[count-1]] 
						= (double)R2(neutpar[0].unic,neutpar[0].k,(*inputp)->config[0],neutpar[0].S);
				else 
					matrix_test[(*inputp)->nloci*NEUTVALUES2+13][listnumbers[count-1]] = -10000;
				matrix_test[(*inputp)->nloci*NEUTVALUES2+14][listnumbers[count-1]] 
					= (double)neutpar[0].S;/*Gxi((*inputp)->config[0],neutpar[0].freq,(double)neutpar[0].S/coef[0][0]);*//*MODIFIED*/
				matrix_test[(*inputp)->nloci*NEUTVALUES2+15][listnumbers[count-1]] 
					= (double)neutpar[0].piw;
				matrix_test[(*inputp)->nloci*NEUTVALUES2+16][listnumbers[count-1]] 
					= (double)neutpar[0].pib;
				/*S, pi, thetaFW*/
				matrix_test[(*inputp)->nloci*NEUTVALUES2+17][listnumbers[count-1]] 
					= (double)neutpar[0].S/(double)coef[0][0];
				matrix_test[(*inputp)->nloci*NEUTVALUES2+18][listnumbers[count-1]] 
					= (double)neutpar[0].k;
				matrix_test[(*inputp)->nloci*NEUTVALUES2+19][listnumbers[count-1]] 
					= (double)neutpar[0].k - (double)(matrix_test[(*inputp)->nloci*NEUTVALUES2+6][listnumbers[count-1]]);
				matrix_test[(*inputp)->nloci*NEUTVALUES2+20][listnumbers[count-1]] 
					= (double)tajima_dvsdmin(neutpar[0].k,neutpar[0].S,coef[0],(*inputp)->nsam);
				matrix_test[(*inputp)->nloci*NEUTVALUES2+21][listnumbers[count-1]] 
					/*= (double)fay_wu_normalized((*inputp)->nsam,neutpar[0].freq,neutpar[0].k,neutpar[0].S);*/
					= (double)fay_wu_normalized2((int)(*inputp)->nsam,(double)neutpar[0].thetaL,(double)neutpar[0].S/(double)coef[0][0],(double)neutpar[0].S,coef[0],neutpar[0].k);
				matrix_test[(*inputp)->nloci*NEUTVALUES2+22][listnumbers[count-1]] 
					= (double)neutpar[0].maxhapl/(double)(*inputp)->nsam;
				matrix_test[(*inputp)->nloci*NEUTVALUES2+23][listnumbers[count-1]] 
					= (double)neutpar[0].maxhapl1/(double)(*inputp)->nsam;
				matrix_test[(*inputp)->nloci*NEUTVALUES2+24][listnumbers[count-1]] 
					= (double)neutpar[0].Rm;
				matrix_test[(*inputp)->nloci*NEUTVALUES2+25][listnumbers[count-1]] 
					= (double)neutpar[0].freq[1];
				matrix_test[(*inputp)->nloci*NEUTVALUES2+26][listnumbers[count-1]] 
					= (double)neutpar[0].thetaL;
				matrix_test[(*inputp)->nloci*NEUTVALUES2+27][listnumbers[count-1]] 
					= (double)E_zeng((int)(*inputp)->nsam,(double)neutpar[0].thetaL,(double)neutpar[0].S/(double)coef[0][0],(double)neutpar[0].S,coef[0]);
				matrix_test[(*inputp)->nloci*NEUTVALUES2+28][listnumbers[count-1]] 
					= (double)EWtest((int)(*inputp)->nsam,neutpar[0].fhapl);
				matrix_test[(*inputp)->nloci*NEUTVALUES2+29][listnumbers[count-1]] 
					= (double)Fst(neutpar[0].withinw,neutpar[0].pib,(*inputp)->npop);
				matrix_test[(*inputp)->nloci*NEUTVALUES2+30][listnumbers[count-1]] 
					= (double)pwh((int)(*inputp)->nsam,neutpar[0].pid)/(double)neutpar[0].k;
				if(neutpar[0].S > 0) {
					matrix_test[(*inputp)->nloci*NEUTVALUES2+30][listnumbers[count-1]] 
						= (double)pwh((int)(*inputp)->nsam,neutpar[0].pid)/(double)neutpar[0].k;
				}
				else {
					matrix_test[(*inputp)->nloci*NEUTVALUES2+30][listnumbers[count-1]] = (double)-10000;
				}
				#if DEBUGSOUT
				fprintf(file_debug,"%d\n",(*inputp)->Sout);
				#endif
				#if PRINTTHETAS
				fprintf(file_debug_sel,"%f\t%f\t%f\n",neutpar[0].k,(double)neutpar[0].S/coef[0][0], 
					(double)neutpar[0].k - (double)(matrix_test[(*inputp)->nloci*NEUTVALUES2+6][listnumbers[count-1]]));
				#endif
			}
		}
	}		
	if((*inputp)->Sfix_alltheta == 1) 
		postp[(*inputp)->nloci][(*inputp)->howmany].thetap = postp[(*inputp)->nloci][(*inputp)->howmany].Ttotp = (double)(*inputp)->howmany/(double)jcount;
	if((*inputp)->Sfix_alltheta == 2) 
		postp[(*inputp)->nloci][(*inputp)->howmany].thetap = postp[(*inputp)->nloci][(*inputp)->howmany].Ttotp = (double)jcount/(double)mcount;
	if((*inputp)->rmfix == 1) 
		postp[(*inputp)->nloci][(*inputp)->howmany].recp = postp[(*inputp)->nloci][(*inputp)->howmany].Ttotp = (double)(*inputp)->howmany/(double)jcount;
	if((*inputp)->rmfix == 2) 
		postp[(*inputp)->nloci][(*inputp)->howmany].recp = postp[(*inputp)->nloci][(*inputp)->howmany].Ttotp = (double)jcount/(double)mcount;

    /* alliberar les matrius i vectors */
    free(posit);
    for(i=0;i<(*inputp)->nsam;i++)
        free(list[i]);
    free(list);
    /*if((*inputp)->neutral_tests) {*/
        free(neutpar[0].freq);
        free(neutpar[0].fhapl);
        free(neutpar[0].unic);
        free(neutpar[0].pid);
    /*}*/
    free(listnumbers);

	free(weightmut);
	free(weightrec);	

    #if DEBUGSOUT
    fclose(file_debug);
	#endif
    #if PRINTTHETAS
    fclose(file_debug_sel);
    #endif
    
    return 0;
}

void print_matrix(long int segsites,struct var2 **inputp,FILE *output,int x,long int count)
{
    long j,k;
    int i,h;
    double randoms(void);
    int ispolnomhit(long int,int,int);
     
    /*x=1 pr_matrix, x=2 print Hudson format (modified for mhits), x=3 print matrix excluding positions with mhits*/
    char **list2;
    
	if(x==3) fputs("\nMatrix of dna excluding all positions with multiple hits\n",output);
	/*
	if((*inputp)->segsitesin != -1) {
		fputs("\nlocus\tmig_r\trec\ttheta1\tdout",output);
		fprintf(output,"\n%d\t%.3G\t%.3G\t%.3G\t%.3G\t", (*inputp)->nloci,(*inputp)->migrate, (*inputp)->r, (*inputp)->theta, (*inputp)->T_out);
	}
	else {
		fputs("\nlocus\tmig_r\trec\tdout",output);
		fprintf(output,"\n%d\t%.3G\t%.3G\t%.3G\t", (*inputp)->nloci,(*inputp)->migrate, (*inputp)->r,(*inputp)->T_out);
	}
	*/
    if(x==2) {
        fprintf(output,"\n//");
        fprintf(output,"\nsegsites: %ld",segsites);
        if(segsites > 0) fprintf(output,"\npositions: ");
        /*for(i=0;i<segsites;i++) fprintf(output,"%ld ",posit[i]);*/
        for(i=0;i<(int)segsites;i++) fprintf(output,"%g ",(double)posit[i]/(double)(*inputp)->nsites);
        fprintf(output,"\n");
        if(segsites > 0) {
			for(i=0;i<(*inputp)->nsam;i++) list[i][segsites] = '\0';
			for(i=0;i<(*inputp)->nsam;i++) fprintf(output,"%s\n",list[i]);
		}
    }
    if(x==1) {
		/*fprintf(output,"\n%d %ld\n",(*inputp)->nsam,(*inputp)->nsites);*/
		fprintf(output,"\nFASTA file: locus %d, nsam %d, nsites %ld, mutations %ld, iteration %ld\n",(*inputp)->nloci,(*inputp)->nsam,(*inputp)->nsites,segsites,count);
        if(!(list2 = (char **)malloc((unsigned)((*inputp)->nsam)*sizeof(char *)))) perror("malloc error in ms.3a");
        for(i=0;i<(*inputp)->nsam;i++) {
            if(!(list2[i] = (char *)malloc((long int)((*inputp)->nsites+1)*sizeof(char))))
                perror("malloc error in ms.3");
            memset(list2[i],'A',(*inputp)->nsites);
            list2[i][(*inputp)->nsites] = '\0';
        }
        /*if(segsites > 0) {*/
            for(i=0;i<(*inputp)->nsam;i++) {
                if(list[i][0] == '1') list2[i][posit[0]] = 'G';
                if(list[i][0] == '2') list2[i][posit[0]] = 'C';
                if(list[i][0] == '3') list2[i][posit[0]] = 'T';
            }
            for(j=1;j<(int)segsites;j++) {
                if(!(posit[j-1] == posit[j])) {
                    for(i=0;i<(*inputp)->nsam;i++) {
                        if(list[i][j] == '1') list2[i][posit[j]] = 'G';
                        if(list[i][j] == '2') list2[i][posit[j]] = 'C';
                        if(list[i][j] == '3') list2[i][posit[j]] = 'T';
                    }
                }
            }
			/*
			if((*inputp)->ifselection==1) {
				if((*inputp)->sel_nt >= 0 && (*inputp)->sel_nt < (long int)(*inputp)->nsites) {
					if((*inputp)->sinit >= (double)0) {
						for(i=0;i<(*inputp)->nsam;i++) {
							if(list2[i][(*inputp)->sel_nt] == 'A') list2[i][(*inputp)->sel_nt] = 'G';
							if(list2[i][(*inputp)->sel_nt] == 'G') list2[i][(*inputp)->sel_nt] = 'C';
							if(list2[i][(*inputp)->sel_nt] == 'C') list2[i][(*inputp)->sel_nt] = 'T';
							if(list2[i][(*inputp)->sel_nt] == 'T') list2[i][(*inputp)->sel_nt] = 'A';
						}
					}
					else {
						j = 0;
						for(i=0;i<(*inputp)->nsam;i++) {
							if(list2[i][(*inputp)->sel_nt] != 'A') {
								j = 1;
								break;
							}
						}
						if(j==0) {
							for(i=0;i<(*inputp)->nsam;i++) {
								list2[i][(*inputp)->sel_nt] = 'G';
							}
						}
					}
				}
			}
			*/
            /*for(i=0;i<(*inputp)->nsam;i++) fprintf(output,"L%-9d%s\n",i,list2[i]);*/
			for(i=0;i<(*inputp)->nsam;i++) fprintf(output,">L%d\n%s\n",i,list2[i]);
        /*}*/
        fprintf(output,"\n");
        for(i=0;i<(*inputp)->nsam;i++) free(list2[i]);
        free(list2);
    }
    if(x==3) { /* exclude mhits from the matrix, but also the outgroup sequence in speciation with mhits */
        if(!(list2 = (char **)malloc((unsigned)((*inputp)->nsam)*sizeof(char *)))) perror("malloc error in ms.3a");
        for(i=0;i<(*inputp)->nsam;i++) {
            if(!(list2[i] = (char *)malloc((long int)((*inputp)->nsites+1)*sizeof(char))))
                perror("malloc error in ms.3");
            memset(list2[i],'A',(*inputp)->nsites);
            list2[i][(*inputp)->nsites] = '\0';
        }
        k = 0;
        if(segsites > 0) {
            if((h=ispolnomhit(0,0,(*inputp)->nsam)) > 0) {
                for(i=0;i<(*inputp)->nsam;i++) {
                    if(list[i][0] == '1') list2[i][posit[0]-k] = 'G';
                    if(list[i][0] == '2') list2[i][posit[0]-k] = 'C';
                    if(list[i][0] == '3') list2[i][posit[0]-k] = 'T';
                }
            }
            else if(h == -2) k++;/*k=0; check positions*/
            for(j=1;j<(int)segsites;j++) {
                if((h=ispolnomhit(j,0,(*inputp)->nsam)) > 0) {
                    for(i=0;i<(*inputp)->nsam;i++) {
                        if(list[i][j] == '1') list2[i][posit[j]-k] = 'G';
                        if(list[i][j] == '2') list2[i][posit[j]-k] = 'C';
                        if(list[i][j] == '3') list2[i][posit[j]-k] = 'T';
                    }
                }
                else if(h == -2) k++;/*k=0; check positions*/
            }
        }
		/*
		if((*inputp)->ifselection==1) {
			if((*inputp)->sel_nt >= 0 && (*inputp)->sel_nt < (long int)(*inputp)->nsites) {
				if((*inputp)->sinit >= (double)0) {
					for(i=0;i<(*inputp)->nsam;i++) {
						if(list2[i][(*inputp)->sel_nt] == 'A') list2[i][(*inputp)->sel_nt] = 'G';
						if(list2[i][(*inputp)->sel_nt] == 'G') list2[i][(*inputp)->sel_nt] = 'C';
						if(list2[i][(*inputp)->sel_nt] == 'C') list2[i][(*inputp)->sel_nt] = 'T';
						if(list2[i][(*inputp)->sel_nt] == 'T') list2[i][(*inputp)->sel_nt] = 'A';
					}
				}
				else {
					j = 0;
					for(i=0;i<(*inputp)->nsam;i++) {
						if(list2[i][(*inputp)->sel_nt] != 'A') {
							j = 1;
							break;
						}
					}
					if(j==0) {
						for(i=0;i<(*inputp)->nsam;i++) {
							list2[i][(*inputp)->sel_nt] = 'G';
						}
					}
				}
			}
		}
		*/
        /*if(segsites > 0) {*/
            /*fprintf(output,"\n%d %ld\n",(*inputp)->nsam,((*inputp)->nsites)-k);*/
			fprintf(output,"\nFASTA file: locus %d, nsam %d, nsites %ld, mutations %ld, iteration %ld\n",(*inputp)->nloci,(*inputp)->nsam,((*inputp)->nsites)-k,segsites,count);
            for(i=0;i<(*inputp)->nsam;i++) {
                list2[i][((*inputp)->nsites)-k] = '\0';
                /*fprintf(output,"L%-9d%s\n",i,list2[i]);*/
				fprintf(output,">L%d\n%s\n",i,list2[i]);
            }
        /*}*/
        fprintf(output,"\n");
        for(i=0;i<(*inputp)->nsam;i++) free(list2[i]);
        free(list2);
    }
}

void mod_mhits(long int segsites, struct var2 **inputp)
{
    long int x,y,z,nsit;
    int r,h,i,j,k;
    char a[1];
    int Sout;
    double rr,ratio,r_transv,r_transc;
    double ran1(void);
    char *mhsout;
	double poissondist(double);
	double p;      
    
	ratio = (*inputp)->ratio_sv;
    r_transc = ratio/(ratio + 1.);
    r_transv = (/*ratio + */0.5)/(ratio + 1.); /*suma de transc + 1nt a transversio, quan sumen els 2nt es total = 1*/
    Sout = (*inputp)->Sout;
    nsit = (*inputp)->nsites;
    x = 0;
	*a = '0';
				        
    if(segsites != 0) {
		while(x<(int)segsites-1) {
			k = 0;
			while(posit[x] == posit[x+1]) {	/* buscar els mhits */
				k++;
				x++;
				if(x == (long int)segsites-1) break;
			}
			if(k) {				/* ordenar mhits de mes antic a mes recent */
				for(y=(x-k);y<x;y++) {
					j = 0;
					for(i=0;i<(*inputp)->nsam;i++) if(list[i][y] != '0') j++;
					for(z=y+1;z<(x+1);z++) {
						h = 0;
						for(i=0;i<(*inputp)->nsam;i++) if(list[i][z] != '0') h++;
						if(j < h) {
							for(i=0;i<(*inputp)->nsam;i++) {
								*a = list[i][y];
								list[i][y] = list[i][z];
								list[i][z] = *a;
							}
							j = h;
						}
					}
				}
				if((*inputp)->segsitesin == -1 || ((*inputp)->segsitesin > 0 && (*inputp)->mhits == 1 && (*inputp)->Sfix_alltheta)) {/*in case no fixed mutations*/
					for(y=(x-k+1);y<x+1;y++) {	/* afegir les mutacions a la posicio mes antiga */
						r = -1;
						for(i=0;i<(*inputp)->nsam;i++) {
							if(list[i][y] != '0') {
								if(r == -1) {
									rr = (double)ran1();	/* inclou ratio trans/transv */
									if(rr < r_transc) r = 0;	/* transicio, la resta transversions */
									else if (rr >= r_transc && rr < (double)1.0 - r_transv) r = 1;
										else r = 2;
									if(list[i][x-k] == '0') {
										if(r==0) list[i][x-k] = *a = '1';
										else if(r==1) list[i][x-k] = *a = '2';
											else if(r==2) list[i][x-k] = *a = '3';
									}
									else
										if(list[i][x-k] == '1') {
											if(r==0) list[i][x-k] = *a = '0';
											else if(r==1) list[i][x-k] = *a = '2';
												else if(r==2) list[i][x-k] = *a = '3';
										}
										else
											if(list[i][x-k] == '2') {
												if(r==0) list[i][x-k] = *a = '3';
												else if(r==1) list[i][x-k] = *a = '0';
													else if(r==2) list[i][x-k] = *a = '1';
											}
											else
												if(list[i][x-k] == '3') {
													if(r==0) list[i][x-k] = *a = '2';
													else if(r==1) list[i][x-k] = *a = '1';
														else if(r==2) list[i][x-k] = *a = '0';
												}
								}
								else list[i][x-k] = *a;
							}
						}
					}
				}
				else {
					for(y=(x+1-k),r=2;y<x+1;y++,r++) {	/* mutacions fix, totes s'han de veure (fins a tres mutacions) */
						for(i=0;i<(*inputp)->nsam;i++) 
							if(list[i][y] != '0') 
								list[i][x-k] = r + '0';
					}
					if(r>4) {
						perror("Error: more than 3 mutations in one position");
						exit(1);
					}
				}
			}
			/* incorporacio per mhits en especiacio.
			 Els mhits al outgroup no s'observen, nomes canviem el nt ancestral, pero no estan indicades!*/
			if((nsit > 0) && ((p = (double)poissondist((double)Sout/(double)nsit)) > (double)0)) {/*nombre de mutacions de Sout a la pos x-k*/
				Sout -= (int)p;
				if(Sout < 0) Sout = 0;
				nsit--;
				/*p -= (double)1;*/
				
				rr = (double)ran1();
				if(rr < r_transc) r = 1;/* transicio, la resta transversions */
				else if (rr >= r_transc && rr < (double)1.0 - r_transv) r = 2;
					else r = 3;
				for(i=0;i<(*inputp)->nsam;i++) {
					if(list[i][x-k] == '0') list[i][x-k] = (char)r + '0';
					else if (list[i][x-k] == (char)r + '0') list[i][x-k] = '0';
				}
				/*cauen al mateix lloc pero no es veuen, nomes la ultima*/
				/*Sout -= (int)p;*/
				/*nsit -= (int)p;*/
			}
			x++;
		}
	}
    /* Treure els mhits de la branca de l'outgroup (no es mostra l'outgroup i per tant es veuen menys). */
    /* Per Sout mutacions en (nsites - segsites) */
    if((*inputp)->segsitesin == -1 || ((*inputp)->segsitesin > 0 && (*inputp)->mhits == 1 && (*inputp)->Sfix_alltheta)) {
        nsit = (*inputp)->nsites-segsites;
        if(nsit <= 0) nsit = 0;
        if(!(mhsout = (char *)calloc((long int)(nsit),sizeof(char)))) {
            perror("calloc error in mhits.0b");
            exit(1);
        }
        for(x=0;x<Sout;x++) {
            y = (long int)(ran1()*(double)nsit);
            rr = (double)ran1();
            if(rr < r_transc) r = 1;	/* transicio, la resta transversions */
            else if (rr >= r_transc && rr < (double)1.0 - r_transv) r = 2;
                else r = 3 ;
            if(mhsout[y] == 0) mhsout[y] = r;
            else if(mhsout[y] == 1) {
                    if(r==1) mhsout[y] = 0;
                    else mhsout[y] = r;
                }
                else if(mhsout[y] == 2) {
                        if(r==1) mhsout[y] = 3;
                        else mhsout[y] = r-2;
                    }else if(mhsout[y] == 3) {
                            if(r==1) mhsout[y] = 2;
                            else mhsout[y] = r-2;
                        }
        }
        Sout = 0;
        for(x=0;x<nsit;x++) if(mhsout[x] > 0) Sout++;
        free(mhsout);
    }
    (*inputp)->Sout = Sout;
}
char **cmatrix(int nsam,long int len)	/* defineix l'espai per col.locar els polimorfismes */
{
    int i;
    char **m;
    
    if(!(m=(char **)malloc((unsigned)nsam*sizeof(char *)))) perror("alloc error in cmatrix");
    for(i=0;i<nsam;i++)
        if(!(m[i] = (char *)malloc((long int)len*sizeof(char)))) perror("alloc error in cmatrix.2");
    return(m);
}

long int gensam(long int npop,int nsam,int inconfig[],long int nsites,double theta,long int segsites,
	double r,double f,double track_len,double mig_rate,int mhits,long int iteration, double *factor, double *lengtht,
	int ifselection, double pop_sel, double sinit,double pop_size,long int sel_nt,double T_out, int *Sout,int nintn,double *nrec,
	double *npast,double *tpast,int split_pop, double time_split, double time_scoal, double factor_anc, double *freq, 
	double tlimit,int iflogistic,double ts, double factor_chrn, double rsv,double *weightmut,double *weightrec)
{
	int i,ii;
    long int nsegs,seg,ns,start,end,len,segsit,k; 
    struct segl *seglst;
    double /*nsinv,*/tseg,tt;
    double *pk,tout,tout2;
    long int *ss;
    long int *len2;
    long int mmax;
	double r_transc,r_transv;
    struct segl *segtre_mig(long int ,int,int *,long int,double,double,double,
        double ,long int *,long int,
        double *,int,double,double,double,long int,int *,int *,int,double *,double *,double *,
        int, double, double, double, double *,double,int,double,double,double *);/* used to be: [MAXSEG]; */
    double ttime(struct node *, int);
    double poissondist(double), ran1(void);
    void biggerlist(int);
    void make_gametes(int,struct node *,double,long int,long int,int,double,double);
    void locate(long int,/*double*/long int,/*double*/ /*long int,*/ /*double **/long int *,int /*nou*/,double *,long int);
    void locate2(long int,/*double*/long int,/*double*/ /*long int,*/ /*double **/long int *,int /*nou*/,
        int/*,struct node *,double,long int*/,double *,long int);
    void mnmial2(long int,long int,double *,long int *,long int *);
    void mnmial2_psel(long int,long int,double *,long int *,long int *,long int);
    /*partial selection*/
    int all_sel,*selnsam; /*number of lines under selection and the vector with the lines (the first all_sel lines)*/
    int segsit_sel=0;/*parameter for partial selection*/    
    /*void make_gametes_psel(int,long int,int,int *);*/
    void locate_psel(long int,long int,/*long int,*/long int *,int,long int,double *,long int);
    void locate2_psel(long int,long int,/*long int,*/long int *,int,int/*,struct node *,double,long int*/,
		long int,double *,long int);
	double wstartm1,ttt;
	long int *len_nozero,nsites_nozero,nz;
	int loopcount;
    /*nsinv = 1./nsites;*/    
    
	if(!(selnsam = (int *)malloc((nsam)*sizeof(int)))) perror("malloc error sel. gensam.");
	all_sel = 0;

    seglst = segtre_mig(npop,nsam,inconfig,nsites,r,f,track_len,mig_rate,&nsegs,
        iteration,factor,ifselection,pop_sel,sinit,pop_size,sel_nt,&all_sel,selnsam,
        nintn, nrec, npast, tpast,split_pop,time_split,time_scoal,factor_anc,freq,tlimit,iflogistic,ts,factor_chrn,weightrec);
    
    r_transc = rsv/(rsv + (double)1.);
    r_transv = ((double)0.5)/(rsv + (double)1.); /*suma de transc + 1nt a transversio, quan sumen els 2nt es total = 1*/

	/*heterogeneity*/
	if(!(len_nozero = (long int *)malloc((nsegs)*sizeof(long int)))) perror("malloc error len. gensam.");
	nsites_nozero = 0;
	for(seg=0,k=0;k<nsegs;seg=seglst[seg].next,k++) {
		end = (k<nsegs-1 ? (long int)seglst[seglst[seg].next].beg -1 : nsites-1); /*next és l'index, beg és el punt físic */
		start = seglst[seg].beg;
		for(len_nozero[k]=0,nz=start;nz<=end;nz++) {
			if(nz==0) wstartm1 = (double)0;
			else wstartm1 = weightmut[nz-1];
			if(weightmut[nz]-wstartm1 != (double)0) {
				len_nozero[k] += 1;
				nsites_nozero += 1;
			}
		}
	}
	
	/*Include mutations*/
    if(segsites == -1) {
        ns = 0;
        *Sout = 0;
		*lengtht = (double)0;
        for(seg=0,k=0;k<nsegs;seg=seglst[seg].next,k++) {
			end = (k<nsegs-1 ? (long int)seglst[seglst[seg].next].beg -1 : nsites-1); /*next és l'index, beg és el punt físic */
			start = seglst[seg].beg;
			/* part de la theta que li correspon al segment *//*
			len = end - start + 1;
			tseg = (double)len*(theta/(double)nsites);*/
			/*including heterogeneity*/
			len = end - start + 1;
			if(start==0) wstartm1 = (double)0;
			else wstartm1 = weightmut[start-1];
			tseg = (double)weightmut[end] - wstartm1;
			tt = ttime(seglst[seg].ptree,nsam); /* Ttot pel segment, en funció de 4No respecte la pob1*/
			*lengtht += (double)tt*((double)len/(double)nsites); /*length in 4N generations, no matter the mutation rate*/
			if(mhits) {/*T es el temps de divergencia entre dos individus de les especies!!!!*/
				tout = T_out;/*T is considered as a fixed value, not a parameter*/
				tout -= (seglst[seg].ptree + 2*nsam-2)->time; /*substract the distance of the sample*/
				if(tout < (double)0) tout = (double)0;	/*Outgroup t can not accumulate negative mutations*/
				tout += -log((double)((double)1-ran1()));/*time after divergence, assuming equal No*/
				*Sout += (int)poissondist((double)(tseg*tout));	/* Sout needed to calculate hka and mhits */
			}
            loopcount = -1;
			do {
				segsit = (long int)poissondist((double)(tseg*tt));		/* nombre de mutacions al segment */
				if(segsit == (long int)0 && all_sel > (int)0 && all_sel < nsam && sel_nt >= (long int)start && sel_nt <= (long int)end)
					segsit = 1;/*we force the selective mut*/
				loopcount += 1;
				if(loopcount > 100) {
					printf("\nSorry, the length of the sequence is too short to include so much mutations: try mhits 1.\n");
					exit(1);
				}
			} while(segsit > len_nozero[k] && mhits == 0); /* mutacions discretes ...: afegit.. */
            if((segsit + ns) >= maxsites) {	/* refem la matriu dels polimorfismes */
                maxsites = segsit + ns + SITESINC;
                posit = (/*double **/long int *)realloc(posit,(long int)maxsites*sizeof(/*double*/long int)); 
                /* canvia mida del vector dels nombres dels polimorfismes */
                if(posit==NULL) perror("realloc error. gensam.1");
                biggerlist(nsam);	/* refem la llista dels polimorfismes */
            }
            /*partial selection*//*not debugged yet*/ 
            if(all_sel > (int)0 && all_sel < (int)nsam && sel_nt >= (long int)start && sel_nt <= (long int)end) {
                /*make_gametes_psel(nsam,ns,all_sel,selnsam);*//*locate the sel mut in the selnsam lines*/
                segsit_sel = 1;
            }
            make_gametes(nsam,seglst[seg].ptree,tt,segsit-segsit_sel,ns+segsit_sel,mhits,r_transc,r_transv);	/* posa les mutacions  a list*/
            free(seglst[seg].ptree);
            /*partial selection*/
            if(segsit_sel == 1) {
                locate_psel(segsit,start,/*len,*/posit+ns,mhits,sel_nt/*-start*/,weightmut,end);
				ii = ns;
				while((long int)posit[ii] != sel_nt) ii++; 
				for(i=0;i<nsam;i++) {
					list[i][ns] = list[i][ii];
					if(all_sel>i) list[i][ii] = '1';
					else list[i][ii] = '0';
				}
                segsit_sel = 0;
            }
            else locate(segsit,start,/*len,*/posit+ns,mhits,weightmut,end);/* posa el nombre de les mutacions a la matriu */
            ns += segsit;
        }
    }
    else {/*THE PARTIAL SELECTION WITH SEG. SITE ARE NOT WELL DEBUGGED YET IN THE Sfix METHOD*/
        /* en cas de S fix */
		pk = (double *)malloc((unsigned)(nsegs*sizeof(double)));
        ss = (long int *)malloc((unsigned)(nsegs*sizeof(long int)));
        len2 = (long int *)malloc((unsigned)(nsegs*sizeof(long int)));
        if((pk==NULL)||(ss==NULL)||(len2 ==NULL)) perror("malloc error. gensam.2");
		/*multiple hits for fixed mutations*/
		if(mhits) mmax = (3 < nsam ? (long int)3 : (long int)(nsam-1)); /* mhits per nsam petites, no accepta en nsam=1, i nomes 2 en nsam=3 */
		else mmax = (long int)1; /*if no multiple hits available*/
        /*set time  and nsites_nozero to zero*/
		tt = ttt = (double)0;
		if(mhits) tout = (double)0;
		/*calcular primer la mida total de tot l'arbre*/
		for(seg=0,k=0;k<nsegs;seg=seglst[seg].next,k++) {		
            end = (k<nsegs-1 ? (long int)seglst[seglst[seg].next].beg -1 : nsites-1);
            start = seglst[seg].beg;
            /*homogeneity. len2[k] afegit, maxim nombre de mutacions per segment *//*
			len = end - start + 1;
			len2[k] = len*mmax;	
			tseg = (double)len/(double)nsites;*/
			/*when it is not possible to fix the specified segsites*/
			if(nsites_nozero*mmax < segsites) {
				printf("\n*****WARNING: It is not possible to fix %ld mutations!!. Instead used 0 mutations.*****\n",segsites);
				segsites = 0;
			}
			
			len = end - start + 1;
			if(start==0) wstartm1 = (double)0;
			else wstartm1 = weightmut[start-1];
			tseg = (double)weightmut[end] - wstartm1;
			len2[k] = len_nozero[k]*mmax;
			/**/
			pk[k] = ttime(seglst[seg].ptree,nsam) * tseg;/*time per chromosome section (in function of mutational rate)*/
			tt += pk[k];
			ttt += pk[k]/tseg * (double)len/(double)nsites;/*time to be counted by lengtht (not in function of mutational rate)*/
			if(mhits) {		/* incorporacio per mhits */
				tout2 = T_out;
				tout2 -= (seglst[seg].ptree + 2*nsam-2)->time; /*substract the distance of the sample*/
				if(tout2 < (double)0) tout2 = (double)0;	/*Outgroup t can not accumulate negative mutations*/
				tout2 += -log((double)((double)1-ran1()));/*time after divergence, assuming equal No*/
				tout += tout2 * tseg;/*time to the outgroup (in function of mutational rate)*/
			}
        }
        *lengtht = (double)ttt; /*afegit per Sfix_allthetas*/
        for(k=0;k<nsegs;k++) pk[k] /= tt;	/* aleshores dividir el temps proporcionalment per situar les mutacions */
        if(mhits) {		/* incorporacio per mhits en especiacio */
            *Sout = (int)poissondist(/**/(double)(theta*tout)/**//*(double)segsites * ((double)tout/((double)tt))*/);
        }
		/*partial selection*/
		if(all_sel > (int)0 && all_sel < (int)nsam)
			mnmial2_psel((long int)segsites,nsegs,pk,ss,len2,sel_nt*mmax);
        else 
			mnmial2((long int)segsites,nsegs,pk,ss,len2);/* afegit, per evitar mes mutacions que posicions, en mhits mes de 3xpos */
        ns = 0;/*mnmial2 distribueix segsites al llarg de la secuencia*/
        for(seg=0,k=0;k<nsegs;seg=seglst[seg].next,k++) {
            end = (k<nsegs-1 ? (long int)seglst[seglst[seg].next].beg -1 : nsites-1);
            start = seglst[seg].beg;
            len = end - start + 1;
			if(start==0) wstartm1 = (double)0;
			else wstartm1 = weightmut[start-1];
			tseg = (double)weightmut[end] - wstartm1;
            /*partial selection*/
			segsit_sel = 0;
            if(all_sel > (int)0 && all_sel < (int)nsam && sel_nt >= (long int)start && sel_nt <= (long int)end && segsites > 0) {
                /*make_gametes_psel(nsam,ns,all_sel,selnsam);*//*locate the sel mut in the selnsam lines*/
                segsit_sel = 1;
            }
            make_gametes(nsam,seglst[seg].ptree,tt*pk[k]/tseg,ss[k]-segsit_sel,ns+segsit_sel,mhits,/*1.0*/r_transc,/*0.0*/r_transv);/*posa a la matriu list les mutacions*/
            /*partial selection*/
            if(segsit_sel == 1) {
                locate2_psel(ss[k],start,/*len,*/posit+ns,mhits,nsam/*,seglst[seg].ptree,tt*pk[k]/tseg,ns*/,sel_nt/*-start*/,weightmut,end);
				ii = ns;
				while((long int)posit[ii] != sel_nt) ii++; 
				for(i=0;i<nsam;i++) {
					list[i][ns] = list[i][ii];
					if(all_sel>i) list[i][ii] = '1';
					else list[i][ii] = '0';
				}
                segsit_sel = 0;
            }
            else locate2(ss[k],start,/*len,*/posit+ns,mhits,nsam/*,seglst[seg].ptree,tt*pk[k]/tseg,ns*/,weightmut,end); /* posa el nombre de les mutacions a la matriu */
            /* modificat per mhits i fix muts */
            free(seglst[seg].ptree);
            ns += ss[k];
        }
        free(pk);
        free(ss);
        free(len2);
    }
    free(selnsam);
	free(len_nozero);
    return(ns);
}
void biggerlist(int nsam)	/* fa més gran la matriu dels polimorfismes */
{
    int i;
    for(i=0;i<nsam;i++) {
        list[i] = (char *)realloc(list[i],maxsites*sizeof(char));
        if(list[i] == NULL) perror("realloc error. biggerlist");
    }
}
double ttime(struct node *ptree, int nsam)	/* la Ttot de l'arbre */
{
    double t;
    int i;
    
    t = (ptree + 2*nsam-2)->time;
    for(i=nsam;i<2*nsam-1;i++) t += (ptree + i)->time;
    return(t);
}
/*
void make_gametes_psel(int nsam,long int ns,int all_sel,int *selnsam, int mhits, double r_transc,double r_transv)
{
    int tip;
	double rr;
	double ran1();
	char r;
	
	if(mhits) {
		rr = (double)ran1();
		if(rr < r_transc) r = '1';
		else if (rr >= r_transc && rr < (double)1.0 - r_transv) r = '2';
			else r = '3';
	}
	else r='1';
	
    for(tip=0;tip<nsam;tip++)
        list[tip][ns] = '0';
    for(tip=0;tip<all_sel;tip++)
        list[selnsam[tip]][ns] = r; 
}
*/
void make_gametes(int nsam, struct node *ptree, double tt,long int newsites,long int ns, int mhits, double r_transc,double r_transv)
{					/* posa les mutacions a la matriu list */
    long int j;
    int tip,node;
    int pickb(int, struct node *,double);
    int tdesn(struct node *,int,int);
       
	double rr;
	double ran1();
	char r;
	

    for(j=ns;j<ns+newsites;j++) {
		if(mhits) {
			rr = (double)ran1();
			if(rr < r_transc) r = '1';/* transicio, la resta transversions */
			else if (rr >= r_transc && rr < (double)1.0 - r_transv) r = '2';
				else r = '3';
		}
		else r='1';
        node = pickb(nsam,ptree,tt);	/* busca una branca al'atzar en funcio de la mida de t*/
        for(tip=0;tip<nsam;tip++) {
            if(tdesn(ptree,tip,node))	{ /* posa mutació si la mostra té relació amb la branca */
                list[tip][j] = r; 
			}
            else
                list[tip][j] = '0';
        }
    }
}
int pickb(int nsam, struct node *ptree,double tt) /* agafa la branca a on ha caigut la mutació */
{
    double x,y,z;
    int i;
    double ran1(void);
    
    x = (double)ran1()*tt;
    for(i=0,y=0;i<2*nsam-2;i++) {
        z = (ptree + (ptree+i)->abv)->time - (ptree+i)->time;
        if(z < 0.) z = 0.; /* en cas d'extincio, per si hi han problemes de precissio */
        y += z;
        if(y >= x) return(i);
    }
    return(i);
}
int tdesn(struct node *ptree, int tip, int node) /* mira si la mostra que mirem esta relacionada amb la branca */
{
    int k;
    
    for(k=tip;k<node;k = (ptree+k)->abv);
    
    if(k==node) return(1);
    else return(0);
}
void mnmial(long int n,long int nclass,double *p,long int *rv)
{
    
    double x,s;
    long int i,j;
    double ran1(void);
    
    for(i=0;i<(long int)nclass;i++) rv[i]=0;	/* inicialitzar */
    for(i=0;i<(long int)n;i++) {	/* posa les n mutacions */
        x = (double)ran1();	/* agafa un nombre [0,1) */
        j = 0;
        s = p[0];	/* p és el temps total de cada segment relatiu a 1 */
        while((x>s) && (j<((long int)nclass-(long int)1))) s += p[++j]; /* busca el segment fins que ran sigui major que s, posa a j+1 */
        rv[j]++; 	/* afegeix la mutació al segment j */
    }
}
void mnmial2(long int n,long int nclass,double *p,long int *rv,long int *len)
{
    double x,s;
    long int i,j;
    double ran1(void);
    
    for(i=0;i<(long int)nclass;i++) rv[i]=0;	/* inicialitzar */
    i = 0;
    while(i<(long int)n) {	/* posa les n mutacions */
        x = (double)ran1();	/* agafa un nombre [0,1) */
        j = 0;
        s = p[0];	/* p és el temps total de cada segment relatiu a 1 */
        while((x>s) && (j<((long int)nclass-(long int)1))) s += p[++j]; /* busca el segment fins que ran sigui major que s, posa a j+1 */
        if(len[j] > 0) {
            len[j]--;	
            rv[j]++; 	/* afegeix la mutació al segment j */
            i++;
        }
    }
}
void mnmial2_psel(long int n,long int nclass,double *p,long int *rv,long int *len,long int sel_nt)
{
    double x,s;
    long int i,j;
    double ran1(void);
	long int sumlen;
    
    for(i=0;i<(long int)nclass;i++) rv[i]=0;	/* inicialitzar */
    j=0;/*locate the selected position*/
	sumlen = len[j];
	while(sumlen<sel_nt) {
		j++;
		sumlen +=len[j];
	}
	rv[j] += 1;
	len[j] -=1;
	i = 1;	
    while(i<(long int)n) {	/* posa les altres n-1 mutacions */
        x = (double)ran1();	/* agafa un nombre [0,1) */
        j = 0;
        s = p[0];	/* p és el temps total de cada segment relatiu a 1 */
        while((x>s) && (j<((long int)nclass-(long int)1))) s += p[++j]; /* busca el segment fins que ran sigui major que s, posa a j+1 */
        if(len[j] > 0) {
            len[j]--;	
            rv[j]++; 	/* afegeix la mutació al segment j */
            i++;
        }
    }
}
void locate_psel(long int n,long int beg,/*long int len,*/long int *ptr,int mhits,long int sel_nt,double *weightmut,long int end)
/* localitza les mutacions en un fragment */
{
     /*long int i;*/
     void ordran_psel(long int,long int *,/*long int,*/int,long int,long int,double *,long int);
     
     ordran_psel(n,ptr,/*len,*/mhits,sel_nt,beg,weightmut,end);	/* mutacions en [0,len) ordenades de major a menor */
     /*for(i=0;i<n;i++) ptr[i] = beg + ptr[i];*/	/* les escala entre inici i final del fragment */
}

void ordran_psel(long int n,long int *pbuf,/*long int len,*/int mhits,long int sel_nt,long int beg,double *weightmut,long int end)
{
    void ranvec_psel(long int,long int *,/*long int,*/int,long int,long int,double *,long int);
    void order(long int,long int *);

    ranvec_psel(n,pbuf,/*len,*/mhits,sel_nt,beg,weightmut,end);
    order(n,pbuf);
}
void ranvec_psel(long int n,long int *pbuf,/*long int len,*/int mhits,long int sel_nt,long int beg,double *weightmut,long int end)
 /* posa un nombre entre [0,len) */
{
    long int i,x;
    double ran1(void);
	double valuer,wstartm1;
	long int localize_positiontop(double *,double,long int,long int);
    
	/*include nsites*/
    pbuf[0] = (long int)sel_nt;
	if(beg==0) wstartm1 = (double)0;
	else wstartm1 = (double)weightmut[beg-1];
    for(i=1;i<(long int)n;i++) {
        valuer = ran1()*(double)(weightmut[end] - wstartm1) + wstartm1;
		pbuf[i] = localize_positiontop(weightmut,valuer,beg,end+1);
		if(pbuf[i] == pbuf[0]) {
			i--;
			continue;
		}   
        if(!mhits) {	/* bucle per no mhits (mhits=0) */
            x = i-1;
            while(x>=0) {
                if(pbuf[i] == pbuf[x]) {
                    valuer = ran1()*(double)(weightmut[end] - wstartm1) + wstartm1;
					pbuf[i] = localize_positiontop(weightmut,valuer,beg,end+1);
                    x = i-1;
                }
                else x--;
            }
        }
    }
}

void locate(long int n,/*double*/long int beg,/*double*/ /*long int len,*/ /*double*/long int *ptr,int mhits/*nou*/,double *weightmut,long int end)
/* localitza les mutacions en un fragment */
{
     /*long int i;*/
     void ordran(long int,/*double **/long int *,/*nou->*/ /*long int,*/int /*nou*/,long int,double *,long int);
     
     ordran(n,ptr,/*len,*/mhits,beg,weightmut,end);	/* mutacions en [0,len) ordenades de major a menor */
     /*for(i=0;i<n;i++) ptr[i] = beg + ptr[i]*//**len*//*;*//* les escala entre inici i final del fragment */
}

void ordran(long int n, /*double */long int *pbuf,/*nou->*/ /*long int len,*/int mhits/*nou*/,long int beg,double *weightmut,long int end)
{
    void ranvec(long int,/*double **/long int *,/*nou->*/ /*long int,*/int/*nou*/,long int,double *,long int);
    /* posa un nombre entre [0,len) */
    void order(long int,/*double **/long int *);/* ordena els valors */

    ranvec(n,pbuf,/*len,*/mhits,beg,weightmut,end);
    order(n,pbuf);
}

void ranvec(long int n, /*double */long int *pbuf,/*nou->*/ /*long int len,*/int mhits /*nou*/,long int beg,double *weightmut,long int end) /* posa un nombre entre [0,len) */
{
    long int i,x;
    double ran1(void);
	double valuer,wstartm1;
	long int localize_positiontop(double *,double,long int,long int);
    
    
	if(beg==0) wstartm1 = (double)0;
	else wstartm1 = (double)weightmut[beg-1];
    for(i=0;i<(long int)n;i++) {
        valuer = ran1()*(double)(weightmut[end] - wstartm1) + wstartm1;
		pbuf[i] = localize_positiontop(weightmut,valuer,beg,end+1);
        if(!mhits) {	/* bucle per no mhits (mhits=0) */
            x = i-1;
            while(x>=0) {
                if(pbuf[i] == pbuf[x]) {
                    valuer = ran1()*(double)(weightmut[end] - wstartm1) + wstartm1;
					pbuf[i] = localize_positiontop(weightmut,valuer,beg,end+1);
                    x = i-1;
                }
                else x--;
            }
        }
    }
}
void locate2(long int n,/*double*/long int beg,/*double*/ /*long int len,*/ /*double*/long int *ptr,int mhits/*nou*/,int nsam/*,struct node *ptree,double tt,long int ns*/,double *weightmut,long int end)/* localitza les mutacions en un fragment */
{
     /*long int i;*/
     void ordran2(long int,/*double **/long int *,/*nou->*/ /*long int,*/ int /*nou*/,int/*,struct node *,double,long int*/,long int,double *,long int);
     
     ordran2(n,ptr,/*len,*/mhits,nsam/*,ptree,tt,ns*/,beg,weightmut,end);	/* mutacions en [0,len) ordenades de major a menor */
     /*for(i=0;i<n;i++) ptr[i] = beg + ptr[i]*//**len*//*;*//* les escala entre inici i final del fragment */
}
void ordran2(long int n, /*double */long int *pbuf,/*nou->*/ /*long int len,*/int mhits/*nou*/,int nsam/*,struct node *ptree,double tt,long int ns*/,long int beg,double *weightmut,long int end)			
{
    void ranvec2(long int,/*double **/long int *,/*nou->*/ /*long int,*/int/*nou*/,int/*,struct node *,double*/,long int,double *,long int);
    /* posa un nombre entre [0,len) */
    void order(long int,/*double **/long int *);/* ordena els valors */
    /*void change_mut(long int,long int *,int,struct node *,double,long int);*/

    ranvec2(n,pbuf,/*len,*/mhits,nsam/*,ptree,tt*/,beg,weightmut,end);
    order(n,pbuf);
    /*if(mhits) change_mut(n,pbuf,nsam,ptree,tt,ns);*/
}
void ranvec2(long int n, /*double */long int *pbuf,/*nou->*/ /*long int len,*/int mhits /*nou*/,int nsam/*,struct node *ptree,double tt*/,long int beg,double *weightmut,long int end)	
/* posa un nombre entre [0,len) */
{
    long int i,x;
    int y;
    int z,f;
	
    double dlen;
    double a;
    double ran1(void);
	double valuer,wstartm1;
	long int localize_positiontop(double *,double,long int,long int);
    
	if(beg==0) wstartm1 = (double)0;
	else wstartm1 = (double)weightmut[beg-1];
    if(mhits) {
        dlen = (double)(weightmut[end] - wstartm1);
        for(i=0;i<(long int)n;i++) {
            a = (double)ran1()*dlen + (double)wstartm1;
            pbuf[i] = localize_positiontop(weightmut,(double)a,beg,end+1);   
            x = i-1;
            y = 1;
            while(x>=0) {
                if(pbuf[i] == pbuf[x]) {
                    y++;
                    if(y == 3 || y == nsam) {/*no mes de 4 nt per posicio*/
                        a = (double)ran1()*dlen + (double)wstartm1;
                        pbuf[i] = localize_positiontop(weightmut,(double)a,beg,end+1);   
                        x = i;
                        y = 1;
                    }
					f=0;
					for(z=0;z<nsam;z++) /*all mutations should be observed. equal pattern in different position*/
						if((list[z][i] == '0' && list[z][x] == '0') ||
						   (list[z][i] != '0' && list[z][x] != '0')) 
							f++;
					if(f==nsam) {
						a = (double)ran1()*dlen + (double)wstartm1;
						pbuf[i] = localize_positiontop(weightmut,(double)a,beg,end+1);   
						x = i;
						y = 1;
						break;
					}
                }
                x--;
            }
		}
    }
    else { /* per no mhits */
        for(i=0;i<(long int)n;i++) {
			valuer = ran1()*(double)(weightmut[end] - wstartm1) + wstartm1;
			pbuf[i] = localize_positiontop(weightmut,valuer,beg,end+1);
            x = i-1;
            while(x>=0) {
                if(pbuf[i] == pbuf[x]) {
                    valuer = ran1()*(double)(weightmut[end] - wstartm1) + wstartm1;
					pbuf[i] = localize_positiontop(weightmut,valuer,beg,end+1);
                    x = i-1;
                }
                else x--;
            }
        }
    }
}

void order(long int n, /*double */long int *pbuf)/* ordena els valors */
{
    long int gap,i,j;
    /*double*/long int temp;
    
    for(gap= n/2; gap>0;gap /= 2)
        for(i=gap;i<(long int)n;i++)
            for(j=i-gap;j>=0 && pbuf[j]>pbuf[j+gap];j -= gap) {
                temp = pbuf[j];
                pbuf[j] = pbuf[j+gap];
                pbuf[j+gap] = temp;
            }
}

void locate2_psel(long int n,long int beg,/*long int len,*/long int *ptr,int mhits,int nsam/*,struct node *ptree,double tt,long int ns*/,long int sel_nt,double *weightmut,long int end)/* localitza les mutacions en un fragment */
{
     /*long int i;*/
     void ordran2_psel(long int,long int *,/*long int,*/int,int/*,struct node *,double,long int*/,long int,long int,double *,long int);
     
     ordran2_psel(n,ptr,/*len,*/mhits,nsam/*,ptree,tt,ns*/,sel_nt,beg,weightmut,end);	/* mutacions en [0,len) ordenades de major a menor */
     /*for(i=0;i<n;i++) ptr[i] = beg + ptr[i];*//* les escala entre inici i final del fragment */
}
void ordran2_psel(long int n,long int *pbuf,/*long int len,*/int mhits,int nsam/*,struct node *ptree,double tt,long int ns*/,long int sel_nt,long int beg,double *weightmut,long int end)			
{
    void ranvec2_psel(long int,long int *,/*long int,*/int,int/*,struct node *,double*/,long int,long int,double *,long int);/* posa un nombre entre [0,len) */
    void order(long int,long int *);/* ordena els valors */
    /*void change_mut(long int,long int *,int,struct node *,double,long int);*/

    ranvec2_psel(n,pbuf,/*len,*/mhits,nsam/*,ptree,tt*/,sel_nt,beg,weightmut,end);
    order(n,pbuf);
    /*if(mhits) change_mut(n,pbuf,nsam,ptree,tt,ns);*/
}
void ranvec2_psel(long int n,long int *pbuf,/*long int len,*/int mhits,int nsam/*,struct node *ptree,double tt*/,long int sel_nt,long int beg,double *weightmut,long int end)	
/* posa un nombre entre [0,len) */
{
    long int i,x;
    int y;
    int z,f;
    double dlen;
    double a;
	double ran1(void);
	double valuer,wstartm1;
	long int localize_positiontop(double *,double,long int,long int);
    
	if(beg==0) wstartm1 = (double)0;
	else wstartm1 = (double)weightmut[beg-1];
    if(mhits) {
        dlen = (double)(weightmut[end] - wstartm1);
        for(i=0;i<(long int)n;i++) {
            if(i) {
				a = (double)ran1()*dlen + (double)wstartm1;
				pbuf[i] = localize_positiontop(weightmut,(double)a,beg,end+1);
			}
			else pbuf[0] = (long int)sel_nt;
			
            x = i-1;
            y = 1;
            while(x>=0) {
                if(pbuf[i] == pbuf[x]) {
                    y++;
                    if(y == 3 || y == nsam || x == 0) {/*no mes de 4 nt per posicio*/
                        a = (double)ran1()*dlen + (double)wstartm1;
                        pbuf[i] = localize_positiontop(weightmut,(double)a,beg,end+1);   
                        x = i;
                        y = 1;
                    }
					f=0;
					for(z=0;z<nsam;z++) /*all mutations should be observed. equal pattern in different position*/
						if(list[z][pbuf[i]] == list[z][pbuf[x]]) f++;
					if(f==nsam) {
						a = (double)ran1()*dlen + (double)wstartm1;
						pbuf[i] = localize_positiontop(weightmut,(double)a,beg,end+1);  
						x = i;
						y = 1;
						break;
					}
                }
                x--;
            }
		}
    }
    else { /* per no mhits */
		pbuf[0] = (long int)sel_nt;
        for(i=1;i<(long int)n;i++) {
			valuer = ran1()*(double)(weightmut[end] - wstartm1) + wstartm1;
			pbuf[i] = localize_positiontop(weightmut,valuer,beg,end+1);
            x = i-1;
            while(x>=0) {
                if(pbuf[i] == pbuf[x]) {
                    valuer = ran1()*(double)(weightmut[end] - wstartm1) + wstartm1;
					pbuf[i] = localize_positiontop(weightmut,valuer,beg,end+1);
                    x = i-1;
                }
                else x--;
            }
        }
    }
}

void order2(long int n, double **pbuf)/* ordena els valors */
{
    long int gap,i,j;
    double temp0;
    double temp1;

    for(gap= n/2; gap>0;gap /= 2)
        for(i=gap;i<(long int)n;i++)
            for(j=i-gap;j>=0 && pbuf[j][0]>pbuf[j+gap][0];j -= gap) {
                temp0 = pbuf[j][0];
                temp1 = pbuf[j][1];
                pbuf[j][0] = pbuf[j+gap][0];
                pbuf[j][1] = pbuf[j+gap][1];
                pbuf[j+gap][0] = temp0;
                pbuf[j+gap][1] = temp1;
            }
}

int print_neuttest(struct var **data,FILE *output,char *file_out)
{
    long int x;
    long int a,b,c;
    double sum,sum2;
    int numloc;
    int totalnloci=1;
    /*int l;*/
    FILE *outputind[NEUTVALUES2];
	char *el;
	char files[NEUTVALUES2+1][420] = {{"\0"},{"\0"},{"\0"},{"\0"},{"\0"},{"\0"},{"\0"},
                                    {"\0"},{"\0"},{"\0"},{"\0"},{"\0"},{"\0"},{"\0"},
                                    {"\0"},{"\0"},{"\0"},{"\0"},{"\0"},
                                    {"\0"},{"\0"},{"\0"},{"\0"},{"\0"},{"\0"},{"\0"},{"\0"},{"\0"},{"\0"},{"\0"},{"\0"},{"\0"}};
    char nfiles[NEUTVALUES2+1][20] = {{"_TD.out"},{"_Fs.out"},{"_FDn.out"},{"_FFn.out"},{"_FD.out"},{"_FF.out"},{"_H.out"},
                                    {"_B.out"},{"_Q.out"},{"_ZnA.out"},{"_Fst.out"},{"_Kw.out"},{"_Hw.out"},{"_R2.out"},
                                    {"_S.out"},{"_piw.out"},{"_pib.out"},{"_thetaWatt.out"},{"_ThetaTaj.out"},
                                    {"_ThetaFW.out"},{"_D_Dmin.out"},{"_H_norm.out"},{"_maxhap.out"},{"_maxhap.out"},{"_Rm.out"},
									{"_thetaFL.out"},{"_thetaL.out"},{"_ZengE.out"},{"_EW.out"},{"_Fstw.out"},{"_Pwh.out"},
									{"_PPercentiles.out"}};
    char names[NEUTVALUES2][20]  = {{"TD"},{"Fs"},{"FDn"},{"FFn"},{"FD"},{"FF"},{"H"},
                                    {"B"},{"Q"},{"ZnA"},{"Fst"},{"Kw"},{"Hw"},{"R2"},
                                    {"S"},{"piw"},{"pib"},{"thetaWatt"},{"ThetaTaj"},
                                    {"ThetaFW"},{"D_Dmin"},{"H_norm"},{"maxhap"},{"maxhap"},{"Rm"},
									{"ThetaFL"},{"thetaL"},{"ZengE"},{"EW"},{"Fstw"},{"Pwh"}};
    /*observed values*/
	FILE *filePvalue;
	double **matrix_avg;
	double **matrix_var;
	double **matrix_Pvalues;
	double **matrix_Pvaluesequal;
	double ***percentiles;
	int n;
	int observed=1;
	int pv=0;
	int compare_(const void *,const void *);
	long int count,countequal,total;
	double fcount;
	double *average,*variance;
	long int **validiter;
	FILE *duplik;
	char file_out1[420];
	
	c = NEUTVALUES2;
    if((*data)->linked == 1 && (*data)->despl > 0 && (*data)->window > 0) {
        totalnloci = (int)ceil(((double)(*data)->nsites[1] - (double)(*data)->window) / ((double)(*data)->despl)) + (int)1;
        /*for(l=(*data)->window;l<(*data)->nsites[1];l += (*data)->despl) totalnloci++;*/
	}	
    else {
        if((*data)->linked > 1) totalnloci = (*data)->linked;
        else totalnloci = (*data)->n_loci;
	}
	
    /*observed values*/
	/*
	for(a=0;a<c;a++) {
		if((*data)->obs_statistics[a][1] == 1) }
			observed = 1;
			break;
		}
	}
	*/
	if(observed == 1) {
		if(!(average = (double *)malloc((unsigned)(c*sizeof(double))))) {
			perror("calloc error ms.matrix_avg*");
			exit(1);
		}
		if(!(variance = (double *)malloc((unsigned)(c*sizeof(double))))) {
			perror("calloc error ms.matrix_var*");
			exit(1);
		}
		if(!(validiter = (long int **)malloc((unsigned)(c*sizeof(long int *))))) {
			perror("calloc error ms.matrix_avg*");
			exit(1);
		}
		if(!(matrix_avg = (double **)malloc((unsigned)(c*sizeof(double *))))) {
			perror("calloc error ms.matrix_avg*");
			exit(1);
		}
		if(!(matrix_var = (double **)malloc((unsigned)(c*sizeof(double *))))) {
			perror("calloc error ms.matrix_var*");
			exit(1);
		}
		if(!(matrix_Pvalues = (double **)malloc((unsigned)(c*sizeof(double *))))) {
			perror("calloc error ms.matrix_var*");
			exit(1);
		}
		if(!(matrix_Pvaluesequal = (double **)malloc((unsigned)(c*sizeof(double *))))) {
			perror("calloc error ms.matrix_var*");
			exit(1);
		}
		if(!(percentiles = (double ***)malloc((unsigned)(c*sizeof(double **))))) {
			perror("calloc error ms.matrix_var*");
			exit(1);
		}
		for(a=0;a<c;a++) {
			if(!(validiter[a] = (long int *)malloc((unsigned)((totalnloci+2)*sizeof(long int))))) {
				perror("calloc error ms.matrix_avg*");
				exit(1);
			}
			if(!(matrix_avg[a] = (double *)malloc((unsigned)((*data)->n_iter*sizeof(double))))) {
				perror("calloc error ms.matrix_avg*");
				exit(1);
			}
			if(!(matrix_var[a] = (double *)malloc((unsigned)((*data)->n_iter*sizeof(double))))) {
				perror("calloc error ms.matrix_var*");
				exit(1);
			}
			if(!(matrix_Pvalues[a] = (double *)malloc((unsigned)((totalnloci+2)*sizeof(double))))) {
				perror("calloc error ms.matrix_var*");
				exit(1);
			}
			if(!(matrix_Pvaluesequal[a] = (double *)malloc((unsigned)((totalnloci+2)*sizeof(double))))) {
				perror("calloc error ms.matrix_var*");
				exit(1);
			}
			if(!(percentiles[a] = (double **)malloc((unsigned)((totalnloci+2)*sizeof(double *))))) {
				perror("calloc error ms.matrix_var*");
				exit(1);
			}
			for(b=0;b<totalnloci+2;b++) {
				if(!(percentiles[a][b] = (double *)malloc((unsigned)(13*sizeof(double))))) {
					perror("calloc error ms.matrix_var*");
					exit(1);
				}
			}
		}
	}
        
    if((*data)->likelihood_line == 0) {
		if((*data)->print_neuttest == 2) {
			for(a=0;a<c;a++) {
				strcat(files[a],file_out);
				el = strrchr(files[a],'.');
				*el = '\0';
				strcat(files[a],nfiles[a]);
				if(!(outputind[a] = fopen(files[a],"w"))) return 1;
				for(b=0;b<totalnloci;b++) fprintf(outputind[a],"%s[%d]\t",names[a],(int)b);
				fputs("\n",outputind[a]);
			}
		}
		
		for(x=0;x<(*data)->n_iter;x++) {
			for(a=0;a<c;a++) {
				sum = sum2 = (double)0;
				numloc = 0;
				for(b=0;b<(long int)totalnloci;b++) {
					if((*data)->print_neuttest == 2) {
						if(matrix_test[(b*c)+a][x] == -10000) fputs("na\t",outputind[a]);
						else fprintf(outputind[a],"%.6g\t",matrix_test[(b*c)+a][x]);
					}
					if(matrix_test[(b*c)+a][x] != -10000) {
						sum  += matrix_test[(b*c)+a][x];
						sum2 += matrix_test[(b*c)+a][x] * matrix_test[(b*c)+a][x];
						numloc++;
					}
				}
				/*if((*data)->print_neuttest == 1) {*/
					if(numloc > 2) {
						sum  = sum/(double)numloc;
						sum2 = (sum2/((double)numloc) - sum*sum)*((double)numloc/((double)numloc-(double)1));
						/**/if(sum2 <= 1E-37) sum2 = (double)0;/**/
						fprintf(output,"%.6g\t%.6g\t",sum,sum2);
						
						/*observed values*/
						if(observed == 1) {
							matrix_avg[a][x] = sum;
							matrix_var[a][x] = sum2;
						}
					}
					else {
						if(numloc > 0) {
							sum  = sum/(double)numloc;
							fprintf(output,"%.6g\t",sum);
							if(totalnloci > 2) fprintf(output,"na\t");

							/*observed values*/
							if(observed == 1) {
								matrix_avg[a][x] = sum;
								matrix_var[a][x] = (double)-10000;
							}
						}
						else {
							fputs("na\t",output);
							if(totalnloci > 2) fputs("na\t",output);

							/*observed values*/
							if(observed == 1) {
								matrix_avg[a][x] = (double)-10000;
								matrix_var[a][x] = (double)-10000;
							}
						}
					}
				/*}*/
			}
			if((*data)->print_neuttest == 2) 
				for(a=0;a<c;a++) fputs("\n",outputind[a]);
			/*else*/ fputs("\n",output);
		}
		if((*data)->print_neuttest == 2) {
			for(a=0;a<c;a++) fclose(outputind[a]);
		}
	}	
	
	

	/*observed values*/
	if(observed == 1) {
		if((*data)->likelihood_line == 0) {
			printf("\n Calculating percentiles and (if required) probabilities for observed values... ");
			/*sort and calculate probabilities*/
			for(a=0;a<c;a++) {
				for(b=0;b<(long int)totalnloci;b++) {
					qsort(matrix_test[(b*c)+a],(long int)(*data)->n_iter,sizeof(double),compare_);
				}
				qsort(matrix_avg[a],(long int)(*data)->n_iter,sizeof(double),compare_);
				qsort(matrix_var[a],(long int)(*data)->n_iter,sizeof(double),compare_);
				/*calculate percentiles and probabilities*/
				for(b=0;b<(long int)totalnloci;b++) {
					count = countequal = total = (long int)0;
					for(x=0;x<(*data)->n_iter;x++) {
						if(matrix_test[(b*c)+a][x] != -10000  && (*data)->obs_statistics[a][b+2] != (double)-10000) {
							total++;
							if((*data)->obs_statistics[a][1] == 1) {
								if((*data)->obs_statistics[a][b+2] >= matrix_test[(b*c)+a][x] + 1e-05 ||
								   (*data)->obs_statistics[a][b+2] >= matrix_test[(b*c)+a][x] - 1e-05) 
								   count++;
							}
							if((*data)->obs_statistics[a][1] == 1) {
								if((*data)->obs_statistics[a][b+2] <= matrix_test[(b*c)+a][x] + 1e-05 &&
								   (*data)->obs_statistics[a][b+2] >= matrix_test[(b*c)+a][x] - 1e-05) 
								   countequal++;
							}
						}
					}
					/*probabilitites*/
					if(total) {
						matrix_Pvalues[a][b] = (double)count/(double)total;
						matrix_Pvaluesequal[a][b] = (double)countequal/(double)total;
					}
					else {
						matrix_Pvalues[a][b] = (double)-10000;
						matrix_Pvaluesequal[a][b] = (double)-10000;
					}
					validiter[a][b] = (long int)total;
					/*percentiles*/
					for(n=0;n<13;n++) {
						switch(n) {
							case 0:
								fcount = (double)total * (double)0.001;
								break;
							case 1:
								fcount = (double)total * (double)0.010;
								break;
							case 2:
								fcount = (double)total * (double)0.025;
								break;
							case 3:
								fcount = (double)total * (double)0.050;
								break;
							case 4:
								fcount = (double)total * (double)0.100;
								break;
							case 5:
								fcount = (double)total * (double)0.250;
								break;
							case 6:
								fcount = (double)total * (double)0.500;
								break;
							case 7:
								fcount = (double)total * (double)0.750;
								break;
							case 8:
								fcount = (double)total * (double)0.900;
								break;
							case 9:
								fcount = (double)total * (double)0.950;
								break;
							case 10:
								fcount = (double)total * (double)0.975;
								break;
							case 11:
								fcount = (double)total * (double)0.990;
								break;
							case 12:
								fcount = (double)total * (double)0.999;
								break;
						}
						if(fcount < (double)1 || fcount > (double)(total-1)) percentiles[a][b][n] = -10000;
						else {
							count=(long int)floor((double)fcount);
							if(fcount == (double)floor((double)fcount)) percentiles[a][b][n] = matrix_test[(b*c)+a][(*data)->n_iter-total+count-1];
							else percentiles[a][b][n] = (matrix_test[(b*c)+a][(*data)->n_iter-total+count-1] + matrix_test[(b*c)+a][(*data)->n_iter-total+count])/(double)2;
						}
					}
				}
				/*average*/
				count = countequal = total = (long int)0;
				numloc = 0;
				average[a] = (double)0;
				for(b=0;b<(long int)totalnloci;b++) {
					if((*data)->obs_statistics[a][b+2] != (double)-10000) {
						average[a] += (*data)->obs_statistics[a][b+2];
						numloc++;
					}
				}
				if(numloc) average[a] /= (double)numloc;
				else average[a] = (double)-10000;
				for(x=0;x<(*data)->n_iter;x++) {
					if(matrix_avg[a][x] != -10000 && average[a] != (double)-10000) {
						total++;
						if((*data)->obs_statistics[a][1] == 1) {
							if(average[a] >= matrix_avg[a][x] + 1e-05 ||
							   average[a] >= matrix_avg[a][x] - 1e-05) 
							   count++;
						}
						if((*data)->obs_statistics[a][1] == 1) {
							if(average[a] <= matrix_avg[a][x] + 1e-05 &&
							   average[a] >= matrix_avg[a][x] - 1e-05) 
							   countequal++;
						}
					}
				}
				/*probabilitites*/
				if(total) {
					matrix_Pvalues[a][totalnloci] = (double)count/(double)total;
					matrix_Pvaluesequal[a][totalnloci] = (double)countequal/(double)total;
				}
				else {
					matrix_Pvalues[a][totalnloci] = (double)-10000;
					matrix_Pvaluesequal[a][totalnloci] = (double)-10000;
				}
				validiter[a][totalnloci] = (long int)total;
				/*percentiles*/
				for(n=0;n<13;n++) {
					switch(n) {
						case 0:
							fcount = (double)total * (double)0.001;
							break;
						case 1:
							fcount = (double)total * (double)0.010;
							break;
						case 2:
							fcount = (double)total * (double)0.025;
							break;
						case 3:
							fcount = (double)total * (double)0.050;
							break;
						case 4:
							fcount = (double)total * (double)0.100;
							break;
						case 5:
							fcount = (double)total * (double)0.250;
							break;
						case 6:
							fcount = (double)total * (double)0.500;
							break;
						case 7:
							fcount = (double)total * (double)0.750;
							break;
						case 8:
							fcount = (double)total * (double)0.900;
							break;
						case 9:
							fcount = (double)total * (double)0.950;
							break;
						case 10:
							fcount = (double)total * (double)0.975;
							break;
						case 11:
							fcount = (double)total * (double)0.990;
							break;
						case 12:
							fcount = (double)total * (double)0.999;
							break;
					}
					if(fcount < (double)1 || fcount > (double)(total-1)) percentiles[a][totalnloci][n] = -10000;
					else {
						count=(long int)floor((double)fcount);
						if(fcount == (double)floor((double)fcount)) percentiles[a][totalnloci][n] = matrix_avg[a][(*data)->n_iter-total+count-1];
						else percentiles[a][totalnloci][n] = (matrix_avg[a][(*data)->n_iter-total+count-1] + matrix_avg[a][(*data)->n_iter-total+count])/(double)2;
					}
				}
				/*variance*/
				count = countequal = total = (long int)0;
				numloc = 0;
				average[a] = variance[a] = (double)0;
				for(b=0;b<(long int)totalnloci;b++) {
					if((*data)->obs_statistics[a][b+2] != (double)-10000) {
						average[a] += (*data)->obs_statistics[a][b+2];
						variance[a] += (*data)->obs_statistics[a][b+2] * (*data)->obs_statistics[a][b+2] ;
						numloc++;
					}
				}
				if(numloc>2) {
					average[a] /= (double)numloc;
					variance[a] = (variance[a]/((double)numloc) - average[a]*average[a])*((double)numloc/((double)numloc-(double)1));
					/**/if(variance[a] <= 1E-37) variance[a] = (double)0;/**/
				}
				else {
					variance[a] = (double)-10000;
					if(numloc==0) average[a] = (double)-10000;
				}
				for(x=0;x<(*data)->n_iter;x++) {
					if(matrix_var[a][x] != -10000 && variance[a] != (double)-10000) {
						total++;
						if((*data)->obs_statistics[a][1] == 1) {
							if(variance[a] >= matrix_var[a][x] + 1e-05 ||
							   variance[a] >= matrix_var[a][x] - 1e-05) 
							   count++;
						}
						if((*data)->obs_statistics[a][1] == 1) {
							if(variance[a] <= matrix_var[a][x] + 1e-05 &&
							   variance[a] >= matrix_var[a][x] - 1e-05) 
							   countequal++;
						}
					}
				}
				/*probabilitites*/
				if(total) {
					matrix_Pvalues[a][totalnloci+1] = (double)count/(double)total;
					matrix_Pvaluesequal[a][totalnloci+1] = (double)countequal/(double)total;
				}
				else {
					matrix_Pvalues[a][totalnloci+1] = (double)-10000;
					matrix_Pvaluesequal[a][totalnloci+1] = (double)-10000;
				}
				validiter[a][totalnloci+1] = (long int)total;
				/*percentiles*/
				for(n=0;n<13;n++) {
					switch(n) {
						case 0:
							fcount = (double)total * (double)0.001;
							break;
						case 1:
							fcount = (double)total * (double)0.010;
							break;
						case 2:
							fcount = (double)total * (double)0.025;
							break;
						case 3:
							fcount = (double)total * (double)0.050;
							break;
						case 4:
							fcount = (double)total * (double)0.100;
							break;
						case 5:
							fcount = (double)total * (double)0.250;
							break;
						case 6:
							fcount = (double)total * (double)0.500;
							break;
						case 7:
							fcount = (double)total * (double)0.750;
							break;
						case 8:
							fcount = (double)total * (double)0.900;
							break;
						case 9:
							fcount = (double)total * (double)0.950;
							break;
						case 10:
							fcount = (double)total * (double)0.975;
							break;
						case 11:
							fcount = (double)total * (double)0.990;
							break;
						case 12:
							fcount = (double)total * (double)0.999;
							break;
					}
					if(fcount < (double)1 || fcount > (double)(total-1)) percentiles[a][totalnloci+1][n] = -10000;
					else {
						count=(long int)floor((double)fcount);
						if(fcount == (double)floor((double)fcount)) percentiles[a][totalnloci+1][n] = matrix_var[a][(*data)->n_iter-total+count-1];
						else percentiles[a][totalnloci+1][n] = (matrix_var[a][(*data)->n_iter-total+count-1] + matrix_var[a][(*data)->n_iter-total+count])/(double)2;
					}
				}
			}	
		}
		else {
			if((*data)->n_iter > 1) {
				printf("\n Calculating likelihoods (2*log(P)) for observed values... ");
				fflush(stdout);
				/*create a new file*/
				file_out1[0] = '\0';
				strncat(file_out1,file_out,420);
				el = strrchr(file_out1,'.');
				*el = '\0';
				strncat(file_out1,"_1.out\0",420);
				if(!(duplik = fopen(file_out1,"w"))) return 1;
				/*init*/
				matrix_Pvaluesequal[0][totalnloci+1] = (double)0;
				/*loop*/
				for(a=0;a<c;a++) {
					matrix_Pvaluesequal[a][totalnloci] = (double)0;
					if((*data)->obs_statistics[a][1] == 1) {
						for(b=0;b<(long int)totalnloci;b++) {
							countequal = total = (long int)0;
							for(x=0;x<(*data)->n_iter;x++) {
								if(matrix_test[(b*c)+a][x] != -10000  && (*data)->obs_statistics[a][b+2] != (double)-10000) {
									total++;
									if((*data)->obs_statistics[a][b+2] - (*data)->likelihood_error[a] <= matrix_test[(b*c)+a][x]  &&
									   (*data)->obs_statistics[a][b+2] + (*data)->likelihood_error[a] >= matrix_test[(b*c)+a][x]) 
									   countequal++;
								}
							}
							/*probabilitites for each locus, each stistic*/
							if(total) {
								if(countequal) matrix_Pvaluesequal[a][b] = (double)2.*(double)log((double)countequal/(double)total);
								else matrix_Pvaluesequal[a][b] = (double)2.*(double)log((double)1/(double)((*data)->n_iter+1));
							}
							else {
								matrix_Pvaluesequal[a][b] = (double)2.*(double)log((double)1/(double)((*data)->n_iter+1));
							}
							/*probabilitites for all loci, each statistic*/
							matrix_Pvaluesequal[a][totalnloci] += matrix_Pvaluesequal[a][b];
							validiter[a][b] = (long int)total;
						}
						/*probabilitites for all loci, all statistics*/
						matrix_Pvaluesequal[0][totalnloci+1] += matrix_Pvaluesequal[a][totalnloci];
						validiter[a][totalnloci] = (long int)total;
					}
				}
				/*PRINT LIKELIHOODS IN A SINGLE LINE in the OUTPUT FILE and in another FILE but FOR A SIGLE LINE*/
				fprintf(output,"%.6g",matrix_Pvaluesequal[0][totalnloci+1]);
				fprintf(duplik,"%.6g",matrix_Pvaluesequal[0][totalnloci+1]);
				for(a=0;a<c;a++) {
					if((*data)->obs_statistics[a][1] == 1) {
						fprintf(output,"\t%.6g",matrix_Pvaluesequal[a][totalnloci]);
						fprintf(duplik,"\t%.6g",matrix_Pvaluesequal[a][totalnloci]);
					}
				}
				for(a=0;a<c;a++) {
					if((*data)->obs_statistics[a][1] == 1) {
						for(b=0;b<(long int)totalnloci;b++) {
							fprintf(output,"\t%.6g",matrix_Pvaluesequal[a][b]);
							fprintf(duplik,"\t%.6g",matrix_Pvaluesequal[a][b]);
						}
					}
				}
				fprintf(output,"\n");
				fclose(duplik);
			}
			else {
				/*FOR NITER==1 PRINT THE OBSERVED VALUES FOR EACH LOCUS IN A SINGLE LINE*/
				for(a=0;a<c;a++) {
					if((*data)->obs_statistics[a][1] == 1) {
						for(b=0;b<(long int)totalnloci;b++) {
							if(matrix_test[(b*c)+a][x] != -10000) 
								fprintf(output,"%.6g\t",matrix_test[(b*c)+a][x]);
							else
								fprintf(output,"na\t");
						}
					}
				}
				fprintf(output,"\n");
			}
		}
		
		if((*data)->likelihood_line == 0) {
			/*printing values in filePvalue*/
			strcat(files[NEUTVALUES2],file_out);
			el = strrchr(files[NEUTVALUES2],'.');
			*el = '\0';
			strcat(files[NEUTVALUES2],nfiles[NEUTVALUES2]);
			if(!(filePvalue = fopen(files[NEUTVALUES2],"w"))) return 1;
			
			for(a=0;a<c;a++) {
				if((*data)->obs_statistics[a][1] == 1) {
					pv = 1;
					break;
				}
			}			
			if(pv==1) {
				fprintf(filePvalue,"OBSERVED VALUES\n");
				fprintf(filePvalue,"\tTD\tFs\tFD*\tFF*\tFD\tFF\tH\tB\tQ\tZA\tFst\tKw\tHw\tR2\tS\tpi_w\tpi_b\tthetaWatt\tthetaTaj\tthetaFW\tD/Dmin\tHnorm\tmaxhap\tmaxhap1\tRm\tThetaFL\tThetaL\tZengE\tEW\tFstw\tPwh");
				for(b=0;b<(long int)totalnloci;b++) {
					fprintf(filePvalue,"\nlocus_%ld",b);
					for(a=0;a<c;a++) {
						if((*data)->obs_statistics[a][1] == 1) {
							/*if((double)(*data)->obs_statistics[a][b+2] == (double)-10000) fprintf(filePvalue,"\tna");
							else */fprintf(filePvalue,"\t%.6g",(*data)->obs_statistics[a][b+2]);
						}
						else fprintf(filePvalue,"\tx");
					}
				}
				fprintf(filePvalue,"\naverage");
				for(a=0;a<c;a++) {
					if((*data)->obs_statistics[a][1] == 1) {
						if(average[a] != (double)-10000) fprintf(filePvalue,"\t%.6g",average[a]);
						else fprintf(filePvalue,"\tna");
					}
					else fprintf(filePvalue,"\tx");
				}
				fprintf(filePvalue,"\nvariance");
				for(a=0;a<c;a++) {
					if((*data)->obs_statistics[a][1] == 1) {
						if(variance[a] != (double)-10000) fprintf(filePvalue,"\t%.6g",variance[a]);
						else fprintf(filePvalue,"\tna");
					}
					else fprintf(filePvalue,"\tx");
				}
				fprintf(filePvalue,"\n\nPROBABILITY THAT SIMULATED VALUES BE SMALLER OR EQUAL THAN THE OBSERVED VALUE. P(Sim <= Obs): ACCURACY OF +- 1e-6\n");
				fprintf(filePvalue,"\tTD\tFs\tFD*\tFF*\tFD\tFF\tH\tB\tQ\tZA\tFst\tKw\tHw\tR2\tS\tpi_w\tpi_b\tthetaWatt\tthetaTaj\tthetaFW\tD/Dmin\tHnorm\tmaxhap\tmaxhap1\tRm\tThetaFL\tThetaL\tZengE\tEW\tFstw\tPwh");

				for(b=0;b<(long int)totalnloci;b++) {
					fprintf(filePvalue,"\nlocus_%ld",b);
					for(a=0;a<c;a++) {
						if((*data)->obs_statistics[a][1] == 1) {
							if(matrix_Pvalues[a][b] != (double)-10000) fprintf(filePvalue,"\t%.6g",matrix_Pvalues[a][b]);
							else fprintf(filePvalue,"\tna");
						}
						else fprintf(filePvalue,"\tna");
					}
				}
				fprintf(filePvalue,"\naverage");
				for(a=0;a<c;a++) {
					if((*data)->obs_statistics[a][1] == 1) {
						if(matrix_Pvalues[a][totalnloci] != (double)-10000) fprintf(filePvalue,"\t%.6g",matrix_Pvalues[a][totalnloci]);
						else fprintf(filePvalue,"\tna");
					}
					else fprintf(filePvalue,"\tna");
				}
				fprintf(filePvalue,"\nvariance");
				for(a=0;a<c;a++) {
					if((*data)->obs_statistics[a][1] == 1) {
						if(matrix_Pvalues[a][totalnloci+1] != (double)-10000) fprintf(filePvalue,"\t%.6g",matrix_Pvalues[a][totalnloci+1]);
						else fprintf(filePvalue,"\tna");
					}
					else fprintf(filePvalue,"\tna");
				}

				fprintf(filePvalue,"\n\nPROBABILITY THAT SIMULATED VALUES BE EQUAL THAN THE OBSERVED VALUE. P(Sim = Obs): ACCURACY OF +- 1e-6\n");
				fprintf(filePvalue,"\tTD\tFs\tFD*\tFF*\tFD\tFF\tH\tB\tQ\tZA\tFst\tKw\tHw\tR2\tS\tpi_w\tpi_b\tthetaWatt\tthetaTaj\tthetaFW\tD/Dmin\tHnorm\tmaxhap\tmaxhap1\tRm\tThetaFL\tThetaL\tZengE\tEW\tFstw\tPwh");
				for(b=0;b<(long int)totalnloci;b++) {
					fprintf(filePvalue,"\nlocus_%ld",b);
					for(a=0;a<c;a++) {
						if((*data)->obs_statistics[a][1] == 1) {
							if(matrix_Pvaluesequal[a][b] != (double)-10000) fprintf(filePvalue,"\t%.6g",matrix_Pvaluesequal[a][b]);
							else fprintf(filePvalue,"\tna");
						}
						else fprintf(filePvalue,"\tna");
					}
				}
				fprintf(filePvalue,"\naverage");
				for(a=0;a<c;a++) {
					if((*data)->obs_statistics[a][1] == 1) {
						if(matrix_Pvaluesequal[a][totalnloci] != (double)-10000) fprintf(filePvalue,"\t%.6g",matrix_Pvaluesequal[a][totalnloci]);
						else fprintf(filePvalue,"\tna");
					}
					else fprintf(filePvalue,"\tna");
				}
				fprintf(filePvalue,"\nvariance");
				for(a=0;a<c;a++) {
					if((*data)->obs_statistics[a][1] == 1) {
						if(matrix_Pvaluesequal[a][totalnloci+1] != (double)-10000) fprintf(filePvalue,"\t%.6g",matrix_Pvaluesequal[a][totalnloci+1]);
						else fprintf(filePvalue,"\tna");
					}
					else fprintf(filePvalue,"\tna");
				}
			}
			fprintf(filePvalue,"\n\nNumber of valid iterations:\n");
			fprintf(filePvalue,"\tTD\tFs\tFD*\tFF*\tFD\tFF\tH\tB\tQ\tZA\tFst\tKw\tHw\tR2\tS\tpi_w\tpi_b\tthetaWatt\tthetaTaj\tthetaFW\tD/Dmin\tHnorm\tmaxhap\tmaxhap1\tRm\tThetaFL\tThetaL\tZengE\tEW\tFstw\tPwh");
			for(b=0;b<(long int)totalnloci;b++) {
				fprintf(filePvalue,"\nlocus_%ld",b);
				for(a=0;a<c;a++) {
					fprintf(filePvalue,"\t%ld",validiter[a][b]);
				}
			}
			fprintf(filePvalue,"\naverage");
			for(a=0;a<c;a++) {
				fprintf(filePvalue,"\t%ld",validiter[a][totalnloci]);
			}
			fprintf(filePvalue,"\nvariance");
			for(a=0;a<c;a++) {
				fprintf(filePvalue,"\t%ld",validiter[a][totalnloci+1]);
			}

			fprintf(filePvalue,"\n\nPERCENTILES:");
			for(b=0;b<(long int)totalnloci;b++) {
				fprintf(filePvalue,"\n\nlocus_%ld\n\n",b);
				fprintf(filePvalue,"\tTD\tFs\tFD*\tFF*\tFD\tFF\tH\tB\tQ\tZA\tFst\tKw\tHw\tR2\tS\tpi_w\tpi_b\tthetaWatt\tthetaTaj\tthetaFW\tD/Dmin\tHnorm\tmaxhap\tmaxhap1\tRm\tThetaFL\tThetaL\tZengE\tEW\tFstw\tPwh");
				for(n=0;n<13;n++) {
					switch(n) {
						case 0:
							fprintf(filePvalue,"\n 0.1%%");
							break;
						case 1:
							fprintf(filePvalue,"\n 1.0%%");
							break;
						case 2:
							fprintf(filePvalue,"\n 2.5%%");
							break;
						case 3:
							fprintf(filePvalue,"\n 5.0%%");
							break;
						case 4:
							fprintf(filePvalue,"\n10.0%%");
							break;
						case 5:
							fprintf(filePvalue,"\n25.0%%");
							break;
						case 6:
							fprintf(filePvalue,"\n50.0%%");
							break;
						case 7:
							fprintf(filePvalue,"\n75.0%%");
							break;
						case 8:
							fprintf(filePvalue,"\n90.0%%");
							break;
						case 9:
							fprintf(filePvalue,"\n95.0%%");
							break;
						case 10:
							fprintf(filePvalue,"\n97.5%%");
							break;
						case 11:
							fprintf(filePvalue,"\n99.0%%");
							break;
						case 12:
							fprintf(filePvalue,"\n99.9%%");
							break;
					}
					for(a=0;a<c;a++) {
						if(percentiles[a][b][n] != (double)-10000) fprintf(filePvalue,"\t%.6g",percentiles[a][b][n]);
						else fprintf(filePvalue,"\tna");
					}
				}
			}
			/*average and variance*/
			fprintf(filePvalue,"\n\naverage\n\n");
			fprintf(filePvalue,"\tTD\tFs\tFD*\tFF*\tFD\tFF\tH\tB\tQ\tZA\tFst\tKw\tHw\tR2\tS\tpi_w\tpi_b\tthetaWatt\tthetaTaj\tthetaFW\tD/Dmin\tHnorm\tmaxhap\tmaxhap1\tRm\tThetaFL\tThetaL\tZengE\tEW\tFstw\tPwh");
			for(n=0;n<13;n++) {
				switch(n) {
					case 0:
						fprintf(filePvalue,"\n 0.1%%");
						break;
					case 1:
						fprintf(filePvalue,"\n 1.0%%");
						break;
					case 2:
						fprintf(filePvalue,"\n 2.5%%");
						break;
					case 3:
						fprintf(filePvalue,"\n 5.0%%");
						break;
					case 4:
						fprintf(filePvalue,"\n10.0%%");
						break;
					case 5:
						fprintf(filePvalue,"\n25.0%%");
						break;
					case 6:
						fprintf(filePvalue,"\n50.0%%");
						break;
					case 7:
						fprintf(filePvalue,"\n75.0%%");
						break;
					case 8:
						fprintf(filePvalue,"\n90.0%%");
						break;
					case 9:
						fprintf(filePvalue,"\n95.0%%");
						break;
					case 10:
						fprintf(filePvalue,"\n97.5%%");
						break;
					case 11:
						fprintf(filePvalue,"\n99.0%%");
						break;
					case 12:
						fprintf(filePvalue,"\n99.9%%");
						break;
				}
				for(a=0;a<c;a++) {
					if(percentiles[a][totalnloci][n] != (double)-10000) fprintf(filePvalue,"\t%.6g",percentiles[a][totalnloci][n]);
					else fprintf(filePvalue,"\tna");
				}
			}
			fprintf(filePvalue,"\n\nvariance\n\n");
			fprintf(filePvalue,"\tTD\tFs\tFD*\tFF*\tFD\tFF\tH\tB\tQ\tZA\tFst\tKw\tHw\tR2\tS\tpi_w\tpi_b\tthetaWatt\tthetaTaj\tthetaFW\tD/Dmin\tHnorm\tmaxhap\tmaxhap1\tRm\tThetaFL\tThetaL\tZengE\tEW\tFstw\tPwh");
			for(n=0;n<13;n++) {
				switch(n) {
					case 0:
						fprintf(filePvalue,"\n 0.1%%");
						break;
					case 1:
						fprintf(filePvalue,"\n 1.0%%");
						break;
					case 2:
						fprintf(filePvalue,"\n 2.5%%");
						break;
					case 3:
						fprintf(filePvalue,"\n 5.0%%");
						break;
					case 4:
						fprintf(filePvalue,"\n10.0%%");
						break;
					case 5:
						fprintf(filePvalue,"\n25.0%%");
						break;
					case 6:
						fprintf(filePvalue,"\n50.0%%");
						break;
					case 7:
						fprintf(filePvalue,"\n75.0%%");
						break;
					case 8:
						fprintf(filePvalue,"\n90.0%%");
						break;
					case 9:
						fprintf(filePvalue,"\n95.0%%");
						break;
					case 10:
						fprintf(filePvalue,"\n97.5%%");
						break;
					case 11:
						fprintf(filePvalue,"\n99.0%%");
						break;
					case 12:
						fprintf(filePvalue,"\n99.9%%");
						break;
				}
				for(a=0;a<c;a++) {
					if(percentiles[a][totalnloci+1][n] != (double)-10000) fprintf(filePvalue,"\t%.6g",percentiles[a][totalnloci+1][n]);
					else fprintf(filePvalue,"\tna");
				}
			}
		}

		for(a=0;a<c;a++) {
			free(validiter[a]);
			free(matrix_avg[a]);
			free(matrix_var[a]);
			free(matrix_Pvalues[a]);
			free(matrix_Pvaluesequal[a]);
			for(n=0;n<(long int)totalnloci+2;n++) 
				free(percentiles[a][n]);
			free(percentiles[a]);
		}
		free(validiter);
		free(average);
		free(variance);
		free(matrix_avg);
		free(matrix_var);
		free(matrix_Pvalues);
		free(matrix_Pvaluesequal);
		free(percentiles);
		if((*data)->likelihood_line == 0) fclose(filePvalue);
	}	
    return 0;
}
void calc_neutpar(int valuep,long int segsit,struct var2 **inputp, struct dnapar *ntpar,double valuer)
{
    long int pi;
    double k_;
    long int S;
    int nhapl,maxhapl;
    int B;
    int A;
    int *haplotype = 0;
    int *piw = 0;
    int *pib = 0;
    long int segsitesm1;
	int pidcount,npw;
    
    int inits,inits1,inits2,initcum,nsam;
    int val10,val20,val21;
    long int j,k;
    int a,b,c,d,h,i,comb,comb2;
    char *hapl;
    int ispolnomhit(long int,int,int);
	int **veca;
	/**/
	int mh1,mut,maxhapl1,*sshg;
    /**/
	int Min_rec(int,int,int,int);

    if(valuep == 0) {
        inits = 0;
        nsam  = (*inputp)->nsam;
    }
    else {
        if(valuep == 1) {
            inits = 0;
            nsam  = (*inputp)->config[0];
        }
        else {
            inits = (*inputp)->config[0];
            nsam  = (*inputp)->config[1];
        } 
    }
    comb = (int)((double)nsam*((double)nsam-(double)1)/(double)2);
    
   if(segsit == 0) {
        (ntpar)->B1 = 0;
        (ntpar)->Q1 = 0;
    	(ntpar)->k = (double)0;
    	(ntpar)->S = 0;
    	(ntpar)->nhapl = 1;
    	for(a=0;a<nsam;a++) {
            (ntpar)->freq[a] = 0;
            (ntpar)->unic[a] = 0;
            (ntpar)->fhapl[a] = 0;
        }
        (ntpar)->fhapl[0] = nsam;
        (ntpar)->piw = (double)0;
        (ntpar)->pib = (double)0;
        (ntpar)->maxhapl = nsam;
        (ntpar)->maxhapl1 = nsam;
		(ntpar)->Rm = (int)0;
    }
    else {
		if((veca = (int **)calloc(segsit,sizeof(int *))) == 0) {
			puts("calloc error veca.1");
			exit(1);
		}
		for(d=0;d<(int)segsit;d++) {
			if((veca[d] = (int *)calloc((*inputp)->nsam,sizeof(int))) == 0) {
				puts("calloc error veca.2");
				exit(1);
			}
		}
    
		if(valuer) (ntpar)->Rm = Min_rec(0,segsit,nsam,inits);
		else (ntpar)->Rm = (int)0;

        /* calcul B' i Q, pero sense dividir per S (Wall) */
        B = 0;
        A = 0;

        a = 0;
        if((segsitesm1 = segsit-1) < 0) segsitesm1 = 0;
        for(j=0;j<(long int)segsitesm1;) {
            k = j;
            while(k+1 < segsit) { /*calcular k*/
                if(ispolnomhit(k,inits,nsam) > 0) break;
                else k++;
            }
            j = k+1;
            while(j < segsit) { /*calcular j*/
                if(ispolnomhit(j,inits,nsam) > 0) break;
                else j++;
            }
            if(j < segsit) {                
                val20 = val21 = -1;
                b = 0;
                for(i=inits;i<inits+nsam;i++) {
                    val10 = (list[i][k] - 48)*4 + (list[i][j]- 48);
                    if(val20 == -1) val20 = val10;
                    else{
                        if(val20 != val10) {
                            if(val21 == -1) val21 = val10;
                            else  if(val21 != val10) {
                                b = 1;
                                break;
                            }
                        }
                    }
                }
				if(!b) {/*congruent*/
					B += 1;
					if(!A) {
						for(i=inits;i<inits+nsam;i++)
							veca[A][i] = list[i][j];
						A += 1;
					}
					else {
						a = 0;
						for(c=0;c<A;c++) {
							d = 0;
							for(i=inits;i<inits+nsam;i++) {
								if(veca[c][i] == list[i][j]) {
									d += 1;
								}
							}
							if(d == nsam || d == 0) {
								a = 1;
								break;
							}
						}
						if(!a) {
							for(i=inits;i<inits+nsam;i++)
								veca[A][i] = list[i][j];
							A += 1;
						}
					}
				}
			}
		}
		
        (ntpar)->Q1 = B+A;
        (ntpar)->B1 = B;
        
		for(d=0;d<(int)segsit;d++) free(veca[d]);
		free(veca);

        /* calcul de S, pi, fr, nhapl i unics/seq (excluding mhits)*/
        if((hapl = (char *)calloc(nsam*segsit,sizeof(char))) == NULL) {
            perror("calloc error calc_neutpar.0");
            exit(1);
        }
        k_ = (double)0;
        for(a=0;a<nsam;a++) {
            (ntpar)->freq[a] = 0;
            (ntpar)->unic[a] = 0;
        }
		pidcount = 0;
		for(a=inits;a<inits+nsam-1;a++) {
			for(b=a+1;b<inits+nsam;b++) {
				(ntpar)->pid[pidcount] = (double)0;
				pidcount++;
			}
		}
        S = 0;
        nhapl = 0;                
        for(j=0;j<segsit;j++) {
            pi = 0;
            while(j < segsit) {
                if((h=ispolnomhit(j,inits,nsam)) > 0) break; /*h is the frequency of the new mutation*/
                else j++;
            }            
            if(j<segsit) {
                (ntpar)->freq[h] += 1;
                for(a=inits;a<inits+nsam;a++) {
                    hapl[(a-inits)*segsit+S] = list[a][j];
                    /*new two lines: for singleton mutations (no outgroup)*/
                    if(h == 1) if(list[a][j] != '0') (ntpar)->unic[a-inits] += 1;
                    if(h == nsam-1) if(list[a][j] == '0') (ntpar)->unic[a-inits] += 1;
                }
                S++;
				pidcount = 0;
                for(a=inits;a<inits+nsam-1;a++) {
                    for(b=a+1;b<inits+nsam;b++) {
						if(list[a][j] != list[b][j]) {
							pi++;
							(ntpar)->pid[pidcount] += (double)1;
						}
						pidcount++;
					}
				}
				k_ += pi;
            }
        }
        (ntpar)->k = k_/(double)comb;
        (ntpar)->S = S;

        if((haplotype = (int *)malloc(nsam*sizeof(int))) == NULL) {
            perror("calloc error calc_neutpar.1");
            exit(1);
        }
        (ntpar)->nhapl = 0;
        for(a=0;a<nsam;a++) haplotype[a] = 1;            
        for(a=0;a<nsam-1;a++) {
            if(haplotype[a]) {
                (ntpar)->nhapl += 1;
                for(b=a+1;b<nsam;b++) {
                    if(haplotype[b]) {
                        if(memcmp(hapl+a*segsit,hapl+b*segsit,S) == 0) { 
                            haplotype[a] += 1;
                            haplotype[b] = 0;
                        }
                    }
                }
            }
        }
        if(haplotype[a]) (ntpar)->nhapl += 1;
        
        for(a=0;a<nsam;a++) (ntpar)->fhapl[a] = haplotype[a]; /*calcular freq. haplotips*/
        
        maxhapl=0;
        for(a=0;a<nsam;a++) 
			if(haplotype[a]>maxhapl)
				(ntpar)->maxhapl = maxhapl = haplotype[a]; /*calcular maxfreq. haplotips*/

		/**/
		if(S>0) {
			if((sshg = (int *)calloc(S,sizeof(int))) == NULL) {
				perror("calloc error calc_neutpar.0");
				exit(1);
			}
			maxhapl1 = 0;
			for(a=0;a<nsam;a++) {
				for(j=0;j<S;j++) {
					sshg[j] = 0;
				}
				for(b=0;b<nsam;b++) {
					mut = 0;
					for(j=0;j<S;j++) {
						if(hapl[a*segsit+j] != hapl[b*segsit+j])
							mut += 1;
					}
					if(mut <= 1) {
						for(j=0;j<S;j++) {
							if(hapl[a*segsit+j] != hapl[b*segsit+j]) {
								sshg[j] += 1;
							}
							if(mut == 0) sshg[j] += 1;
						}
					}
				}
				mh1 = 0;
				for(j=0;j<S;j++) {
					if(mh1 < sshg[j]) mh1 = sshg[j];
				}
				if(maxhapl1 < mh1)
					maxhapl1 = mh1;
			}
			(ntpar)->maxhapl1 = maxhapl1;
			free(sshg);
		}
		else {
			(ntpar)->maxhapl1 = nsam;
		}
		/**/
/*
        b = 0;
        for(a=0;a<nsam;a++) if(haplotype[b] < haplotype[a]) b=a;
        printf("%d\r",haplotype[b]);
*/
        free(hapl);
        free(haplotype);

        /*calcular thetaL*/
		(ntpar)->thetaL = 0;
		for(a=0;a<nsam;a++) {
            (ntpar)->thetaL += (double)a * (double)ntpar->freq[a];
		}
		(ntpar)->thetaL *= (double)1/(double)(nsam-1);
        /* Calcular pi_within i pi_between: nomes es permet per poblacions amb config > 1*/
        if((*inputp)->npop_sampled > 1) {
            if((piw = (int *)calloc((*inputp)->npop_sampled,sizeof(int))) == NULL) {
                perror("calloc error calc_neutpar.0");
                exit(1);
            }
            if((pib =(int *)calloc(((*inputp)->npop_sampled*((*inputp)->npop_sampled-1)/2),sizeof(int))) ==NULL) {
                perror("calloc error calc_neutpar.0");
                exit(1);
            }        
            (ntpar)->piw = 0.;
            (ntpar)->withinw = 0.;
            (ntpar)->pib = 0.;
            for(j=0;j<segsit;j++) {
                while(j < segsit) {
                    if((c=ispolnomhit(j,inits,nsam)) > 0) break;
                    else j++;
                }                    
                if(j<segsit) {
                    inits1 = 0;
                    for(h=0;h<(*inputp)->npop_sampled;h++) {
                        for(a=inits1;a<inits1+(*inputp)->config[h]-1;a++)
                            for(b=a+1;b<inits1+(*inputp)->config[h];b++)
                                if(list[a][j] != list[b][j]) piw[h]++;
                        inits1 += (*inputp)->config[h];
                    }
                    comb2 = 0;
                    inits1 = initcum = inits2 = 0;
                    for(h=0;h<(*inputp)->npop_sampled-1;h++) {
                        initcum = inits2 = initcum + (*inputp)->config[h];
                        for(i=h+1;i<(*inputp)->npop_sampled;i++) {
                            if((*inputp)->config[h] > 1 && (*inputp)->config[i] > 1) {
                                for(a=inits1;a<inits1+(*inputp)->config[h];a++)
                                    for(b=inits2;b<inits2+(*inputp)->config[i];b++)
                                        if(list[a][j] != list[b][j]) pib[comb2]++;
                                comb2++;
                            }
                            inits2 += (*inputp)->config[i];
                        }
                        inits1 += (*inputp)->config[h];
                    }
                }
            }
            comb2 = 0;
            for(h=0;h<(*inputp)->npop_sampled;h++) {
                if((*inputp)->config[h] > 1) {
					comb = (int)((double)(*inputp)->config[h] * ((double)(*inputp)->config[h] - (double)1.0)/(double)2.0);
                    (ntpar)->piw += (double)piw[h]/(double)comb; /*es recull el valor a cada subpop indiviadualment*/
					(ntpar)->withinw += (double)(*inputp)->config[h] * (double)piw[h]/(double)comb;
                    comb2++;
                }
            }
            if(comb2) {
				(ntpar)->piw = (ntpar)->piw/(double)comb2;/*despres es divideix pel nombre de subpoblacions*/
				npw=0;
				for(h=0;h<(*inputp)->npop_sampled;h++) npw += (*inputp)->config[h];
				(ntpar)->withinw /= (double)npw;
			}
            comb2 = 0;
            for(h=0;h<(*inputp)->npop_sampled-1;h++) {
                for(i=h+1;i<(*inputp)->npop_sampled;i++) {
                    if((*inputp)->config[h] > 1 && (*inputp)->config[i] > 1) {
                        comb = (*inputp)->config[h] * (*inputp)->config[i];
                        (ntpar)->pib += (double)pib[comb2]/(double)comb;
                        comb2++;
                    }
                }
            }
            if(comb2) (ntpar)->pib = (ntpar)->pib/(double)(comb2);
            free(piw);
            free(pib);
        }
    }
}


double Zns(int valuep,long int segsites,struct var2 **inputp)
{
    double ZnS;
    long k,j,comb;
    int i,inits,nsam;
    int vala,valb,val00,a,b,na,nb;
    double A,B,C;
    int ispolnomhit(long int,int,int);
    
    if(valuep == 0) {
        inits = 0;
        nsam  = (*inputp)->nsam;
    }
    else {
        if(valuep == 1) {
            inits = 0;
            nsam  = (*inputp)->config[0];
        }
        else {
            inits = (*inputp)->config[0];
            nsam  = (*inputp)->config[1];
        } 
    }
    
    ZnS = (double)0;
    j = 0;
    comb = 0;
    while(j+1 < (long)segsites) {
        while(j < (long)segsites) {
            if(ispolnomhit(j,inits,nsam) > 0) break;
            else j++;
        }
        k = j+1;
        while(k < (long)segsites) {
            while(k < (long)segsites) {
                if(ispolnomhit(k,inits,nsam) > 0) break;
                else k++;
            }
            if(k < (long)segsites) {
                /*calcular freqs p1,q1*/
                vala = valb = 1;
                a = list[inits][j];
                b = list[inits][k];
                for(i=1+inits;i<inits+nsam;i++) {
                    if(list[i][j] == a) vala++;
                    else na = list[i][j];
                    if(list[i][k] == b) valb++;                                
                    else nb = list[i][k];
                }
                if(nsam - vala > vala) {
                    a = na;
                    vala = nsam - vala;
                }
                if(nsam - valb > valb) {
                    b = nb;
                    valb = nsam - valb;
                }
                /*calcular p1q1*/
                val00 = 0;
                for(i=0+inits;i<inits+nsam;i++) if(list[i][j] == a && list[i][k] == b) val00++;
                
                /*calcular R2 = (p1q1 - (p1*q1))"2 / (p1*(1-p1) * (q1*(1-q1))) */
                A = (double)vala/(double)nsam;
                B = (double)valb/(double)nsam;
                C = (double)val00/(double)nsam;
                /*if(vala != ((double)nsam-1.0) && valb != ((double)nsam-1.0)) if only informative */
                ZnS += ((C - (A*B)) * (C - (A*B))) / (A*(1.0 - A)*B*(1.0 - B));
                comb++;
            }
            k++;
        }
        j++;
    }
    if(comb) ZnS = ZnS/(double)comb;
    else return(-10000);
    return(ZnS);
}

double ZnA_(int valuep,long int segsites,struct var2 **inputp)
{
    double ZnA;
    long k,j,comb;
    int i,inits,nsam;
    int vala,valb,val00,a,b,na,nb;
    double A,B,C;
    int ispolnomhit(long int,int,int);
    
    if(valuep == 0) {
        inits = 0;
        nsam  = (*inputp)->nsam;
    }
    else {
        if(valuep == 1) {
            inits = 0;
            nsam  = (*inputp)->config[0];
        }
        else {
            inits = (*inputp)->config[0];
            nsam  = (*inputp)->config[1];
        } 
    }
    
    ZnA = (double)0;
    j = 0;
    comb = 0;
    while(j+1 < (long)segsites) {
        while(j < (long)segsites) {
            if(ispolnomhit(j,inits,nsam) > 0) break;
            else j++;
        }
        k = j+1;
        while(k < (long)segsites) {
            if(ispolnomhit(k,inits,nsam) > 0) break;
            else k++;
        }
        if(k < (long)segsites) {
            /*calcular freqs p1,q1*/
            vala = valb = 1;
            a = list[inits][j];
            b = list[inits][k];
            for(i=1+inits;i<inits+nsam;i++) {
                if(list[i][j] == a) vala++;
                else na = list[i][j];
                if(list[i][k] == b) valb++;                                
                else nb = list[i][k];
            }
            if(nsam - vala > vala) {
                a = na;
                vala = nsam - vala;
            }
            if(nsam - valb > valb) {
                b = nb;
                valb = nsam - valb;
            }
            /*calcular p1q1*/
            val00 = 0;
            for(i=0+inits;i<inits+nsam;i++) if(list[i][j] == a && list[i][k] == b) val00++;
            
            /*calcular R2 = (p1q1 - (p1*q1))"2 / (p1*(1-p1) * (q1*(1-q1))) */
            A = (double)vala/(double)nsam;
            B = (double)valb/(double)nsam;
            C = (double)val00/(double)nsam;
            /*if(vala != ((double)nsam-1.0) && valb != ((double)nsam-1.0)) if only informative */
            ZnA += ((C - (A*B)) * (C - (A*B))) / (A*(1.0 - A)*B*(1.0 - B));
            comb++;
        }
        j++;
    }
    if(comb) ZnA = ZnA/(double)comb;
    else return(-10000);
    return(ZnA);
}

int ispolnomhit(long int j,int init, int nsam)
{
    int i,h,g;
    int a1;
    
    if(j) if(posit[j-1] == posit[j]) return(-1);/*estem a la segona mutació o més*/
    
    h = g = a1 = 0;
    for(i=init;i<init+nsam;i++) {
        if(list[i][j] == '0') h++;
        else {
            if(!a1) a1 = list[i][j];
            if(list[i][j] == a1) g++;
            else return(-2);/*en el cas de la primera mutacio i hagin almenys tres variants*/
        }
    }
    if(h==nsam || h == 0) return(0);
    return(nsam-h);
}

double koutg(int nsam, int Sout, int *freq, long int nsites)
{
    int x;
    double div;
    
    if(nsam) {
        div = (double)Sout;
        for(x=0;x<nsam;x++) div += (double)freq[x]/(double)nsam;
        div/= (double)nsites;
        return -(double)0.75*(double)log((double)((double)1.0-(double)4.0/(double)3.0*div));/*Jukes and Cantor correction...*/
    }
    else return -10000;
}


void calc_neutpar_window(struct var2 **inputp,struct dnapar *ntpar,long int  s0,long int  s1,double valuer)
{
    long int segsit;    
    long int pi;
    double k_;
    long int S;
    int nhapl,maxhapl;
    int B;
    int A;
    int *haplotype = 0;
    int *piw = 0;
    int *pib = 0;
    
    int inits1,inits2,initcum,nsam;
    int val10,val20,val21;
    long int j,k;
    int a,b,c,d,h,i,comb,comb2;
    char *hapl;
    int ispolnomhit(long int,int,int);
	int **veca;
	/**/
	int mut,mh1,maxhapl1,*sshg;
    /**/
	int Min_rec(int,int,int,int);
	
	int pidcount,npw;

    nsam  = (*inputp)->nsam;
    comb = (int)((double)nsam*((double)nsam-(double)1)/(double)2);
    
    /*define segsit first*/
    segsit = s1 - s0;
    
    if(segsit <= 0) {
        (ntpar)->B1 = 0;
        (ntpar)->Q1 = 0;
    	(ntpar)->k = (double)0;
    	(ntpar)->S = 0;
    	(ntpar)->nhapl = 1;
    	for(a=0;a<nsam;a++) {
            (ntpar)->freq[a] = 0;
            (ntpar)->unic[a] = 0;
            (ntpar)->fhapl[a] = 0;
        }
        (ntpar)->fhapl[0] = nsam;
        (ntpar)->piw = (double)0;
        (ntpar)->pib = (double)0;
        (ntpar)->maxhapl = nsam;
        (ntpar)->maxhapl1 = nsam;
		(ntpar)->Rm = (int)0;
    }
    else {
		if((veca = (int **)calloc(segsit,sizeof(int *))) == 0) {
			puts("calloc error veca.1");
			exit(1);
		}
		for(d=0;d<(int)segsit;d++) {
			if((veca[d] = (int *)calloc((*inputp)->nsam,sizeof(int))) == 0) {
				puts("calloc error veca.2");
				exit(1);
			}
		}
		
		if(valuer) (ntpar)->Rm = Min_rec(s0,s1,nsam,0);
		else (ntpar)->Rm = (int)0;

        /* calcul B' i Q, pero sense dividir per S (Wall) */
        B = 0;
        A = 0;

        a = 0;
        for(j=s0;j<s1;) {
            k = j;
            while(k+1 < s1) { /*calcular k*/
                if(ispolnomhit(k,0,nsam) > 0) break;
                else k++;
            }
            j = k+1;
            while(j < s1) { /*calcular j*/
                if(ispolnomhit(j,0,nsam) > 0) break;
                else j++;
            }
            if(j < s1) {                
                val20 = val21 = -1;
                b = 0;
                for(i=0;i<0+nsam;i++) {
                    val10 = (list[i][k] - 48)*4 + (list[i][j]- 48);
                    if(val20 == -1) val20 = val10;
                    else{
                        if(val20 != val10) {
                            if(val21 == -1) val21 = val10;
                            else  if(val21 != val10) {
                                b = 1;
                                break;
                            }
                        }
                    }
                }
				if(!b) {/*congruent*/
					B += 1;
					if(!A) {
						for(i=0;i<0+nsam;i++)
							veca[A][i] = list[i][j];
						A += 1;
					}
					else {
						a = 0;
						for(c=0;c<A;c++) {
							d = 0;
							for(i=0;i<0+nsam;i++) {
								if(veca[c][i] == list[i][j]) {
									d += 1;
								}
							}
							if(d == nsam || d == 0) {
								a = 1;
								break;
							}
						}
						if(!a) {
							for(i=0;i<0+nsam;i++)
								veca[A][i] = list[i][j];
							A += 1;
						}
					}
				}
            }
        }
		
        (ntpar)->Q1 = B+A;
        (ntpar)->B1 = B;
        
		for(d=0;d<(int)segsit;d++) free(veca[d]);
		free(veca);
        
        /* calcul de S, pi, fr, nhapl i unics/seq (excluding mhits)*/
        if((hapl = (char *)calloc(nsam*segsit,sizeof(char))) == NULL) {
            perror("calloc error calc_neutpar.0");
            exit(1);
        }
        k_ = (double)0;
        for(a=0;a<nsam;a++) {
            (ntpar)->freq[a] = 0;
            (ntpar)->unic[a] = 0;
        }
		pidcount = 0;
		for(a=0;a<nsam-1;a++) {
			for(b=a+1;b<nsam;b++) {
				(ntpar)->pid[pidcount] = (double)0;
				pidcount++;
			}
		}
        S = 0;
        nhapl = 0;                
        for(j=s0;j<s1;j++) {
            pi = 0;
            while(j < s1) {
                if((h=ispolnomhit(j,0,nsam)) > 0) break; /*h is the frequency of the new mutation*/
                else j++;
            }            
            if(j<s1) {
                (ntpar)->freq[h] += 1;
                for(a=0;a<0+nsam;a++) {
                    hapl[(a-0)*segsit+S] = list[a][j];
                    /*new two lines: for singleton mutations (no outgroup)*/
                    if(h == 1) if(list[a][j] != '0') (ntpar)->unic[a-0] += 1;
                    if(h == nsam-1) if(list[a][j] == '0') (ntpar)->unic[a-0] += 1;
                }
                S++;
				pidcount = 0;
                for(a=0;a<0+nsam-1;a++) {
                    for(b=a+1;b<0+nsam;b++) {
						if(list[a][j] != list[b][j]) {
							pi++;
							(ntpar)->pid[pidcount] += (double)1;
						}
						pidcount++;
					}
				}
                k_ += pi;
            }
        }
        (ntpar)->k = k_/(double)comb;
        (ntpar)->S = S;

        if((haplotype = (int *)malloc(nsam*sizeof(int))) == NULL) {
            perror("calloc error calc_neutpar.1");
            exit(1);
        }
        (ntpar)->nhapl = 0;
        for(a=0;a<nsam;a++) haplotype[a] = 1;            
        for(a=0;a<nsam-1;a++) {
            if(haplotype[a]) {
                (ntpar)->nhapl += 1;
                for(b=a+1;b<nsam;b++) {
                    if(haplotype[b]) {
                        if(memcmp(hapl+a*segsit,hapl+b*segsit,S) == 0) { 
                            haplotype[a] += 1;
                            haplotype[b] = 0;
                        }
                    }
                }
            }
        }
        if(haplotype[a]) (ntpar)->nhapl += 1;
        
        for(a=0;a<nsam;a++) (ntpar)->fhapl[a] = haplotype[a]; /*calcular freq. haplotips*/
        
        maxhapl=0;
        for(a=0;a<nsam;a++) 
			if(haplotype[a]>maxhapl)
				(ntpar)->maxhapl = maxhapl = haplotype[a]; /*calcular maxfreq. haplotips*/
		/**/
		if(S>0) {
			if((sshg = (int *)calloc(S,sizeof(int))) == NULL) {
				perror("calloc error calc_neutpar.0");
				exit(1);
			}
			maxhapl1 = 0;
			for(a=0;a<nsam;a++) {
				for(j=0;j<S;j++) {
					sshg[j] = 0;
				}
				for(b=0;b<nsam;b++) {
					mut = 0;
					for(j=0;j<S;j++) {
						if(hapl[a*segsit+j] != hapl[b*segsit+j])
							mut += 1;
					}
					if(mut <= 1) {
						for(j=0;j<S;j++) {
							if(hapl[a*segsit+j] != hapl[b*segsit+j]) {
								sshg[j] += 1;
							}
							if(mut == 0) sshg[j] += 1;
						}
					}
				}
				mh1 = 0;
				for(j=0;j<S;j++) {
					if(mh1 < sshg[j]) mh1 = sshg[j];
				}
				if(maxhapl1 < mh1)
					maxhapl1 = mh1;
			}
			(ntpar)->maxhapl1 = maxhapl1;
			free(sshg);
		}
		else {
			(ntpar)->maxhapl1 = nsam;
		}
		/**/
/*
        b = 0;
        for(a=0;a<nsam;a++) if(haplotype[b] < haplotype[a]) b=a;
        printf("%d\r",haplotype[b]);
*/
        free(hapl);
        free(haplotype);

        /*calcular thetaL*/
		(ntpar)->thetaL = 0;
		for(a=0;a<nsam;a++) {
            (ntpar)->thetaL += (double)a * (double)ntpar->freq[a];
		}
		(ntpar)->thetaL *= (double)1/(double)(nsam-1);

        /* Calcular pi_within i pi_between: nomes es permet per poblacions amb config > 1 */
        if((*inputp)->npop_sampled > 1) {
            if((piw = (int *)calloc((*inputp)->npop_sampled,sizeof(int))) == NULL) {
                perror("calloc error calc_neutpar.0");
                exit(1);
            }
            if((pib =(int *)calloc(((*inputp)->npop_sampled*((*inputp)->npop_sampled-1)/2),sizeof(int))) ==NULL) {
                perror("calloc error calc_neutpar.0");
                exit(1);
            }        
            (ntpar)->piw = (double)0;
            (ntpar)->pib = (double)0;
            for(j=s0;j<s1;j++) {
                while(j < s1) {
                    if((c=ispolnomhit(j,0,nsam)) > 0) break;
                    else j++;
                }                    
                if(j<s1) {
                    inits1 = 0;
                    for(h=0;h<(*inputp)->npop_sampled;h++) {
                        for(a=inits1;a<inits1+(*inputp)->config[h]-1;a++)
                            for(b=a+1;b<inits1+(*inputp)->config[h];b++)
                                if(list[a][j] != list[b][j]) piw[h]++;
                        inits1 += (*inputp)->config[h];
                    }
                    comb2 = 0;
                    inits1 = initcum = inits2 = 0;
                    for(h=0;h<(*inputp)->npop_sampled-1;h++) {
                        initcum = inits2 = initcum + (*inputp)->config[h];
                        for(i=h+1;i<(*inputp)->npop_sampled;i++) {
                            if((*inputp)->config[h] > 1 && (*inputp)->config[i] > 1) {
                                for(a=inits1;a<inits1+(*inputp)->config[h];a++)
                                    for(b=inits2;b<inits2+(*inputp)->config[i];b++)
                                        if(list[a][j] != list[b][j]) pib[comb2]++;
                                comb2++;
                            }
                            inits2 += (*inputp)->config[i];
                        }
                        inits1 += (*inputp)->config[h];
                    }
                }
            }
            comb2 = 0;
            for(h=0;h<(*inputp)->npop_sampled;h++) {
                if((*inputp)->config[h] > 1) {
					comb = (int)((double)(*inputp)->config[h] * ((double)(*inputp)->config[h] - (double)1)/(double)2);
                    (ntpar)->piw += (double)piw[h]/(double)comb; /*es recull el valor a cada subpop indiviadualment*/
					(ntpar)->withinw += (double)(*inputp)->config[h] * (double)piw[h]/(double)comb;
                    comb2++;
                }
            }
            if(comb2) {
				(ntpar)->piw = (ntpar)->piw/(double)comb2;/*despres es divideix pel nombre de subpoblacions*/
				npw=0;
				for(h=0;h<(*inputp)->npop_sampled;h++) npw += (*inputp)->config[h];
				(ntpar)->withinw /= (double)npw;
			}
            
            comb2 = 0;
            for(h=0;h<(*inputp)->npop_sampled-1;h++) {
                for(i=h+1;i<(*inputp)->npop_sampled;i++) {
                    if((*inputp)->config[h] > 1 && (*inputp)->config[i] > 1) {
                        comb = (*inputp)->config[h] * (*inputp)->config[i];
                        (ntpar)->pib += (double)pib[comb2]/(double)comb;
                        comb2++;
                    }
                }
            }
            if(comb2) (ntpar)->pib = (ntpar)->pib/(double)(comb2);
            free(piw);
            free(pib);
        }
    }
}

double Zns_window(struct var2 **inputp, long int s0, long int s1)
{
    double ZnS;
    long k,j,comb;
    int i,nsam;
    int vala,valb,val00,a,b,na,nb;
    double A,B,C;
    int ispolnomhit(long int,int,int);
    
    nsam  = (*inputp)->nsam;
    
    ZnS = (double)0;
    j = s0;
    comb = 0;
    while(j+1 < (long)s1) {
        while(j < (long)s1) {
            if(ispolnomhit(j,0,nsam) > 0) break;
            else j++;
        }
        k = j+1;
        while(k < (long)s1) {
            while(k < (long)s1) {
                if(ispolnomhit(k,0,nsam) > 0) break;
                else k++;
            }
            if(k < (long)s1) {
                /*calcular freqs p1,q1*/
                vala = valb = 1;
                a = list[0][j];
                b = list[0][k];
                for(i=1+0;i<0+nsam;i++) {
                    if(list[i][j] == a) vala++;
                    else na = list[i][j];
                    if(list[i][k] == b) valb++;                                
                    else nb = list[i][k];
                }
                if(nsam - vala > vala) {
                    a = na;
                    vala = nsam - vala;
                }
                if(nsam - valb > valb) {
                    b = nb;
                    valb = nsam - valb;
                }
                /*calcular p1q1*/
                val00 = 0;
                for(i=0+0;i<0+nsam;i++) if(list[i][j] == a && list[i][k] == b) val00++;
                
                /*calcular R2 = (p1q1 - (p1*q1))"2 / (p1*(1-p1) * (q1*(1-q1))) */
                A = (double)vala/(double)nsam;
                B = (double)valb/(double)nsam;
                C = (double)val00/(double)nsam;
                /*if(vala != ((double)nsam-1.0) && valb != ((double)nsam-1.0)) if only informative */
                ZnS += ((C - (A*B)) * (C - (A*B))) / (A*(1.0 - A)*B*(1.0 - B));
                comb++;
            }
            k++;
        }
        j++;
    }
    if(comb) ZnS = ZnS/(double)comb;
    else return(-10000);
    return(ZnS);
}

double ZnA_window(struct var2 **inputp, long int s0, long int s1)
{
    double ZnA;
    long k,j,comb;
    int i,nsam;
    int vala,valb,val00,a,b,na,nb;
    double A,B,C;
    int ispolnomhit(long int,int,int);
    
    nsam  = (*inputp)->nsam;
    
    ZnA = (double)0;
    j = s0;
    comb = 0;
    while(j+1 < (long)s1) {
        while(j < (long)s1) {
            if(ispolnomhit(j,0,nsam) > 0) break;
            else j++;
        }
        k = j+1;
        while(k < (long)s1) {
            if(ispolnomhit(k,0,nsam) > 0) break;
            else k++;
        }
        if(k < (long)s1) {
            /*calcular freqs p1,q1*/
            vala = valb = 1;
            a = list[0][j];
            b = list[0][k];
            for(i=1+0;i<0+nsam;i++) {
                if(list[i][j] == a) vala++;
                else na = list[i][j];
                if(list[i][k] == b) valb++;                                
                else nb = list[i][k];
            }
            if(nsam - vala > vala) {
                a = na;
                vala = nsam - vala;
            }
            if(nsam - valb > valb) {
                b = nb;
                valb = nsam - valb;
            }
            /*calcular p1q1*/
            val00 = 0;
            for(i=0+0;i<0+nsam;i++) if(list[i][j] == a && list[i][k] == b) val00++;
            
            /*calcular R2 = (p1q1 - (p1*q1))"2 / (p1*(1-p1) * (q1*(1-q1))) */
            A = (double)vala/(double)nsam;
            B = (double)valb/(double)nsam;
            C = (double)val00/(double)nsam;
            /*if(vala != ((double)nsam-1.0) && valb != ((double)nsam-1.0)) if only informative */
            ZnA += ((C - (A*B)) * (C - (A*B))) / (A*(1.0 - A)*B*(1.0 - B));
            comb++;
        }
        j++;
    }
    if(comb) ZnA = ZnA/(double)comb;
    else return(-10000);
    return(ZnA);
}

int print_matrix_sfixalltheta(struct var **data,FILE *file_thetap,FILE *file_ttotp,FILE *file_Srange) 
{
	int x;
	long int j;
	
	if(file_Srange) x=0;
	else x=0;
	
    if((*data)->ifgamma == 1 || (*data)->range_thetant) {
		for(x=0;x<(*data)->n_loci;x++) {                                
			#if ADDITIONAL_PRINTS == 1 || ADDITIONAL_PRINTS == 12 || ADDITIONAL_PRINTS == 13 || ADDITIONAL_PRINTS == 123
			fprintf(file_thetap,"thetap[%d]\t",x);
			#endif
			#if ADDITIONAL_PRINTS == 2 || ADDITIONAL_PRINTS == 12 || ADDITIONAL_PRINTS == 23 || ADDITIONAL_PRINTS == 123
			fprintf(file_ttotp,"Ttot[%d]\t",x);
			#endif
			#if ADDITIONAL_PRINTS == 3 || ADDITIONAL_PRINTS == 13 || ADDITIONAL_PRINTS == 23 || ADDITIONAL_PRINTS == 123
			fprintf(file_Srange,"thetae[%d]\t",x);
			#endif
		}
		#if ADDITIONAL_PRINTS == 1 || ADDITIONAL_PRINTS == 12 || ADDITIONAL_PRINTS == 13 || ADDITIONAL_PRINTS == 123
		fprintf(file_thetap,"\n");
		#endif
		#if ADDITIONAL_PRINTS == 2 || ADDITIONAL_PRINTS == 12 || ADDITIONAL_PRINTS == 23 || ADDITIONAL_PRINTS == 123
		fprintf(file_ttotp,"\n");
		#endif
		#if ADDITIONAL_PRINTS == 3 || ADDITIONAL_PRINTS == 13 || ADDITIONAL_PRINTS == 23 || ADDITIONAL_PRINTS == 123
		fprintf(file_Srange,"\n");
		#endif
		for(j=0;j<(*data)->n_iter;j++) {                                
			for(x=0;x<(*data)->n_loci;x++) {                                
				#if ADDITIONAL_PRINTS == 1 || ADDITIONAL_PRINTS == 12 || ADDITIONAL_PRINTS == 13 || ADDITIONAL_PRINTS == 123
				fprintf(file_thetap,"%G\t",postp[x][j].thetap);
				#endif
				#if ADDITIONAL_PRINTS == 2 || ADDITIONAL_PRINTS == 12 || ADDITIONAL_PRINTS == 23 || ADDITIONAL_PRINTS == 123
				fprintf(file_ttotp,"%G\t",postp[x][j].Ttotp);
				#endif
				#if ADDITIONAL_PRINTS == 3 || ADDITIONAL_PRINTS == 13 || ADDITIONAL_PRINTS == 23 || ADDITIONAL_PRINTS == 123
				fprintf(file_Srange,"%G\t",thetafromS[x][j]);
				#endif
			}
			#if ADDITIONAL_PRINTS == 1 || ADDITIONAL_PRINTS == 12 || ADDITIONAL_PRINTS == 13 || ADDITIONAL_PRINTS == 123
			fprintf(file_thetap,"\n");
			#endif
			#if ADDITIONAL_PRINTS == 2 || ADDITIONAL_PRINTS == 12 || ADDITIONAL_PRINTS == 23 || ADDITIONAL_PRINTS == 123
			fprintf(file_ttotp,"\n");
			#endif
			#if ADDITIONAL_PRINTS == 3 || ADDITIONAL_PRINTS == 13 || ADDITIONAL_PRINTS == 23 || ADDITIONAL_PRINTS == 123
			fprintf(file_Srange,"\n");
			#endif
		}
		#if ADDITIONAL_PRINTS == 1 || ADDITIONAL_PRINTS == 12 || ADDITIONAL_PRINTS == 13 || ADDITIONAL_PRINTS == 123
		fprintf(file_thetap,"\n");
		#endif
		#if ADDITIONAL_PRINTS == 2 || ADDITIONAL_PRINTS == 12 || ADDITIONAL_PRINTS == 23 || ADDITIONAL_PRINTS == 123
		fprintf(file_ttotp,"\n");
		#endif
		#if ADDITIONAL_PRINTS == 3 || ADDITIONAL_PRINTS == 13 || ADDITIONAL_PRINTS == 23 || ADDITIONAL_PRINTS == 123
		fprintf(file_Srange,"\n");
		#endif
		if((*data)->sfix_allthetas > 0) {
			for(x=0;x<(*data)->n_loci;x++) {                                
				#if ADDITIONAL_PRINTS == 1 || ADDITIONAL_PRINTS == 12 || ADDITIONAL_PRINTS == 13 || ADDITIONAL_PRINTS == 123
				fprintf(file_thetap,"ratio[%d]\t",x);
				#endif
				#if ADDITIONAL_PRINTS == 2 || ADDITIONAL_PRINTS == 12 || ADDITIONAL_PRINTS == 23 || ADDITIONAL_PRINTS == 123
				fprintf(file_ttotp,"ratio[%d]\t",x);
				#endif
			}
			#if ADDITIONAL_PRINTS == 1 || ADDITIONAL_PRINTS == 12 || ADDITIONAL_PRINTS == 13 || ADDITIONAL_PRINTS == 123
			fprintf(file_thetap,"\n");
			#endif
			#if ADDITIONAL_PRINTS == 2 || ADDITIONAL_PRINTS == 12 || ADDITIONAL_PRINTS == 23 || ADDITIONAL_PRINTS == 123
			fprintf(file_ttotp,"\n");
			#endif
			for(x=0;x<(*data)->n_loci;x++) {                                
				#if ADDITIONAL_PRINTS == 1 || ADDITIONAL_PRINTS == 12 || ADDITIONAL_PRINTS == 13 || ADDITIONAL_PRINTS == 123
				fprintf(file_thetap,"%G\t",postp[x][(*data)->n_iter].thetap);
				#endif
				#if ADDITIONAL_PRINTS == 2 || ADDITIONAL_PRINTS == 12 || ADDITIONAL_PRINTS == 23 || ADDITIONAL_PRINTS == 123
				fprintf(file_ttotp,"%G\t",postp[x][(*data)->n_iter].Ttotp);
				#endif
			}
			#if ADDITIONAL_PRINTS == 1 || ADDITIONAL_PRINTS == 12 || ADDITIONAL_PRINTS == 13 || ADDITIONAL_PRINTS == 123
			fprintf(file_thetap,"\n");
			#endif
			#if ADDITIONAL_PRINTS == 2 || ADDITIONAL_PRINTS == 12 || ADDITIONAL_PRINTS == 23 || ADDITIONAL_PRINTS == 123
			fprintf(file_ttotp,"\n");
			#endif
		}
	}
	return 0;
}
int print_matrix_rmfix(struct var **data,FILE *file_recp,FILE *file_ttotp/*,FILE *file_recrange */,int Sfix) 
{
	int x;
	long int j;
	
    if((*data)->ifgammar == 1 || (*data)->range_rnt) {
		for(x=0;x<(*data)->n_loci;x++) {                                
			#if ADDITIONAL_PRINTS == 1 || ADDITIONAL_PRINTS == 12 || ADDITIONAL_PRINTS == 13 || ADDITIONAL_PRINTS == 123
			fprintf(file_recp,"recp[%d]\t",x);
			#endif
			#if ADDITIONAL_PRINTS == 2 || ADDITIONAL_PRINTS == 12 || ADDITIONAL_PRINTS == 23 || ADDITIONAL_PRINTS == 123
			if(Sfix == 0) fprintf(file_ttotp,"Ttot[%d]\t",x);
			#endif
			/*
			#if ADDITIONAL_PRINTS == 3 || ADDITIONAL_PRINTS == 13 || ADDITIONAL_PRINTS == 23 || ADDITIONAL_PRINTS == 123
			fprintf(file_recrange,"rece[%d]\t",x);
			#endif
			*/
		}
		#if ADDITIONAL_PRINTS == 1 || ADDITIONAL_PRINTS == 12 || ADDITIONAL_PRINTS == 13 || ADDITIONAL_PRINTS == 123
		fprintf(file_recp,"\n");
		#endif
		#if ADDITIONAL_PRINTS == 2 || ADDITIONAL_PRINTS == 12 || ADDITIONAL_PRINTS == 23 || ADDITIONAL_PRINTS == 123
		if(Sfix == 0) fprintf(file_ttotp,"\n");
		#endif
		/*
		#if ADDITIONAL_PRINTS == 3 || ADDITIONAL_PRINTS == 13 || ADDITIONAL_PRINTS == 23 || ADDITIONAL_PRINTS == 123
		fprintf(file_recrange,"\n");
		#endif
		*/
		for(j=0;j<(*data)->n_iter;j++) {                                
			for(x=0;x<(*data)->n_loci;x++) {                                
				#if ADDITIONAL_PRINTS == 1 || ADDITIONAL_PRINTS == 12 || ADDITIONAL_PRINTS == 13 || ADDITIONAL_PRINTS == 123
				fprintf(file_recp,"%G\t",postp[x][j].recp);
				#endif
				#if ADDITIONAL_PRINTS == 2 || ADDITIONAL_PRINTS == 12 || ADDITIONAL_PRINTS == 23 || ADDITIONAL_PRINTS == 123
				if(Sfix == 0) fprintf(file_ttotp,"%G\t",postp[x][j].Ttotp);
				#endif
				/*
				#if ADDITIONAL_PRINTS == 3 || ADDITIONAL_PRINTS == 13 || ADDITIONAL_PRINTS == 23 || ADDITIONAL_PRINTS == 123
				fprintf(file_recrange,"%G\t",thetafromS[x][j]);
				#endif
				*/
			}
			#if ADDITIONAL_PRINTS == 1 || ADDITIONAL_PRINTS == 12 || ADDITIONAL_PRINTS == 13 || ADDITIONAL_PRINTS == 123
			fprintf(file_recp,"\n");
			#endif
			#if ADDITIONAL_PRINTS == 2 || ADDITIONAL_PRINTS == 12 || ADDITIONAL_PRINTS == 23 || ADDITIONAL_PRINTS == 123
			if(Sfix == 0) fprintf(file_ttotp,"\n");
			#endif
			/*
			#if ADDITIONAL_PRINTS == 3 || ADDITIONAL_PRINTS == 13 || ADDITIONAL_PRINTS == 23 || ADDITIONAL_PRINTS == 123
			fprintf(file_recrange,"\n");
			#endif
			*/
		}
		#if ADDITIONAL_PRINTS == 1 || ADDITIONAL_PRINTS == 12 || ADDITIONAL_PRINTS == 13 || ADDITIONAL_PRINTS == 123
		fprintf(file_recp,"\n");
		#endif
		#if ADDITIONAL_PRINTS == 2 || ADDITIONAL_PRINTS == 12 || ADDITIONAL_PRINTS == 23 || ADDITIONAL_PRINTS == 123
		if(Sfix == 0) fprintf(file_ttotp,"\n");
		#endif
		/*
		#if ADDITIONAL_PRINTS == 3 || ADDITIONAL_PRINTS == 13 || ADDITIONAL_PRINTS == 23 || ADDITIONAL_PRINTS == 123
		fprintf(file_recrange,"\n");
		#endif
		*/
		if((*data)->rmfix > 0) {
			for(x=0;x<(*data)->n_loci;x++) {                                
				#if ADDITIONAL_PRINTS == 1 || ADDITIONAL_PRINTS == 12 || ADDITIONAL_PRINTS == 13 || ADDITIONAL_PRINTS == 123
				fprintf(file_recp,"ratio[%d]\t",x);
				#endif
				#if ADDITIONAL_PRINTS == 2 || ADDITIONAL_PRINTS == 12 || ADDITIONAL_PRINTS == 23 || ADDITIONAL_PRINTS == 123
				if(Sfix == 0) fprintf(file_ttotp,"ratio[%d]\t",x);
				#endif
			}
			#if ADDITIONAL_PRINTS == 1 || ADDITIONAL_PRINTS == 12 || ADDITIONAL_PRINTS == 13 || ADDITIONAL_PRINTS == 123
			fprintf(file_recp,"\n");
			#endif
			#if ADDITIONAL_PRINTS == 2 || ADDITIONAL_PRINTS == 12 || ADDITIONAL_PRINTS == 23 || ADDITIONAL_PRINTS == 123
			if(Sfix == 0) fprintf(file_ttotp,"\n");
			#endif
			for(x=0;x<(*data)->n_loci;x++) {                                
				#if ADDITIONAL_PRINTS == 1 || ADDITIONAL_PRINTS == 12 || ADDITIONAL_PRINTS == 13 || ADDITIONAL_PRINTS == 123
				fprintf(file_recp,"%G\t",postp[x][(*data)->n_iter].recp);
				#endif
				#if ADDITIONAL_PRINTS == 2 || ADDITIONAL_PRINTS == 12 || ADDITIONAL_PRINTS == 23 || ADDITIONAL_PRINTS == 123
				if(Sfix == 0) fprintf(file_ttotp,"%G\t",postp[x][(*data)->n_iter].Ttotp);
				#endif
			}
			#if ADDITIONAL_PRINTS == 1 || ADDITIONAL_PRINTS == 12 || ADDITIONAL_PRINTS == 13 || ADDITIONAL_PRINTS == 123
			fprintf(file_recp,"\n");
			#endif
			#if ADDITIONAL_PRINTS == 2 || ADDITIONAL_PRINTS == 12 || ADDITIONAL_PRINTS == 23 || ADDITIONAL_PRINTS == 123
			if(Sfix == 0) fprintf(file_ttotp,"\n");
			#endif
		}
	}
	return 0;
}
void free_matrix_sfixalltheta(struct var **data)
{
	int i;
	if((*data)->sfix_allthetas > 0) {
		for(i=0;i<(*data)->n_loci;i++) {
			free(postp[i]);
			free(thetafromS[i]);
		}
		free(postp);
		free(thetafromS);
	}
	return;
}

/*compare two double numbers in a long int list*/
int compare_(const void *i,const void *j) 
{
    if(*(double *)i < *(double *)j) return -1;
    if(*(double *)i > *(double *)j) return  1;
    return 0;
}

/*compare two numbers*/
int compare(const void *i,const void *j) 
{
    return (*((int *)i) - (*(int *)j));
}

double logPPoisson2(long int Si, double lambda)
{
    double value;
	double factln(long int);
    
	value = ((double)Si*(double)log((double)lambda) - lambda - factln(Si));
    return value;
}


/*Wall's program for calculating minimum recombination events*/
int Min_rec(int x, int segsit, int nsam, int inits)
{  /* Calculate min # rec. events */
  int a, b, c, e, gtest, flag = 0;
  int h;
  int t11,t12,t21,t22;
  int ispolnomhit(long int,int,int);
  
	if (segsit<2 || x >= (segsit-1)) return (0);
	
	for (a=x+1; a<segsit; ++a) {
		while(a < segsit) {
			if((h=ispolnomhit((long int)a,inits,nsam)) > 0) break; /*h is the frequency of the new mutation*/
			else a++;
		}            
		if(a < segsit) {
			for (b=x; b<a; ++b) {
				while(b < a) {
					if((h=ispolnomhit((long int)b,inits,nsam)) > 0) break; /*h is the frequency of the new mutation*/
					else b++;
				}
				if(b < a) {
					t21 = list[0][b];
					t11 = list[0][a];
					for (e=inits+1; e<inits+nsam; ++e) {
						if(list[e][b] != t21) t22 = list[e][b];
						if(list[e][a] != t11) t12 = list[e][a];
					}
					
					gtest = 0;
					for (e=inits; e<inits+nsam; ++e) {
						if (list[e][b] == t21 && list[e][a] == t11) {
							++gtest;
						break;
						}
					}
					for (e=inits; e<inits+nsam; ++e) {
						if (list[e][b] == t21 && list[e][a] == t12) {
							++gtest;
							break;
						}
					}
					for (e=inits; e<inits+nsam; ++e) {
						if (list[e][b] == t22 && list[e][a] == t11) {
							++gtest;
							break;
						}
					}
					for (e=inits; e<inits+nsam; ++e) {
						if (list[e][b] == t22 && list[e][a] == t12) {
							++gtest;
							break;
						}
					}
					if (gtest == 4) {
						flag = 1;
						break;
					}
				}
			}
		}
		if (flag == 1) break;
	}
	if (a >= segsit) return (0);
	else {
		c = Min_rec(a,segsit,nsam,inits);
		return (1+c);
	}
	return 0;
}

double correction_rec(double r,double f,int m) 
{
	if(f == (double)1) {
		if(m == 1) 
			return (r * (double)0.5);
		else 
			return r;
	}
	if(f == (double)1/(double)0.75) 
		return (r * (double)2/(double)3);
	if(f <= (double)1/(double)0.5) 
		return (double)0;
	
	return r;
}

int do_heter_gamma_sitesrec(double *categories,double gammashape,double poppar,long int nsites) 
{
	/*do a cummulative vector*/
	long int i;
	double gammadist(double);
	
	categories[0] = (double)0;
	
	if(gammashape <= (double)0) {/*we do simply all equal*/
		/*categories[0] = poppar/(double)nsites;*/
		for (i=1;i<nsites;i++) {
			categories[i] = (double)poppar/(double)(nsites-1);
			categories[i] += categories[i-1];
		}
	}
	else {/*gamma distribution*/
		/*categories[0] = gammadist(gammashape)/gammashape;
		categories[0] *= poppar/(double)nsites;*/
		for (i=1;i<nsites;i++) {
			categories[i] = (double)gammadist(gammashape)/(double)gammashape;
			categories[i] *= (double)poppar/(double)(nsites-1);
			categories[i] += categories[i-1];
		}
	}
	return 0;
}

int do_heter_gamma_sites(double *categories,double gammashape,double poppar,long int nsites,long int invariable) 
{
	/*do a cummulative vector*/
	long int i;
	double gammadist(double),newpoppar;
	double ran1();
		
	newpoppar = poppar * (double)nsites/(double)(nsites-invariable);
	
	if(gammashape <= (double)0) {/*we do simply all equal*/
		if(invariable > 0) {
			if(ran1() > (double)(nsites-invariable)/(double)nsites) categories[0] = (double)0;
			else categories[0] = (double)newpoppar/(double)nsites;
		}
		else categories[0] = (double)newpoppar/(double)nsites;
		for (i=1;i<nsites;i++) {
			if(invariable > 0) {
				if(ran1() > (double)(nsites-invariable)/(double)nsites) categories[i] = (double)0;
				else categories[i] = (double)newpoppar/(double)nsites;
			}
			else categories[i] = (double)newpoppar/(double)nsites;
			categories[i] += categories[i-1];
		}
	}
	else {/*gamma distribution*/
		if(invariable > 0) {
			if(ran1() > (double)(nsites-invariable)/(double)nsites) categories[0] = (double)0;
			else {
				categories[0] = (double)gammadist(gammashape)/(double)gammashape;
				categories[0] *= (double)newpoppar/(double)nsites;
			}
		}
		else {
			categories[0] = (double)gammadist(gammashape)/(double)gammashape;
			categories[0] *= (double)newpoppar/(double)nsites;
		}
		for (i=1;i<nsites;i++) {
			if(invariable > 0) {
				if(ran1() > (double)(nsites-invariable)/(double)nsites) categories[i] = (double)0;
				else {
					categories[i] = (double)gammadist(gammashape)/(double)gammashape;
					categories[i] *= (double)newpoppar/(double)nsites;
				}
			}
			else {
				categories[i] = (double)gammadist(gammashape)/(double)gammashape;
				categories[i] *= (double)newpoppar/(double)nsites;
			}
			categories[i] += categories[i-1];
		}
	}
	
	return 0;
}

long int localize_positiontop(double *categories,double valuer,long int start,long int end)
{
	long int half;
	
	half = (long int)floor((double)(start+end)/(double)2);
	if(half == start) {
		if((double)valuer < categories[half]) return half;
		else return half+1;
	}
	
	if((double)valuer < categories[half]) half = localize_positiontop(categories,valuer,start,half);
	else if((double)valuer > categories[half]) half = localize_positiontop(categories,valuer,half,end);
	
	return half;
}

