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

#include <math.h>
#include <stdlib.h>
#include <stdio.h>

void init_coef(double *p,int sample_size)
{
/*
 Tajima and Fu coefficients.
*/
	int x,n;
	double an,bn,cn;
		
	an = bn = 0.0;
	n = sample_size;

	for(x=1;x<n;x++) {
            an += 1.0/(double)x;
            bn += 1.0/((double)x*(double)x);
	}
        	
	p[0] = an;
	p[1] = bn;

	/* vt */
	p[2] = (2.0*((double)n*(double)n + (double)n + 3.0)/(9.0*(double)n*((double)n-1.0))
	        -((double)n+2.0)/(an * (double)n) + bn/(an * an)) / (an*an + bn);
	        
	/* ut */
	p[3] = (( ((double)n+1.0)/(3.0*((double)n-1.0)) - 1.0/an)/an) - p[2];
	
	/* vd* */
	p[4] = (bn/(an*an) - 2.0/(double)n *(1.0 + 1.0/an - an + an/(double)n) - 1.0/((double)n*(double)n))
	       / (an*an + bn);
	
	/* ud* */
	p[5] = ((((double)n-1.0)/(double)n - 1.0/an) / an) - p[4];
	
	/* vf* */
	p[6] = (((2.0*(double)n*(double)n*(double)n + 110.0 * (double)n*(double)n - 255.0 * (double)n + 153.0)
	        / (9.0 * (double)n * (double)n * ((double)n-1.0)) + (2.0*((double)n-1.0) *an)/ ((double)n*(double)n) 
	        - (8.0 * bn)/(double)n)) / (an*an + bn);
	
	/* uf* */
	p[7] = (((4.0*(double)n*(double)n + 19.0*(double)n + 3.0 - 12.0 * ((double)n +1.0) * (an + 1.0/(double)n))
	      / (3.0*(double)n *((double)n-1.0))) / an) - p[6];

	/* cn */
	cn = 2.0 * ((double)n*an-2.0*((double)n-1.0)) / (((double)n-1.0)*((double)n-2.0));
	
	/* vd */
	p[8] = 1.0 + (an*an/(bn+an*an)) * (cn - (((double)n+1.0)/((double)n-1.0)));
	
	/* ud* */
	p[9] = an -1.0 - p[8];
	
	/* vf */
	p[10] = (cn + (2.0*((double)n*(double)n+(double)n+3.0))/(9.0*(double)n*((double)n-1.0)) - 2.0
	       /((double)n-1.0)) / (an*an + bn);
	
	/* uf */
	p[11] = (1.0 + ((double)n+1.0)/(3.0*((double)n-1.0)) - 4*((double)n+1.0)
	       /(((double)n-1.0)*((double)n-1.0)) * (an+(1.0/(double)n) - 2.0*(double)n/((double)n+1.0)))
	       / an  -  p[10];
}
/*Tajima's D*/
double tajima_d(double k_, int S_, double *coef_taj)
{
	double an,ut,vt;
	double S_D = -10000;
        
        if(S_ == 0 || *(coef_taj+0) < 1.51) return(-10000); 
        
	an = *(coef_taj+0);
	ut = *(coef_taj+3);
	vt = *(coef_taj+2);

	S_D = (k_ - ((double)S_/an)) / (sqrt((ut*(double)S_) + (vt*(double)S_*((double)S_))));
	
	if (fabs(S_D) < 1.0E-15)
		S_D = 0.0;

	return S_D;
}
/*Tajima's D/Dmin*/
double tajima_dvsdmin(double k_, int S_, double *coef_taj,int sample_size)
{
	double an,kmin;
	double D_Dmin = -10000;
        
        if(S_ == 0 || *(coef_taj+0) < 1.51) return(-10000); 
        
	an = *(coef_taj+0);
        /*kmin = (double)S_ * (2.*(double)sample_size-1.)/((double)sample_size*(double)sample_size);*/
        kmin = (double)S_ * ((double)2/(double)sample_size);

	D_Dmin = (k_ - ((double)S_/an)) / (double)fabs(kmin - ((double)S_/an));
	
	return D_Dmin;
}
/*Fu and Li's D*/
double fl_d(int sample_size,int fr1,int S, double *coef) /* amb outgroup */
{
	double an;
	double ud,vd;
	int re,n;
	double D = -10000;
	n = sample_size;
	
	if(S == 0 || *(coef+0) < 1.5) return(-10000);
                
	re = fr1;	
	an = *(coef+0);
	
	vd = *(coef+8);
	ud = *(coef+9);
	D  = ((double)S - an*(double)re) / 
	     sqrt(ud*(double)S + vd*(double)S*(double)S);

	return D;
}
/* Fu and Li's D* */
double fl_d2(int sample_size,int fr1w,int S, double *coef) /* NO outgroup */
{
	double an;
	int n;
	double ud2,vd2;
	int rs;
	double D2 = -10000;
        
        if(S == 0 || *(coef+0) < 1.51) return(-10000);
	
	rs = fr1w;

	n = sample_size;
	an = *(coef+0);
	
	vd2 = *(coef+4);
	ud2 = *(coef+5);
	D2  = ((double)S/an - (double)rs*(((double)n-(double)1)/(double)n)) /
	      sqrt(ud2*(double)S + vd2*(double)S*(double)S);

	return D2;
}
/*Fu and Li's F*/
double fl_f(int sample_size,int fr1, int S, double pi, double *coef) /* amb outgroup */
{
	double uf,vf;
	int re,n;
	double F;
	n = sample_size;
	
	if(S == 0 || *(coef+0) < 1.5) return(-10000);

	re = fr1;		
	vf = *(coef+10);
	uf = *(coef+11);

	F  = (pi - (double)re) / sqrt(uf*(double)S + vf*(double)S*(double)S);

	return F;
}
/* Fu and Li's F* */
double fl_f2(int sample_size,int fr1w, int S, double pi, double *coef) /* NO outgroup */
{
	int n;
	double uf2,vf2;
	int rs;
	double F2;
        
        if(S == 0 || *(coef+0) < 1.51) return(-10000);
	
	rs = fr1w;
	
	n   = sample_size;
	vf2 = *(coef+6);
	uf2 = *(coef+7);
	
	F2  = (pi - ((((double)n-1.0)/(double)n)*(double)rs)) / 
	        sqrt(uf2*(double)S + vf2*(double)S*(double)S);

	return F2;
}

/*Fay and Wu H*/
double fay_wu(int sample_size,int *fr,double pi) /* nomes outgroup */
{
    int i;
    double Th,H;
    
    if(sample_size < 2) return(-10000);
    
    Th = 0.;
    for(i=1;i<sample_size;i++) Th += ((double)*(fr+i))*((double)i*(double)i);
    Th *= (double)2/((double)sample_size*(double)(sample_size-1));
    
    H = pi - Th;

    return H;
}

/*Fay and Wu H divided by their minimum given S*/
double fay_wuvsminH(int sample_size,int *fr,double pi,int S) /* nomes outgroup */
{
    int i;
    double Th,Hmin,Thmax,pimin;
    
    if(sample_size < 2) return(-10000);
    
    Th = 0.;
    for(i=1;i<sample_size;i++) Th += ((double)*(fr+i))*((double)i*(double)i);
    Th *= 2.0/((double)sample_size*(sample_size-1));
    
    Thmax  = (double)S * ((double)(sample_size-1)*(double)(sample_size-1)) * ((double)2/((double)sample_size*(double)(sample_size-1)));
    /*pimin  = (double)S * (2./((double)sample_size*(double)sample_size) * ((double)sample_size-1.);*/
    pimin  = (double)S * ((double)2/((double)sample_size));
    
    if(pimin == Thmax) return(-10000);
    else Hmin = (pi - Th)/fabs((pimin - Thmax));

    return Hmin;
}

double fay_wu_normalized(int n,int *fr,double pi) /* Fay and Wu H nomes outgroup NORMALIZED (Zeng et al. Genetics 2006 174: 1431-9)*/
{
    int i;
    double TL,H,varpiTL,thetaw,an,bn,S;
    
    if(pi == (double)0.0 || n < 2) return(-10000);
    
    TL = thetaw = an = bn = (double)0;
	for(i=1;i<n;i++) {
		TL += ((double)*(fr+i))*((double)i);
		thetaw += (double)*(fr+i); 
		an += (double)1/(double)i;
		bn += (double)1.0/((double)i*(double)i);
	}
    TL *= (double)1.0/((double)(n-(double)1));
    S = thetaw;
	thetaw = thetaw/an;
	varpiTL = thetaw * ((double)(n-(double)2))/((double)6*((double)(n-(double)1))) + 
	          S*(S-(double)1)/(an*an+bn) * 
			  ((double)18*(double)n*(double)n*((double)3*(double)n+(double)2)*(bn+(double)1.0/((double)n*(double)n)) - 
			  ((double)88*(double)n*(double)n*(double)n + (double)9*(double)n*(double)n - 13*(double)n + (double)6)) /
			  ((double)9*((double)n*(n-(double)1)*(n-(double)1)));
	
	H = (pi - TL)/(double)sqrt(varpiTL);

    return H;
}

double fay_wu_normalized2(int n,double thetaL,double thetaw,double S,double *coef,double pi) /* Fay and Wu H nomes outgroup NORMALIZED (Zeng et al. Genetics 2006 174: 1431-9)*/
{
    double H,varpiTL,an,bn;
    
    if(pi == (double)0.0 || n < 2) return(-10000);
    an = coef[0];
	bn = coef[1];
	varpiTL = thetaw * ((double)(n-(double)2))/((double)6*((double)(n-(double)1))) + 
	          S*(S-(double)1)/(an*an+bn) * 
			  ((double)18*(double)n*(double)n*((double)3*(double)n+(double)2)*(bn+(double)1.0/((double)n*(double)n)) - 
			  ((double)88*(double)n*(double)n*(double)n + (double)9*(double)n*(double)n - 13*(double)n + (double)6)) /
			  ((double)9*((double)n*(n-(double)1)*(n-(double)1)));
	
	H = (pi - thetaL)/(double)sqrt(varpiTL);

    return H;
}


double E_zeng(int n,double thetaL,double thetaw,double S,double *coef) /* (Zeng et al. Genetics 2006 174: 1431-9)*/
{
    double E,varLW,an,bn;
    
    if(thetaw == (double)0.0 || n < 2) return(-10000);
    an = coef[0];
	bn = coef[1];
	varLW = thetaw * ((double)n/((double)2*(double)(n-1)) - (double)1/an) +
			S*(S-(double)1)/(an*an+bn) * 
			(bn/(an*an) + (double)2*bn*((double)n/(double)(n-1))*((double)n/(double)(n-1)) - 
			 (double)2*((double)n*bn-(double)n+(double)1)/((double)(n-1)*an) - 
			 ((double)3*(double)n+(double)1)/((double)(n-1)));
	
	E = (thetaL - thetaw)/(double)sqrt(varLW);

    return E;
}

double Fst(double piwithin, double pibetween,int ntotpop)
{
    double fst; /*equation 3 or 6 in Hudson, Slatkin and Maddison 1992*/

    if((pibetween == 0.0 && piwithin == 0.0) || ntotpop < 2) return(-10000);
    /*fst = 1. - piwithin/( piwithin/ntotpop + (1.-1./ntotpop)*pibetween);*//* eq. 6 */
    fst = (double)1 - (piwithin/pibetween); /* eq. 3 */
    return(fst);
}

/*Ewens-Watterson test*/
double EWtest(int Nsample, int *Freqhap)
{
    int i;
    double H;
    
    H = 0;
    for(i=0;i<Nsample;i++)
        H += Freqhap[i] * Freqhap[i];
    H = (double)H/((double)Nsample*(double)Nsample);
    
    return H;
}

double pwh(int n,double *pid)
{/*ty to find two  divergent groups in the mismatch distribution with low variation within in at least one*/
	int ncomp;
	int i;
	double within,between;
	double diff=(double)0;
	int compare_(const void *,const void *);
	
	ncomp = n * (n-1) / 2;
	/*the pairwise vector has to be sorted*/
	qsort(pid,(int)ncomp,sizeof(double),compare_);/*all values are sorted*/
	/*divide sorted vector in two. find the highest difference*/
	for(i=0;i<ncomp-1;i++) {
		between= pid[i+1]-pid[i];
		within = fabs((pid[ncomp-1]-pid[i+1]) - (pid[i]-pid[0]));
		if((between > within) && (between + within > diff)) 
			diff = between + within;
	}
	return diff;
}

/*Haplotype tests from Depaulis et al.*/
double testHap(int Nsample, int *Freqhap)
{/*Depaulis statistics*/
    int i;
    double H;
    
    H = 0;
    for(i=0;i<Nsample;i++)
        H += Freqhap[i] * Freqhap[i];
    H = (double)1 - H/((double)Nsample*(double)Nsample);
    /*and weighted: haplotype diversity*/
    H = H*(double)Nsample/(double)(Nsample-1);
    
    return H;
}

/*Fs from Fu */
double Fs(int Nsample, double pi, int NumAlelos)
{
    /* Rozas program */
	
    double SumaP;
    double RestaP;
    int AleloI;
    long int i;       
    double ValorFs;
    double *qew;
    double est_var;
    double FunEq23Ewens(int, int, double, double *);

    if(pi == 0.0 || Nsample < 2) return(-10000);	
    est_var = pi;
    qew  = (double *)malloc((long int)Nsample*(long int)Nsample*sizeof(double));
    
    for(i=0;i<(long int)Nsample*(long int)Nsample;i++)
    	qew[i] = -1.0;
            
    SumaP=RestaP=0.0;
    for (AleloI=1;AleloI<NumAlelos;AleloI++) {
        /* calculo q(n,aleloI)   ecuacion 21 (recurrente con eq. 19 y 20) */
        SumaP += FunEq23Ewens(Nsample, AleloI, est_var, qew);
    }

    if(SumaP > 1.-1E-37) {
    	for (AleloI = NumAlelos;AleloI <= Nsample; AleloI++)
            RestaP += FunEq23Ewens(Nsample, AleloI, est_var, qew);	 	
        if(RestaP < 1E-37)
            return -10000;
        ValorFs = log((double)RestaP) - log((double)1.0-RestaP);
    }
    else {
        if(SumaP < 1E-37)
            return +10000;
        else
            ValorFs = log((double)1.0-(double)SumaP) - log((double)SumaP);
    }
    if (fabs(ValorFs) < 1.0E-15)
        ValorFs = 0.0;	    
    
    free(qew);
    return ValorFs;
}

double FunEq23Ewens(int N,int i,double theta, double *qew_)
{                  
    long int acceso; 
    int jj;
    double ValorN;  /* log del numerador */
    double ValorD;  /* log del denominador */

    acceso= (long int)(N-1) * (long int)N + (long int)i - (long int)1;    
    ValorN=0.0;
    ValorD=0.0;        
    if (qew_[acceso] < 0.0) {   
        if (i==1) {
            /* calculo de qj,1   (i = 1)   Antigua equacion 19  */
            if(N > 2) {             
                for (jj=2;jj<N;jj++)
                    ValorN = ValorN + log((double)jj);  
            }
            ValorN = ValorN + log(theta);
            for(jj=0;jj<N;jj++)
                ValorD  = ValorD + log((double)theta + (double)jj);      
            qew_[acceso] = exp((double)(ValorN - ValorD)); 
        }    
        if(i==N) {          
            /* calculo de qj,j   (n = i)   antigua equacion 20 */
            ValorN = log((double)theta) * (double)N;
            for(jj=0;jj<N;jj++)     
                ValorD  = ValorD + log((double)theta + (double)jj);
            qew_[acceso] = exp((double)(ValorN - ValorD));
	}
	if(i>1 && i<N) {    
            /*  recursividad  */
            qew_[acceso] = FunEq23Ewens(N-1,i,  theta,qew_) * ((double)(N-1)/(theta + (double)N-1.0))
                         + FunEq23Ewens(N-1,i-1,theta,qew_) *         (theta/(theta + (double)N-1.0));
        }    
    }  
    return(qew_[acceso]);
}

/* Programa de Ying */
double estnm(int npop,int nsam,int *config,long int segsites,char **list)
{
    int i;
    double gst,within,between;
    void diff_within_between(int,long int,int,int *,char **,double *,double *);

    nsam = 0 ;
    for(i=0;i<npop; i++) nsam += config[i] ;

    diff_within_between(nsam,segsites,npop,config,list,&within,&between);
    if(within == 0) gst = 0.;
    else gst = 1. - within/(within/npop + (1.-1./npop)*between);
    return (gst);
}

void diff_within_between(int nsam,long int ns,int npop,int *config,char **list,double *pwithin,double *pbetween)
{
    int pop,ind,ind2,indpop,endpop,count ;
    double diff(long int,char *,char *);
	
    *pwithin = *pbetween = 0. ;
    for(pop=indpop=ind=endpop=count=0; pop<npop; pop++,ind++)
        for(endpop += config[pop];ind < endpop-1; ind++) 
            for(indpop = 1; (indpop+ind)<endpop; indpop++) {
                *pwithin += diff(ns,list[ind],list[ind+indpop]);
                count++;
            }
    *pwithin /= count ; /*Es divideix pel conjunt de totes les combinacions. no es separa el valor de cada subpoblacio*/
    for(pop=ind=count=endpop=0;pop<npop-1;pop++)
        for(endpop += config[pop]; ind<endpop; ind++)
            for(ind2 = endpop ; ind2<nsam ; ind2++) {
                *pbetween += diff(ns,list[ind],list[ind2]);
                count++;
            }
    *pbetween /= count ;
}

double diff(long int ns,char *gam1,char *gam2)
{
    long int  i;
    double count;
	
    count = 0.;
    for(i=0;i<ns;i++) if(gam1[i] != gam2[i]) count += 1.;
    return(count);
}

/*R2 Ramos & Rozas: "*unic" is the number of singletons in each sequence (in comparison to the sample studied)*/	
double R2(long int *unic,double pi,int sample_size,long int S)
{
    double sm2 = 0.0;
    int i;
    
    if(S == 0 || sample_size == 0) return(-10000);
    for (i=0;i<sample_size;i++)
            sm2 += ((double)unic[i] - pi/2.0)*((double)unic[i] - pi/2.0);
    
    sm2 = sqrt(sm2/((double)sample_size))/(double)S;
            
    if (sm2 < 1.0E-15)
            sm2 = 0.0;

    return (double)sm2;
}

/*Gxi(outgroup) test from Fu(1995) based on frequency of segregating sites, but using theta from watterson*/
double Gxi(int sample_size,int *fr, double thetaw) /* with outgroup */
{
    int i;
    double G;
    double varxi(int,int,double);
    
    if(thetaw == 0.0 || sample_size < 2) return(-10000);

    G = 0.;
    for(i=1;i<sample_size;i++) 
        G += ((fr[i] - thetaw/(double)i) * (fr[i] - thetaw/(double)i)) / varxi(i,sample_size,thetaw);
    G /= ((double)sample_size - (double)1);

    return G;
}

double varxi(int i,int sample_size,double theta)
{
    double sigma(int,int);
    return (double)1/((double)i)*theta + sigma(i,sample_size)*theta*theta;
}

double sigma(int i,int sample_size)
{
    double ai(int);
    double bn(int,int);
    
    if(i <  (double)sample_size/(double)2) return bn(i+1,sample_size);
    if(i == (double)sample_size/(double)2) return (double)2*(ai(sample_size) - ai(i))/((double)(sample_size - i)) - (double)1/((double)(i*i));
    if(i >  (double)sample_size/(double)2) return bn(i,sample_size) - (double)1/((double)(i*i));
    
    return -10000;
}

double ai(int i)
{
    int j;
    double a = 0;
    for(j=1;j<i;j++) a += (double)1/(double)j;
    return a;
}

double bn(int i,int sample_size)
{
    double ai(int);
    return ((double)2*(double)sample_size * (ai(sample_size + 1) - ai(i))) / (((double)sample_size - (double)i + (double)1) * ((double)sample_size - (double)i)) - 
            (double)2/((double)sample_size - (double)i);
}

/*Gximod (outgroup) test from Fu(1995) based on frequency of segregating sites, but using theta from watterson,
 and not variance...*/
double Gximod(int sample_size,int *fr, double thetaw) /* with outgroup */
{
    int i;
    double G;
    
    if(thetaw == 0.0 || sample_size < 2) return(-10000);

    G = 0.;
    for(i=1;i<sample_size;i++) 
        G += ((fr[i] - thetaw/(double)i) * (fr[i] - thetaw/(double)i))/((double)ceil(thetaw/(double)i));

    return G;
}

double frabs(int sample_size,int *fr, double thetaw) /* with outgroup */
{
    int i;
    double G;
    
    if(thetaw == 0.0 || sample_size < 2) return(-10000);

    G = 0.;
    for(i=1;i<sample_size;i++) 
        G += (double)fabs((double)fr[i] - thetaw/(double)i);

    return G/(double)sample_size;
}

