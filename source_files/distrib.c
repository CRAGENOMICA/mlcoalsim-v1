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
#define PI 3.14159265359

#include <stdio.h>
#include <stdlib.h>

/*Routines based on Numerical Recipes in C*/

double gammalogn(double zz)
{
	/*Based on Numerical Recipes in C, Press et al. 1992. p. 213. and on 
	Lanczos 1964 J. SIAM Numer. Anal. Ser. B, Vol. 1 pp 86-96.*/
	
	/*gamma distribution for a z integer*/
	double loggammaz;
	double z,logg,h,sumc;
	static double gamma = 5.0;
	static double c0 =  1.000000000178;
	static double c1 = 76.180091729406;
	static double c2 = 86.505320327112;
	static double c3 = 24.014098222230;
	static double c4 =  1.231739516140;
	static double c5 =  0.001208580030;
	static double c6 =  0.000005363820;
	
	if(zz <= 0.) {
		puts("Error gamma");
		return (double)-10000.;
	}
	
	z = (double)zz;
	h = (double)sqrt(2. * PI);
	sumc = c0 + c1/(z+1.) - c2/(z+2.) + c3/(z+3.)  - c4/(z+4.) + c5/(z+5.) - c6/(z+6.);
	logg = (z + 0.5)*(double)log((double)(z + gamma + 0.5)) - (z + gamma + 0.5);
	loggammaz = log((double)h);
	loggammaz += logg + log((double)sumc);
	loggammaz -= log((double)z);
	
	return (double)loggammaz;
}

double factln(long int x)
{
	/*Based on Numerical Recipes in C, Press et al. 1992 and on
	Lanczos 1964 J. SIAM Numer. Anal. Ser. B, Vol. 1 pp 86-96.*/
	
	/*G(n+1) = n!*/
	/*do the log(n!)*/
	double gammalogn(double);
	static double factlog[120];
		
	/*
	if(x < 0) { 
		puts("Error factln");
		return (double)-10000.;
	}
	*/
	if(x == 0) return 0.;
	
	if(x < 120) {
		if(factlog[x] == (double)0) {
			factlog[x] = gammalogn((double)x+(double)1.0);
			return factlog[x];
		}
		else return factlog[x];
	}
	return (gammalogn((double)x+(double)1.0));
}

double gammadist(double alfa) 
{
	/*Based on Numerical Recipes in C, Press et al. 1992, on
	Cheng anf Feast 1979, Appl. Statist. 28, No. 3, pp. 290-295, and on
	PSeq-Gen v1.1. Nicholas C. Grassly, Jun Adachi and Andrew Rambaut*/
	
	double ran1();
	double rndgamma1 (double);
	double rand0,rand1,rand2;
	double a,b,c,d,f,W;
	
	if(alfa <= (double)0) {
		return (double)-10000.0;
	}
	if(alfa < (double)1) {
		a = rndgamma1(alfa);
		return a;
	}
	if(alfa == (double)1.) return (-(double)log((double)ran1()));
	
	a = (double)alfa - (double)1.0;
	b = (alfa - (double)1./((double)6.*(double)alfa))/a;
	c = (double)2./a;
	d = c + (double)2.0;
	f = (double)sqrt((double)alfa);
		
	do {
		if(alfa < (double)3) {
			rand1 = ran1();
			rand2 = ran1();
		}
		else {
			do {
				rand0 = ran1();
				rand1 = ran1();
				rand2 = rand1 + (double)1/f * ((double)1-(double)1.86*rand0);
			}while(rand2 < (double)0 || rand2 > (double)1.0);
		}
		W = b * rand1/rand2;
		if(c*rand2-d+W+(double)1./W <= (double)0.) break;
	}while(c*log((double)rand2)-log((double)W)+W-(double)1. >= (double)0.);
	
	return a*W;
}

/*From PSeq-Gen v1.1. Nicholas C. Grassly, Jun Adachi and Andrew Rambaut */
double rndgamma1 (double s)
{

	double ran1(void);
	double			r, x=0.0, small=1e-37, w;
	static double	a, p, uf, ss=10.0, d;
	
	if (s!=ss) {
		a  = 1.0-s;
		p  = a/(a+s*exp(-a));
		uf = p*pow(small/a,s);
		d  = a*log(a);
		ss = s;
	}
	for (;;) {
		r = (double)ran1();
		if (r > p) {
			x = a-log((1.0-r)/(1.0-p));
			w=a*log(x)-d;
		}
		else {
			if (r>uf) {
				x = a*pow(r/p,1/s);
				w=x;
			}
			else return ((double)0.0);
		}
		r = (double)ran1();
		if (1.0-r <= w && r > 0.0)
			if (r*(w+1.0) >= 1.0 || -log(r) <= w)  
				continue;
		break;
	}
	return ((double)x);
}

double poissondist(double lambda) 
{
	/*Based on Atkinson 1979 Appl. Statist. 28: No. 1, pp, 29-35.*/
	double ran1(void);	
	double r,s;
	int N;
	
	double factln(long int);
	double alfa,beta,k,X;
	double rand1,rand2;
	static double c = (double)0.6;
	
	if(lambda < (double)0) {
		puts("Error poissondist");
		return (double)-10000.;
	}
	
	if(lambda == (double)0) return (double)0;
	if(lambda <= (double)20) {
		r = (double)exp(-(double)lambda);
		N = (int)0;
		s = (double)1;
		do {
			s *= ran1();
			if(s >= r) N += (int)1;
			else break;
		}while(1);
	}
	else {
		beta = (double)PI * (double)1./(double)sqrt((double)3*(double)lambda);
		alfa = beta * (double)lambda;
		k = (double)log((double)c) - (double)lambda - (double)log((double)beta);
		do{
			rand1 = ran1();
			X = (alfa-(double)log(((double)1-(double)rand1)/(double)rand1))/beta;
			if(X >= (double)-0.5) {
				N = (int)(X + (double)0.5);
				rand2 = ran1();
				if(alfa - beta*X +
				  (double)log((double)(rand2/(((double)1+(double)exp((double)(alfa-beta*X)))*((double)1+(double)exp((double)(alfa-beta*X)))))) 
				  <= k + (double)N*(double)log((double)lambda) - (double)factln((long int)N)) 
					break;
			}
		}while(1);
	}
	return (double)N;
}

double binomialdist(double pp, int n) 
{
	/*Based on Numerical Recipes in C, Press et al. 1992 and on
	Fishman 1979 J. American Statistical Association Vol. 74, No. 366, pp 418-423*/
	
	double ran1(void);	
	double p,np;
	int N;
	int nn;
	double r;
	double A,B,C,D,V,s;
	double m,mu;	
	static double *f=0;
	double poissondist(double);
	static int max = 200;
	
	if(f == 0) {
		if((f=(double *)calloc(max+1,sizeof(double))) == 0) {
			printf("Error allocation memory");
			return -10000;
		}
		f[1] = (double)0;
		for(N=1;N<max;N++)
			f[N+1] = f[N] + (double)log((double)N);
	}
	if(n > max) {
		if((f=(double *)realloc(f,n*sizeof(double))) == 0) {
			printf("Error allocation memory");
			return -10000;
		}
		for(N=max;N<n;N++)
			f[N+1] = f[N] + (double)log((double)N);
		max = n;
	}
	
	if(pp > 0.5) p = (double)1.-pp;
	else p = pp;
	
	np = n * p;
	
	if(n==0) {
		puts("Error bindist");
		return (double)-10000.;
	}
	if(p==(double)0) {
		if(pp > 0.5) return (double)n;
		return (double)0;
	}
	
	if(n < 20) {
		/*Bernouilli Method*/
		nn = n;
		N=0;
		while(nn--) 
			if(ran1()<p) N++;		
	}
	else {
		if(np < (double)10) {
			/*Rejection Method: BI Algorithm*/
			s = (double)1- p;
			A = (double)1;
			B = p/s;
			C = ((double)n+(double)1)*B;
			D = A;
			N = 0;
			V = ran1()/(double)pow(s,(double)n);
			while(V > A) {
				N++;
				D *= (C/(double)N - B);
				A += D;
				if(N > n) break;
			}
		}
		else {
			/*Poisson method: BP Algorithm*/
			mu = n - (double)floor((double)(n*((double)1 - p)));
			if(n*((double)1-p) - (double)floor((double)(n*((double)1-p))) > p)
				mu = p*((double)floor((double)(n*((double)1-p))) + (double)1) / ((double)1-p);
			r = ((double)1/p - (double)1) * mu;
			s = (double)log((double)r);
			m = (double)floor((double)(r));
			do {
				do {
					N = (int)poissondist(mu);
				}while((int)N > n);
				V = -(double)log((double)ran1());
			}while(V < (m-(double)(n - N))*s - f[(int)m+1] + f[(int)(n-N)+1]);
		}
	}
	if(pp > 0.5) N = n - N;
	return (double)N;
}

int zbracn(double (*func)(double,double,double,double,double,double,double,double,double,double),double *x1,double *x2,double a0,double a1,double a2,double a3,double a4,double a5,double a6,double a7,double a8)
{
	/* Based on Numerical Recipes in C. Press et al. 1992. 
    We need *x1 and *x2 be the range where a root is within them. We expand geometrically the range until finding
	(one a positive and one a negative value). If not, return 0.
	*/

    double f1,f2;
    int k=60;
    
    if(*x1 == *x2) return 0;

    f1 = (*func)(*x1,a0,a1,a2,a3,a4,a5,a6,a7,a8);
    f2 = (*func)(*x2,a0,a1,a2,a3,a4,a5,a6,a7,a8);
	
	if(f1*f2 < (double)0) return 1;

    while(k--) {
        if(fabs(f1) < fabs(f2)) {
            *x1 += (double)1.5 * (*x1 - *x2);
            f1 = (*func)(*x1,a0,a1,a2,a3,a4,a5,a6,a7,a8);
        }
        else {
            *x2 += (double)1.5 * (*x2 - *x1);
            f2 = (*func)(*x2,a0,a1,a2,a3,a4,a5,a6,a7,a8);
        }
        if(f1*f2 < (double)0) return 1;
    }
    return 0;
}

double zriddrn(double (*func)(double,double,double,double,double,double,double,double,double,double),double xlow,double xhigh,double xacc,double a0,double a1,double a2,double a3,double a4,double a5,double a6,double a7, double a8)
{
	/* Based on Numerical Recipes in C. Press et al. 1992., p. 358 an on
	Ridders, 1979, IEEE Transactions on Circuits and systems, Vol. Cas-26, No. 11, pp. 979-980.
	*/
    int k=60;
	double flow,fhigh;
	double f1,f2,f3,f4;
	double x1,x2,x3,x4;
	double den,num,nsign;
	

    flow  = (*func)(xlow ,a0,a1,a2,a3,a4,a5,a6,a7,a8);
    fhigh = (*func)(xhigh,a0,a1,a2,a3,a4,a5,a6,a7,a8);

	if(flow  == (double)0) return xlow;
	if(fhigh == (double)0) return xhigh;
	if(flow*fhigh > (double)0) 
		return (double)-1e32;
	
	x1 = xlow;
	x2 = xhigh;
	f1 = flow;
	f2 = fhigh;
		
	while(k--) {
		x3 = (x1+x2)/(double)2;
		f3 = (*func)(x3,a0,a1,a2,a3,a4,a5,a6,a7,a8);
		if(f1 - f2 < (double)0) nsign = (double)-1;
		else nsign = (double)1;
		num = (x3-x1) * f3 * nsign;
		den = (double)sqrt((double)f3*(double)f3 - (double)f1*(double)f2);
		if(den <= xacc && -den <= xacc) return x3;
		x4 = x3 + num/den;
		f4 = (*func)(x4,a0,a1,a2,a3,a4,a5,a6,a7,a8);
		if(f4 <= xacc && -f4 <= xacc) return x4;
		if(f3*f4<(double)0) {
			x1 = x3;
			f1 = f3;
			x2 = x4;
			f2 = f4;
		}
		else {
			if(f1*f4<(double)0) {
				x2 = x4;
				f2 = f4;
			}
			else {
				if(f2*f4<(double)0) {
					x1 = x4;
					f1 = f4;
				}
			}
		}
		if(fabs(x1-x2) <= xacc) return x1;
	}	
	return (double)-1e32;
}



