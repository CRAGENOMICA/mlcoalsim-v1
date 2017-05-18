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
#include <math.h>
#include <time.h>
#include "mlsp_sm.h"

int main (int argc,char *argv[])
{	    
    void input_data(FILE *,struct var **);
    void output_data(FILE *, char *, struct var **);
    void getpars_fix(struct var **, struct var2 **);
    void getpars_mod(struct var **, struct var2 **,int);
    void free_getpars_fix(struct var **, struct var2 **);
    void free_inputdata(struct var **);
	void free_matrix_sfixalltheta(struct var **);
	int print_matrix_sfixalltheta(struct var **,FILE *,FILE *,FILE *);
	int print_matrix_rmfix(struct var **,FILE *,FILE * /*,FILE * */,int);
    
    int print_neuttest(struct var **,FILE *,char *);
    int ms(struct var2 **,FILE *);    
    
    struct var *data;
    struct var2 *inputp;
    FILE *file_input;
    FILE *file_output;
	FILE **file_outputm;
	char numchar[5];
    char file_in[400];
    char file_out[400];
    int x,y,j;
	#if ADDITIONAL_PRINTS
	char namefile_[420];
	char *el;
	#endif
	FILE *file_thetap=0;
	FILE *file_ttotp=0;
	FILE *file_Srange=0;
	FILE *file_recp=0;
	/*
	FILE *file_recrange=0;
	*/
	#if TIMEDURATION
    time_t start;
    time_t end;
    struct tm *date;
    char s[80];
    int hour,min,sec;
	#endif
	int a,b,c,totalnloci;
    char names[NOBS_STATISTICS][20]  = {{"TD"},{"Fs"},{"FDn"},{"FFn"},{"FD"},{"FF"},{"H"},
                                    {"B"},{"Q"},{"ZnA"},{"Fst"},{"Kw"},{"Hw"},{"R2"},
                                    {"S"},{"piw"},{"pib"},{"thetaWatt"},{"ThetaTaj"},
                                    {"ThetaFW"},{"D_Dmin"},{"Hnorm"},{"maxhap"},{"maxhap1"},{"Rm"},
									{"ThetaFL"},{"ThetaL"},{"ZengE"},{"EWtest"},{"Fstw"},{"Pwh"}};

	puts(MLCOALSIM);
    /******************************* INPUT/OUTPUT FILES ***************************************/
    if(argc != 1 && argc != 3) {
        puts("usage: input_file output_file");
        exit(1);
    }
    if(argc == 3) {
        strcpy(file_in,argv[1]);
        strcpy(file_out,argv[2]);
    }    
    if(argc == 1) {
        puts("\nInput file?");
        scanf("%s",file_in);
    }

    if (!(file_input = fopen (file_in,"r"))) perror("Error in input/output");
    input_data(file_input,&data);
    fclose(file_input);

    if(argc == 1 /*&& (*data).likelihood_line == 0*/) {
        puts("\nOutput file?");
        scanf("%s",file_out);
		if(strrchr(file_out,'.') == 0) {
			printf("\n Sorry. The output file must contain an extension. \n");
			exit(1);
		}
    }

	#if TIMEDURATION
    /*Starting time*/
    time(&start);
    date = localtime(&start);
    strftime(s,80,"%c",date);
	#endif

	if((*data).neutral_tests == 0 && (*data).n_loci > 1 && (*data).likelihood_line == 0) {
		if(!(file_outputm = (FILE **)calloc((*data).n_loci,sizeof(FILE *)))) perror("Error in input/output");
		for(x=0;x<(*data).n_loci;x++) {
			namefile_[0] = '\0';
			strncat(namefile_,file_out,420);
			el = strrchr(namefile_,'.');
			*el = '\0';
			strncat(namefile_,"locus\0",420);
			j = x;
			numchar[0] = (int)(j/1000) + 48; 
			if((int)(j/1000) > 0) j -= ((int)(j/1000))*1000;
			numchar[1] = (int)(j/100) + 48;
			if((int)(j/100) > 0) j -= ((int)(j/100))*100;
			numchar[2] = (int)(j/10) + 48;
			if((int)(j/10) > 0) j -= ((int)(j/10))*10;
			numchar[3] = (int)(j/1) + 48;
			numchar[4] = '\0';
			strncat(namefile_,numchar,420);
			strncat(namefile_,".out\0",420);
			if(!(file_outputm[x] = fopen(namefile_,"w"))) perror("Error in input/output");
			fputs(MLCOALSIM,file_outputm[x]);
			output_data(file_outputm[x],file_in,&data);
		}
	}
	else {
		if((*data).likelihood_line == 0) {
			if(!(file_output = fopen (file_out,"w"))) perror("Error in input/output");
			fputs(MLCOALSIM,file_output);
			output_data(file_output,file_in,&data);
		}
		else {
			if(!(file_output = fopen (file_out,"r"))) {
				if(!(file_output = fopen (file_out,"a+"))) perror("Error in input/output");
				c = NOBS_STATISTICS;
				if((*data).linked == 1 && (*data).despl > 0 && (*data).window > 0) {
					totalnloci = (int)ceil(((double)(*data).nsites[1] - (double)(*data).window) / ((double)(*data).despl)) + (int)1;
				}	
				else {
					if((*data).linked > 1) totalnloci = (*data).linked;
					else totalnloci = (*data).n_loci;
				}
				if((*data).n_iter > 1) {
					fprintf(file_output,"Total\t");
					for(a=0;a<c;a++) {
						if((*data).obs_statistics[a][1] == 1) {
							fprintf(file_output,"%s[Total]\t",names[a]);
						}
					}
				}
				for(a=0;a<c;a++) {
					if((*data).obs_statistics[a][1] == 1) {
						for(b=0;b<(int)totalnloci;b++) {
							fprintf(file_output,"%s[%d]\t",names[a],(int)b);
						}
					}
				}
				fprintf(file_output,"\n");
			}
			else {
				fclose(file_output);
				if(!(file_output = fopen (file_out,"a"))) perror("Error in input/output");
			}
			fflush(file_output);
		}
	}
	/******************************** OUTPUT HEADER ********************************************/
    if((*data).pr_matrix == 0 && (*data).likelihood_line == 0) {
        /*/if(!(*data).neutral_tests) fputs("\nlocus\tmig_rate\text_rate\t",file_output);*/
        if((*data).neutral_tests) {
            fputs("Neutral tests (excluding multiple hits positions)\n",file_output);
			if((*data).n_loci > 2 || (*data).linked > 2)
				fputs("avg(TD)\tvar(TD)\tavg(Fs)\tvar(Fs)\tavg(FD*)\tvar(FD*)\tavg(FF*)\tvar(FF*)\tavg(FD)\tvar(FD) \tavg(FF)\tvar(FF)\tavg(H)\tvar(H)\tavg(B)\tvar(B)\tavg(Q)\tvar(Q)\tavg(ZA)\tvar(ZA)\tavg(Fst)\tvar(Fst)\tavg(Kw)\tvar(Kw) \tavg(Hw)\tvar(Hw)\tavg(R2)\tvar(R2)\tavg(S)\tvar(S)\tavg(pi_w)\tvar(pi_w)\tavg(pi_b)\tvar(pi_b)\tavg(ThetaWatt)\tvar(thetaWatt)\tavg(ThetaTaj)\tvar(ThetaTaj)\tavg(ThetaFW)\tvar(ThetaFW)\tavg(D/Dmin)\tvar(D/Dmin)\tavg(Hnorm)\tvar(Hnorm)\tavg(maxhap)\tvar(maxhap)\tavg(maxhap1)\tvar(maxhap1)\tavg(Rm)\tvar(Rm)\tavg(theta_fl)\tvar(theta_fl)\tavg(theta_l)\tvar(theta_l)\tavg(zeng_E)\tvar(zeng_E)\tavg(EW-test)\tvar(EW-test)\tavg(Fstw)\tvar(Fstw)\tavg(Pwh)\tvar(Pwh)",file_output);
			else {
				if((*data).n_loci < 2 && (*data).linked < 2)
					fputs("value(TD)\tvalue(Fs)\tvalue(FD*)\tvalue(FF*)\tvalue(FD)\tvalue(FF)\tvalue(H)\tvalue(B)\tvalue(Q)\tvalue(ZA)\tvalue(Fst)\tvalue(Kw)\tvalue(Hw)\tvalue(R2)\tvalue(S)\tvalue(pi_w)\tvalue(pi_b)\tvalue(ThetaWatt)\tvalue(ThetaTaj)\tvalue(ThetaFW)\tvalue(D/Dmin)\tvalue(Hnorm)\tvalue(maxhap)\tvalue(maxhap1)\tvalue(Rm)\tvalue(theta_fl)\tvalue(theta_l)\tvalue(zeng_E)\tvalue(EW-test)\tvalue(Fstw)\tvalue(Pwh)",file_output);
				else 
					fputs("avg(TD)\tavg(Fs)\tavg(FD*)\tavg(FF*)\tavg(FD)\tavg(FF)\tavg(H)\tavg(B)\tavg(Q)\tavg(ZA)\tavg(Fst)\tavg(Kw)\tavg(Hw)\tavg(R2)\tavg(S)\tavg(pi_w)\tavg(pi_b)\tavg(ThetaWatt)\tavg(ThetaTaj)\tavg(ThetaFW)\tavg(D/Dmin)\tavg(Hnorm)\tavg(maxhap)\tavg(maxhap1)\tavg(Rm)\tavg(theta_fl)\tavg(theta_l)\tavg(zeng_E)\tavg(EW-test)\tavg(Fstw)\tavg(Pwh)",file_output);
			}
			fputs("\n",file_output);            
        }
    }
	/******************************** OUTPUT Sfixall_thetas ***************************************/
    if(((*data).ifgamma == 1 && (*data).likelihood_line == 0) || ((*data).range_thetant && (*data).likelihood_line == 0)) {
        #if ADDITIONAL_PRINTS == 1 || ADDITIONAL_PRINTS == 12 || ADDITIONAL_PRINTS == 13 || ADDITIONAL_PRINTS == 123
        namefile_[0] = '\0';
		strncat(namefile_,file_out,420);
		el = strrchr(namefile_,'.');
		*el = '\0';
		strncat(namefile_,"_thetapost.out",420);
		if (!(file_thetap = fopen(namefile_,"w"))) perror("Error in input/output");
        #endif
        #if ADDITIONAL_PRINTS == 3 || ADDITIONAL_PRINTS == 13 || ADDITIONAL_PRINTS == 23 || ADDITIONAL_PRINTS == 123
        namefile_[0] = '\0';
		strncat(namefile_,file_out,420);
		el = strrchr(namefile_,'.');
		*el = '\0';
		strncat(namefile_,"_thetaprior.out",420);
		if (!(file_Srange = fopen(namefile_,"w"))) perror("Error in input/output");
        #endif		
	}
	/******************************** OUTPUT rmfix ***************************************/
    if(((*data).ifgammar == 1 && (*data).likelihood_line == 0) || ((*data).range_rnt && (*data).likelihood_line == 0)) {
        #if ADDITIONAL_PRINTS == 1 || ADDITIONAL_PRINTS == 12 || ADDITIONAL_PRINTS == 13 || ADDITIONAL_PRINTS == 123
        namefile_[0] = '\0';
		strncat(namefile_,file_out,420);
		el = strrchr(namefile_,'.');
		*el = '\0';
		strncat(namefile_,"_recpost.out",420);
		if (!(file_recp = fopen(namefile_,"w"))) perror("Error in input/output");
        #endif
		/*
        #if ADDITIONAL_PRINTS == 3 || ADDITIONAL_PRINTS == 13 || ADDITIONAL_PRINTS == 23 || ADDITIONAL_PRINTS == 123
        namefile_[0] = '\0';
		strncat(namefile_,file_out,420);
		el = strrchr(namefile_,'.');
		*el = '\0';
		strncat(namefile_,"_recprior.out",420);
		if (!(file_recrange = fopen(namefile_,"w"))) perror("Error in input/output");
        #endif
		*/		
	}
    if((((*data).ifgamma == 1  && (*data).likelihood_line == 0) || (((*data).range_thetant)  && (*data).likelihood_line == 0)) || (((*data).ifgammar == 1  && (*data).likelihood_line == 0) || ((*data).range_rnt && (*data).likelihood_line == 0))) {
        #if ADDITIONAL_PRINTS == 2 || ADDITIONAL_PRINTS == 12 || ADDITIONAL_PRINTS == 23 || ADDITIONAL_PRINTS == 123
        namefile_[0] = '\0';
		strncat(namefile_,file_out,420);
		el = strrchr(namefile_,'.');
		*el = '\0';
		strncat(namefile_,"_treelength.out",420);
		if (!(file_ttotp = fopen(namefile_,"w"))) perror("Error in input/output");
        #endif
	}
    /**************** COALESCENT SIMULATIONS FOR EACH INDEPENDENT LOCI **************************/
    if(!(inputp = (struct var2 *)calloc(1,sizeof(struct var2)))) perror("calloc error.main.0");    

	#if SHOWPROGRESS == 1
    printf("\n Starting coalescent simulations...");
    /*printf("\n Warning: in case conditioning for mutation parameter, the first dot may take long time.");*/
    printf("\n Each dot indicate that aproximately 2%% of the simulation is done.");
    printf("\n         1    2    3    4    5    6    7    8    9  100%%");
    printf("\n RUN ");
    fflush(stdout);
	#endif

    getpars_fix(&data,&inputp);
    y = 0;
	for(x=0;x<(*data).n_loci;x++) {                                
        getpars_mod(&data,&inputp,x);        
		if((*data).neutral_tests == 0 && (*data).n_loci > 1) {
			if(ms(&inputp,file_outputm[x])) {
				y = 1;
				break;
			}
		}
		else {
			if(ms(&inputp,file_output)) {
				y = 1;
				break;
			}
		}
    }

	#if SHOWPROGRESS == 1
    printf(" Simulation finished.");
    fflush(stdout);
	#endif

    /******************************** OUTPUT Sfixall_thetas  ************************************/
    #if SHOWPROGRESS == 1
	printf("\n Saving results in output file/s...");
    fflush(stdout);
	#endif

    if(((*data).ifgamma == 1  && (*data).likelihood_line == 0) || ((*data).range_thetant && (*data).likelihood_line == 0)) {
		if(print_matrix_sfixalltheta(&data,file_thetap,file_ttotp,file_Srange)) return 1;  /*exit*/
		#if ADDITIONAL_PRINTS == 1 || ADDITIONAL_PRINTS == 12 || ADDITIONAL_PRINTS == 13 || ADDITIONAL_PRINTS == 123
		fclose(file_thetap);
		#endif
		#if ADDITIONAL_PRINTS == 3 || ADDITIONAL_PRINTS == 13 || ADDITIONAL_PRINTS == 23 || ADDITIONAL_PRINTS == 123
		fclose(file_Srange);
		#endif
	}
    /******************************** OUTPUT Sfixall_thetas  ************************************/
    if(((*data).ifgammar == 1  && (*data).likelihood_line == 0) || ((*data).range_rnt && (*data).likelihood_line == 0)) {
		if(print_matrix_rmfix(&data,file_recp,file_ttotp/*,file_recrange*/,((*data).ifgamma == 1 || (*data).range_thetant))) return 1;  /*exit*/
		#if ADDITIONAL_PRINTS == 1 || ADDITIONAL_PRINTS == 12 || ADDITIONAL_PRINTS == 13 || ADDITIONAL_PRINTS == 123
		fclose(file_recp);
		#endif
		/*
		#if ADDITIONAL_PRINTS == 3 || ADDITIONAL_PRINTS == 13 || ADDITIONAL_PRINTS == 23 || ADDITIONAL_PRINTS == 123
		fclose(file_recrange);
		#endif
		*/
	}
    /******************************** OUTPUT trees  ************************************/
    if((((*data).ifgamma == 1 || (*data).range_thetant)  && (*data).likelihood_line == 0) || (((*data).ifgammar == 1 || (*data).range_rnt) && (*data).likelihood_line == 0)) {
		#if ADDITIONAL_PRINTS == 2 || ADDITIONAL_PRINTS == 12 || ADDITIONAL_PRINTS == 23 || ADDITIONAL_PRINTS == 123
		fclose(file_ttotp);
		#endif
	}
    /****************** EXIT IF ERROR, AFTER PRINTING ADDITIONAL VALUES  ************************/
    if(y == 1) return 1; /*exit*/
    /****************************** PRINT STATS IN OUTPUT FILE **********************************/
	if((*data).neutral_tests) 
        if(print_neuttest(&data,file_output,file_out)) return 1; /*exit*/
    
	if((*data).neutral_tests == 0 && (*data).n_loci > 1 && (*data).likelihood_line == 0) {
		for(x=0;x<(*data).n_loci;x++) {
			fclose(file_outputm[x]);
		}
		free(file_outputm);
	}
	else fclose(file_output);
    /*************************************** FREE VECTORS ***************************************/
	free_matrix_sfixalltheta(&data);
    free_getpars_fix(&data,&inputp);
    free(inputp);
    free_inputdata(&data);
    free(data);
	
	#if TIMEDURATION
    /*time and duration*/
	if((*data).likelihood_line == 0) {
		time(&end);
		date = localtime(&end);
		strftime(s,80,"%c",date);
		fprintf(file_output,"\n\nDate of completion: %s\n",s);
		fprintf(file_output,"\nDuration of the process: %.0f seconds.",difftime(end,start));
		hour =  difftime(end,start)/3600.;
		min  = (difftime(end,start) - (double)hour*3600.)/60.;
		sec  =  difftime(end,start) - (double)hour*3600. - (double)min*60.;
		fprintf(file_output," (%dh:%dm:%ds)\n",hour,min,sec);
	}
	#endif

    printf("\n %s exited succesfully.\n\n",MLCOALSIM);
    return 0;
}
