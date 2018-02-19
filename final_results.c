/* 
 * final_results.c
 * 
 * compute final results: E_PDF,E_CDF, EX2
 *
 * to compile: gcc amoeba.c  funcomp.c final_results.c -o amoeba -lm
 *      to run:     ./amoeba
 *		to show fit: gnuplot f_fit.gp
*/


#include <math.h>
#include <stdio.h>
#define	MAX	100000


void final_results(int nloop,long double results[][8]){
	
int i,j,k;
double c[MAX],t[MAX],n[MAX],p[MAX],f[MAX];
double t_end;             /* dt = time step; t_end = final time */
int n_steps,dt;                 /* n_steps = number of iterations     */
double kappa, lambda,mu,N0,beta,alpha,c_star,delta,A; //parameters
long double Epdf,Ecdf,EX2;	
FILE *fp;

double fexp[121]={0.000000,0.000000,0.000000,0.005000,0.000000,0.000000,0.010000,0.015000,0.020000,0.034000,0.068000,0.078000,0.127000,0.210000,0.264000,0.381000,0.537000,0.649000,0.913000,1.094000,1.382000,1.470000,1.611000,1.748000,2.173000,2.241000,2.437000,2.612000,2.715000,2.930000,3.027000,3.179000,3.208000,3.311000,3.359000,3.306000,3.262000,3.013000,3.057000,2.988000,2.725000,2.695000,2.481000,2.412000,2.246000,2.051000,1.992000,1.768000,1.631000,1.582000,1.533000,1.499000,1.328000,1.108000,1.162000,1.123000,1.021000,0.903000,0.830000,0.752000,0.752000,0.693000,0.723000,0.645000,0.469000,0.581000,0.498000,0.464000,0.508000,0.439000,0.361000,0.371000,0.381000,0.356000,0.288000,0.342000,0.254000,0.303000,0.254000,0.225000,0.249000,0.190000,0.176000,0.171000,0.156000,0.151000,0.151000,0.132000,0.103000,0.059000,0.073000,0.088000,0.063000,0.049000,0.039000,0.054000,0.049000,0.029000,0.039000,0.024000,0.020000,0.015000,0.010000,0.010000,0.010000,0.020000,0.010000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000};

double pexp[121]={0,0,0,0.005,0.005,0.005,0.015,0.03,0.05,0.084,0.152,0.23,0.357,0.567,0.831,1.212,1.749,2.398,3.311,4.405,5.787,7.257,8.868,10.616,12.789,15.03,17.467,20.079,22.794,25.724,28.751,31.93,35.138,38.449,41.808,45.114,48.376,51.389,54.446,57.434,60.159,62.854,65.335,67.747,69.993,72.044,74.036,75.804,77.435,79.017,80.55,82.049,83.377,84.485,85.647,86.77,87.791,88.694,89.524,90.276,91.028,91.721,92.444,93.089,93.558,94.139,94.637,95.101,95.609,96.048,96.409,96.78,97.161,97.517,97.805,98.147,98.401,98.704,98.958,99.183,99.432,99.622,99.798,99.969,100.125,100.276,100.427,100.559,100.662,100.721,100.794,100.882,100.945,100.994,101.033,101.087,101.136,101.165,101.204,101.228,101.248,101.263,101.273,101.283,101.293,101.313,101.323,101.323,101.323,101.323,101.323,101.323,101.323,101.323,101.323,101.323,101.323,101.323,101.323,101.323,101.323};

for(k=0;k<nloop;k++){
		

	kappa=results[k][1];

	lambda=results[k][2];
	mu=results[k][3];
	alpha=results[k][4];
	beta=results[k][5];

	delta=results[k][6]; 
	A=results[k][7];
	
	fp = fopen("result.dat","w");

	dt=1;
	t[0] = 0.0;
	t_end = 72000;
	    n_steps = (t_end+1.0E-6)/dt;

	c[0]=p[0]=0.0; 
	n[0]=pow(10,kappa*exp(-exp(1+((exp(1)*mu/kappa)*(lambda-t[1]))))); ;

	
	i=0;
	t[1] = t[0] + dt;   //update time
	  n[1]=pow(10,kappa*exp(-exp(1+((exp(1)*mu/kappa)*(lambda-t[1])))));  //number(concentration) of bacteria
	  
	  c[1] = c[0] - dt*(beta*c[0] -alpha*n[0]);                           //number of autoinducers
	  	  

	  p[1]=A*n[1]/(1+pow((1/c[1]),delta)); 
	  	
	  		  f[0]= ((p[1]-p[0])/dt);
	  		
	
    fprintf(fp,"%lf\t%lf\t%lf\t%lf\t%lf\n",t[0],n[0],c[0],p[0],f[0]);
	  		  
	  		  for(i=1;i<n_steps;i++){

	  t[i+1] = t[i] + dt;   //update time
	  n[i+1]=pow(10,kappa*exp(-exp(1+((exp(1)*mu/kappa)*(lambda-t[i+1])))));  //number(concentration) of bacteria
	  c[i+1] = c[i] - dt*(beta*c[i] -alpha*n[i]);                           //number of autoinducers
	  p[i+1]=A*n[i+1]/(1+pow((1/c[i+1]),delta));  
	  f[i]=  (0.5*(p[i+1]-p[i-1])/dt);
	  
	  fprintf(fp,"%f\t%f\t%f\t%f\t%f\n",t[i],n[i],c[i],p[i],f[i]);

		}
		
     fprintf(fp,"%f\t%f\t%f\t%f\t%f\n",t[n_steps],n[n_steps],c[n_steps],p[n_steps],f[n_steps]);   /* final time output    */
     fclose(fp);
	
	long double S[MAX];
	
	
	for(i=0;i<n_steps;i++){S[i]=0.0;
	for(j=0;j<i+1;j++){ 
		S[i]+=(f[j]);}

	}
	
	Epdf=0.0;
	 for (i=0;i<n_steps;i=i+(600/dt)){
		
		
		
	    Epdf += ((f[i]-fexp[dt*i/600])*(f[i]-fexp[dt*i/600]));
	
	}
	
	Ecdf=0.0;
	
	for (i=0;i<n_steps;i=i+(600/dt)){
		
		

	    Ecdf += (((S[i]/600)-pexp[dt*i/600])*((S[i]/600)-pexp[dt*i/600]));
	
	}
	
	EX2=0.0;
	
	for (i=0;i<n_steps;i=i+(600/dt)){
		
		

	     EX2 += (((f[i]-fexp[dt*i/600])*(f[i]-fexp[dt*i/600]))/f[i]);
	
	}
	printf("result %d E= %Lf\n",k,results[k][0]);
	for(j=1;j<8;j++){printf("%.15Lf\n",results[k][j]);}
	printf("E_PDF= %Lf\nE_CDF= %Lf\nE_X2= %Lf\n",Epdf,Ecdf,EX2);   
	
	//plot fit with gnuplot
	system("gnuplot f_fit.gp");
	
}	
	return;
	
}

