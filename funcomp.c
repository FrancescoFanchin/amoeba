/*
 *	funcomp.c
 *
 *	compute cost-function
 *	parameters have been taken by amoeba.c
 *      Equations are solved by means of Eulero method in iterate()
 * 
 * to compile: gcc amoeba.c  funcomp.c final_results.c -o amoeba -lm
 *      to run:     ./amoeba
 *		to show fit: gnuplot f_fit.gp
 *      
 */

#include <math.h>
#include <stdio.h>
#define	MAX	100000



long double funcomp(long double param[],int endf){
	double num=0.0;
	double E=0.00;
	double c[MAX],t[MAX],n[MAX],p[MAX],f[MAX];
int dt;             /* dt = time step; t_end = final time */
int n_steps;
	int t_end = 72000;
                 /* n_steps = number of iterations     */
double kappa, lambda,mu,beta,alpha,delta,A; //parameters
FILE *ofile,*end_param;
	
	
	double fexp[121]={0.000000,0.000000,0.000000,0.005000,0.000000,0.000000,0.010000,0.015000,0.020000,0.034000,0.068000,0.078000,0.127000,0.210000,0.264000,0.381000,0.537000,0.649000,0.913000,1.094000,1.382000,1.470000,1.611000,1.748000,2.173000,2.241000,2.437000,2.612000,2.715000,2.930000,3.027000,3.179000,3.208000,3.311000,3.359000,3.306000,3.262000,3.013000,3.057000,2.988000,2.725000,2.695000,2.481000,2.412000,2.246000,2.051000,1.992000,1.768000,1.631000,1.582000,1.533000,1.499000,1.328000,1.108000,1.162000,1.123000,1.021000,0.903000,0.830000,0.752000,0.752000,0.693000,0.723000,0.645000,0.469000,0.581000,0.498000,0.464000,0.508000,0.439000,0.361000,0.371000,0.381000,0.356000,0.288000,0.342000,0.254000,0.303000,0.254000,0.225000,0.249000,0.190000,0.176000,0.171000,0.156000,0.151000,0.151000,0.132000,0.103000,0.059000,0.073000,0.088000,0.063000,0.049000,0.039000,0.054000,0.049000,0.029000,0.039000,0.024000,0.020000,0.015000,0.010000,0.010000,0.010000,0.020000,0.010000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000};

 
    n_steps = (t_end+1.0E-6)/dt;
	dt=600;
	
	t[0] = 0.0;
	kappa=param[0];
	lambda=param[1];
	mu=param[2];
	c[0]=p[0]=0.0;
	n[0]=pow(10,kappa*exp(-exp(1+((exp(1)*mu/kappa)*(lambda-t[0])))));
	alpha=param[3];
	beta=param[4];
	delta=param[5];
	A=param[6];
	         /* parte intera! */
	
	if (n_steps>MAX-1) n_steps=MAX-1;    /* avoiding too big arrays */

	int i,j;
	
	

	i=0;
	t[1] = t[0] + dt;   //update time
	  n[1]=pow(10,kappa*exp(-exp(1+((exp(1)*mu/kappa)*(lambda-t[1])))));  //number(concentration) of bacteria

	  c[1] = c[0] - dt*(beta*c[0] -alpha*n[0]);                           //number of autoinducers
		  

	  p[1]=A*n[1]/(1+pow((1/c[1]),delta)); 
	    
	  f[0]= ((p[1]-p[0])/dt);

	  
	  if(isnan(f[0])||((isinf(f[0]))==1)){ 
	  E=sqrt(-1);
	  return E;}
	
	for(i=1;i<n_steps;i=i+(600/dt)){

	  t[i+1] = t[i] + dt;   //update time

	  n[i+1]=pow(10,kappa*exp(-exp(1+((exp(1)*mu/kappa)*(lambda-t[i+1])))));  //number(concentration) of bacteria

	  c[i+1] = c[i] - dt*(beta*c[i] -alpha*n[i]);                           //number of autoinducers

	  p[i+1]=A*n[i+1]/(1+pow((1/c[i+1]),delta));                         //number of activated cells
	  
	  f[i]=  (0.5*(p[i+1]-p[i-1])/dt);                                    //radiant flux
	  if(isnan(f[i+1])||((isinf(f[i+1]))==1)){ 
	  E=sqrt(-1);
	  return E;}
	   	 
                                        
	}
    
   
    for (i=0;i<n_steps;i=i+(600/dt)){
		
		
	    E += ((f[i]-fexp[dt*i/600])*(f[i]-fexp[dt*i/600]));
	}
	

	
	if(endf==1) {
		
		
	
		
		
	end_param=fopen("end_param","w");
	for(i=0;i<7;i++){
	fprintf(end_param,"%.12Lf\n",param[i]);
	}
	
	return E;
}
long long loge=log(E);
return E;
}  

	

	

