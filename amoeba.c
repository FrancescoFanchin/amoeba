/*
 *	amoeba.c
 *
 *	A combination of Nelder-Mead method and simulated annealing  to find cost-function global minimum 
 *
 * to compile: gcc amoeba.c  funcomp.c final_results.c -o amoeba -lm
 *      to run:     ./amoeba
 *		to show fit: gnuplot f_fit.gp
 */      
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#define MAX 99999999999999

long double funcomp(long double param[],int endf);
double amoeba(int row,int col,long double simpts[][col]);
long double dist (int row,int col,long double simpts[][col],long double simval[]);
double ran2(long *idum);  
void final_results(int nloop,long double results[][8]);


main(){
	//define variables 

int i,j,k;
int row=8;
int col=7;
long double simpts1[8][7];
 int nloop=1; 
long double par[7]={1.01,
10000.16,
1.1,
0.000059,
0.000023798,
0.20079500646484246754,
100.09198731440190129061
};

long double dev1[7]={1.2,100,0.002,0.000002,0.00001,1.01,100.1};

long double E[nloop];
long double results[nloop][row];
long int idum=58634760;
//begin loop
k=0;


for (i=0;i<nloop;i++){
for (j=0;j<col+1;j++){
		results[i][j]=0.0;
			}}	



do{

printf("loop number %d\n",k);
	//randomize initial conditions

for(j=0;j<col;j++){
	simpts1[0][j]=(ran2(&idum))*par[j]+(ran2(&idum))*dev1[j];
	printf("row0 %Lf\n",simpts1[0][j]);
}
		idum=rand()/100;
printf("idum %ld\n",idum);
for(i=1;i<row;i++){
	simpts1[i][i-1]=par[i-1]+ran2(&idum)*dev1[i-1];
	for(j=0;j<col;j++){
		if(j!=i-1) simpts1[i][j]=ran2(&idum)*par[j] ;
	}
}
		idum=rand()/1000;
printf("idum %ld\n",idum);
//inf exception

if((k!=0)&&((E[k-1]>=220)||(E[k-1]<=-220))){

for(j=0;j<col;j++){
	simpts1[0][j]=par[j];
	printf("rowinf0 %Lf\n",simpts1[0][j]);
}
		idum=rand();
printf("idum %ld\n",idum);
for(i=1;i<row;i++){
	simpts1[i][i-1]=par[i-1]+dev1[i-1];
	for(j=0;j<col;j++){
		if(j!=i-1) simpts1[i][j]=par[j] ;
	}
}	
	for(i=0;i<row;i++){
		for(j=0;j<col;j++){ printf("simpts1[%d][%d] %Lf\n",i,j,simpts1[i][j]);}}
}


   //run amoeba

E[k]=	amoeba(row,col,simpts1);

	//compute final amoeba barycenter and store it in "results" with its own energy 

	long double sumbar[]={0.0,0.0,0.0,0.0,0.0,0.0,0.0};
	for(i=0;i<row;i++){
		for(j=0;j<col;j++){
		sumbar[j] += simpts1[i][j];
					}
				}
			long double baric[col];
	for(j=0;j<col;j++){
		baric[j]=sumbar[j]/row;
	}
		funcomp(baric,1);
results[k][0]=E[k];
for(j=0;j<col;j++){results[k][j+1]=baric[j];}		
		 k++;
	}while(k<nloop);	
		
		//find best fit and dump it on the screen

	long double min=0.0;
		
		for(i=0;i<nloop;i++){
		min=E[i];
			for(j=0;j<nloop;j++){
				if(E[j]<min){
					min=E[j];}
			}
		}
		int g=0;
		for(i=0;i<nloop;i++){if (min==E[i]) g=i;} 
		
printf("best fit:E= %Lf\n",results[g][0]);
	for (i=1;i<row;i++){printf("%.20Lf\n",results[g][i]);}		
	

		for (i=0;i<nloop;i++){printf("E%d = %Lf\n",i,results[i][0]);
			for (j=1;j<col+1;j++){
			printf("%.20Lf\n",results[i][j]);}}	
			
			final_results(nloop,results);
		return 0;
		}
		
			
 	
			//Nelder-Mead algorithm
	double amoeba(int row,int col,long double simpts[][col]){
long int idum=1230975; 
float T=10;
		int i,j,r;
		int u=0;
		int n=0;
		const int m=30;
		const float cool=0.99;
		long double temp[col];
		long double baric[col];
		long double sum[col];
		long double simval[row];
		long double reflpt[col];
		long double exppt[col];
		long double contrpt[col];
		long double contr_val=0.0000000;
		long double exp_val=0.00000000;
		long double refl_val=0.00000000;
		long double first_max=0.00000000;
		long double sec_max=0.0000000000;
		long double min=0.00000000000;
		int reflected=1;// maxpt is replaced by reflected pt
		int changemax2=1;//maxpt can change for reflection or expansion
		int reducted=1;
		int contracted=1;
		double epsilon=0.001;
	FILE* ET;
	ET=fopen("ET","w");
		//array inzialization
		for(j=0;j<col;j++){
			temp[j]=baric[j]=sum[j]=reflpt[j]=exppt[j]=contrpt[j]=0.0;}
			
			
		for(i=0;i<row;i++){
			simval[i]=0.0;
			}
			
		//begin cycle
		
		
	while(reducted==1){

		reducted=0;
		
	while(contracted==1){
		contracted=0;
		
		
	while(changemax2==1){
	changemax2=0;
	
	while(reflected==1){
	reflected=0;
	if (n%5==0){printf("T %f n %d %Lf\n",T,n,min);
		fprintf(ET,"%Lf %d %f\n",min,n,T);
		}
	for(i=0;i<row;i++){
	for(j=0;j<col;j++){
	temp[j]=simpts[i][j];
		}
		
		simval[i]=funcomp(temp,0)-T*(log(ran2(&idum)));
		if (isnan(simval[i])||((isinf(simval[i]))==1)||(simval[i]>MAX)) {
			for(j=0;j<col;j++){simpts[i][j]=ran2(&idum)*simpts[0][j];}
		}
		for(j=0;j<col;j++){
	temp[j]=simpts[i][j];
		}
		simval[i]=funcomp(temp,0)-T*(log(ran2(&idum)));
		if (isnan(simval[i])||((isinf(simval[i]))==1)||(simval[i]>MAX)) {
			simval[i]=MAX;}
		}
	


	//find max value
	
	for(i=0;i<row;i++){
		first_max=simval[i];
		for(j=0;j<row;j++){
			if(simval[j]>first_max){
				first_max=simval[j];}
				}
			}
for(i=0;i<row;i++){
		if(first_max==simval[i]){
			u=i;}
			}


			//find second max value
	
	for(i=0;i<row;i++){
		if(i!=u) {sec_max=simval[i];
		for(j=0;j<row;j++){
			if((simval[j]>sec_max)&&(j!=u)){
				sec_max=simval[j];}
			}
			}
		}
		
		
		//find min
	for(i=0;i<row;i++){
		min=simval[i];
			for(j=0;j<row;j++){
				if(simval[j]<min){
					min=simval[j];}
			}
		}
		int g=0;
		for(i=0;i<row;i++){if (min==simval[i]) g=i;} 
			
		//calculating barycenter ...
	
	for(j=0;j<col;j++){sum[j]=0.0;}
		

				for(j=0;j<col;j++){
					for(i=0;i<row;i++){
						if (i!=u){
							sum[j] +=simpts[i][j];
					}
				}
			}
		
	

	for(j=0;j<col;j++){
		baric[j]=sum[j]/(row-1);
		}

	//reflex
			
	for (j=0;j<col;j++){
		reflpt[j]=baric[j]+1*(baric[j]-simpts[u][j]);
		if (reflpt[j]<0)reflpt[j]=0;
	}


	//compute function value in reflexes point
	
	refl_val=funcomp(reflpt,0)+T*(log(ran2(&idum)));
	if (isnan(refl_val)||((isinf(refl_val))==1)||(refl_val>MAX)) {
			refl_val=MAX;
			
		}

	if((min<=refl_val)&&(refl_val<sec_max)){
		for(j=0;j<col;j++){
			simpts[u][j]=reflpt[j];
		}
		
		reflected=1;
		reducted=1;
		contracted=1;
		changemax2=1;
	
		if(first_max-min<epsilon) {printf("endrefl");
			fclose(ET);
			return min;
			} 
	}
n++;
if(n%m==0)T*=cool;
	}
		
		
// eventually "expand"
	if (refl_val<min){
		for (j=0;j<col;j++){
				exppt[j]=(3*baric[j])-(2*simpts[u][j]);
			if (exppt[j]<0)exppt[j]=0;
			}
		
			exp_val=funcomp(exppt,0)+T*(log(ran2(&idum)));
			if (isnan(exp_val)||((isinf(exp_val))==1)||(exp_val>MAX)) {
			exp_val=MAX;
			
		}
	if(exp_val<refl_val){
		for (j=0;j<col;j++){
			simpts[u][j]=exppt[j];
			
		}
		reducted=1;
		contracted=1;
		changemax2=1;
		reflected=1;

		if(first_max-min<epsilon){printf("endchangemax2A");
			fclose(ET);
			return min;
			} 
	}

			
		if(exp_val>=refl_val){
		for (j=0;j<col;j++){
			simpts[u][j]=reflpt[j];
			
		}
		reducted=1;
		contracted=1;
		changemax2=1;
		reflected=1;
	

		if(first_max-min<epsilon){printf("endchangemax2B\n");
			fclose(ET);
			return min;
			} 
	}
	}	

n++;
if(n%m==0)T*=cool;
}		
		
		
		//compute contrpt e contr_val
		
	for (j=0;j<col;j++){
		contrpt[j]=baric[j]-0.5*(baric[j]-simpts[u][j]);
		if (contrpt[j]<0)contrpt[j]=0;
	}		
			
		
		contr_val=funcomp(contrpt,0)+T*(log(ran2(&idum)));
		if (isnan(contr_val)||((isinf(contr_val))==1)||(contr_val>MAX)) {
			contr_val=MAX;
			
		}	

			//eventual contraction
	if(contr_val<first_max){
			for(j=0;j<col;j++){
			simpts[u][j]=contrpt[j];
		}

		reducted=1;
		contracted=1;
		changemax2=1;
		reflected=1;
		if(first_max-min<epsilon){printf("endcontr %Lf max %Lf\n",min,first_max);
			fclose(ET);
			return min;
			} 
	}
	


n++;
if(n%m==0)T*=cool;
}	
			//find min row
	
	for(i=0;i<col;i++){
		if(min==simval[i]) u=i;}
		
		
		
		//reduction
		for(i=0;i<row;i++){
			for(j=0;j<col;j++){	
				if(i!=u) {
					simpts[i][j]= simpts[u][j] +0.5*(simpts[i][j]-simpts[u][j]);
				}
			}	
		}
			
			reducted=1;
			contracted=1;
		changemax2=1;
		reflected=1;

			if(first_max-min<epsilon) {
				printf("endred min %Lf max %Lf\n",min,first_max);
				fclose(ET);
				return min;}
			
n++;
if(n%m==0)T*=cool;
		}
			fclose(ET);
		return min;
		}
		
			
long double dist(int row,int col,long double simpts[][col],long double simval[]){
	int i,j;
	long double first_max=0.0;
	long double min=0.0;
	long double distance=0.0;
	
	
	
	
	
	//find max
	
	
	for(i=0;i<row;i++){
		first_max=simval[i];
		for(j=0;j<row;j++){
			if(simval[j]>first_max){
				first_max=simval[j];}
				}
			}
			int u=0;
for(i=0;i<row;i++){
		if(first_max==simval[i]){
			u=i;}
			}
		
//find min
	for(i=0;i<row;i++){
		min=simval[i];
			for(j=0;j<row;j++){
				if(simval[j]<min){
					min=simval[j];}
			}
		}
		int g=0;
		for(i=0;i<row;i++){if (min==simval[i]) g=i;} 

for (j=0;j<col;j++){
	distance += (((simpts[u][j])-(simpts[g][j]))*((simpts[u][j])-(simpts[g][j])));
}
distance=sqrt(distance);
printf("distance %Lf\n",distance);
return distance;

}



/* Random number generator */

#define IM1 2147483563
#define IM2 2147483399
#define AM (1.0/IM1)
#define IMM1 (IM1-1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 3791
#define NTAB 32
#define NDIV (1+IMM1/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)

double ran2(long *idum)
{
        int j;
        long k;
        static long idum2=123456789;
        static long iy=0;
        static long iv[NTAB];
        double temp;

        if (*idum <= 0) {
                if (-(*idum) < 1) *idum=1;
                else *idum = -(*idum);
                idum2=(*idum);
                for (j=NTAB+7;j>=0;j--) {
                        k=(*idum)/IQ1;
                        *idum=IA1*(*idum-k*IQ1)-k*IR1;
                        if (*idum < 0) *idum += IM1;
                        if (j < NTAB) iv[j] = *idum;
                }
                iy=iv[0];
        }
        k=(*idum)/IQ1;
        *idum=IA1*(*idum-k*IQ1)-k*IR1;
        if (*idum < 0) *idum += IM1;
        k=idum2/IQ2;
        idum2=IA2*(idum2-k*IQ2)-k*IR2;
        if (idum2 < 0) idum2 += IM2;
        j=iy/NDIV;
        iy=iv[j]-idum2;
        iv[j] = *idum;
        if (iy < 1) iy += IMM1;
        if ((temp=AM*iy) > RNMX) return RNMX;
        else return temp;
}
#undef IM1
#undef IM2
#undef AM
#undef IMM1
#undef IA1
#undef IA2
#undef IQ1
#undef IQ2
#undef IR1
#undef IR2
#undef NTAB
#undef NDIV
#undef EPS
#undef RNMX
