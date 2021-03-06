#include <stdio.h> 
#include <stdlib.h>
#include <time.h>
#include <math.h> // for the exp function 
#include <string.h> // for concatenate chars

/*Function declarations-- */
double rand_double(); /* 0.0 to 1.0, unif. dist. */
int mod (int a, int b);
void state0(int N1, int N2, int* p);/* random initial state*/
double calc_Energy(int *p,int *q,int N1,int N2,double J, double J2);
double calc_mag(int *p,int *q,int N1,int N2);
void montecarlo_step(int N1, int N2 ,double J,double J2, double T,int *p,int *q,double* e,double*m);
void plot_termalizacion(int N1,int N2,double J,double J2,double T,int* p,int* q,double* e,double* m,int tau);
void export_data(int* p, int* q, int N1, int N2, double T,FILE* d);
void iniciar(int*p, int*q,int N1,int N2,double J, double J2,double*e,double* m);
void AnnealingT(double* e, double* m, int Nsites, int Neq,
double T0, int N1, int N2,int* p, int* q, double J,double J2,FILE* d,
double*ep,double*e2p,double*mp,double*m2p,int exportflag);
void exportThermo(double T0, double NT, int Nsites, int Nprom,double* ep,
double* e2p, double*mp,double*m2p,FILE*f);


/* -------------------------- Main -------------------------- */
int main() {
clock_t start, end;
double cpu_time_used;

int N2=15;
int N1=2*N2;
int Nsites=2*N1*N2;
int Nprom=1000;
int Neq=1000;// number of montecarlosteps / sites that i wait in each temperature to thermalize
//int tau=10000;// time in units of  montecarlo steps / site 
int lattA[N1][N2]; //uninitialized array, sublattice A
int lattB[N1][N2]; //uninitialized array sublattice B
for (int i=0;i<N1;i++){//set array to zero
	for (int j=0;j<N2;j++){
	lattA[i][j]=0;
	lattB[i][j]=0;
	}
}
double mag=0.; //total magnetization 
double Energy=0.; //total energy
int *p =&lattA[0][0]; /* p points to lattA first element */
int *q =&lattB[0][0];/* q points to lattB first element */
double *m=&mag;
double *e=&Energy;

double J=-1.; // NN coupling, positive means ferromagnetic .
double T0=4.53;// initial temperature
//double dT=-0.02265; // temperature step
int NT=200; // number of temperatures
double eprom[NT]; //thermodinamic variables
double e2prom[NT];
double mprom[NT];
double m2prom[NT];
double* ep=eprom;
double* e2p=e2prom;
double* mp=mprom;
double* m2p=m2prom;
//export check
/*
for(int i=0;i<N1;i++){ 
	int j=mod(-2*i,N2); 
	for(int k=0;k<N2;k++){ // repito N2 veces 
		printf("%d ",mod(j,N2) );
		j=mod(j+1,N2);
	}
	printf(" ");
}
*/ 

/*
//thermalization check
double J2=0.;
double T=.5;
iniciar(p,q,N1,N2,J,J2,e,m);
plot_termalizacion( N1, N2, J, J2, T,p,q,e,m,tau);
*/



double J2list[ ]={-1.0}; //positive means ferromagnetic 
double* J2p=J2list;
int size = sizeof J2list / sizeof J2list[0];
printf("J2list has %d elements\n",size);

for (int k=0;k<size;k++){
	printf("J2=%.2f\n",*(J2p+k));
	start = clock();
	for ( int j=0; j< NT;j++){ //initialize counters
		eprom[j]=0;
		e2prom[j]=0;
		mprom[j]=0;
		m2prom[j]=0;
	} 
	char filenameA [100];
	sprintf (filenameA,"/home/s/Desktop/arrays_%d_%.2f.txt",N1,*(J2p+k));
	FILE *d = fopen(filenameA, "w"); // open file to export arrays

	char filenameT [100];
	sprintf (filenameT,"dataAF/termodinamica_%d_%.2f.txt",N1,*(J2p+k));
	FILE *f = fopen(filenameT, "w"); // open file to write thermodynamic variables

	int exportflag=1;
	for(int ii=1;ii<=Nprom;ii++){// Nprom independent annealings
		srand((long)time(NULL)); /* initialize RNG rand() */
		iniciar(p,q,N1,N2,J,*(J2p+k),e,m);
		AnnealingT(e,m,Nsites,Neq,T0,N1,N2,p,q,J,*(J2p+k),d,ep,e2p,mp,m2p,exportflag);
		printf("%d, m=%f, e=%f\n",ii,*m/(double)Nsites,*e/(double)Nsites);
	}
	//export mean thermodinamic quantities 
	exportThermo(T0,NT,Nsites,Nprom,ep,e2p, mp,m2p,f);
	printf("\n");
	fclose(d); //close configuration file 
	fclose(f);//close thermodynamics file 
	end = clock();
	cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
	printf("this task took %.1f seconds\n", cpu_time_used); 
}


return 0;
} /* ---------------------- End Main ---------------------- */

void exportThermo(double T0, double NT, int Nsites, int Nprom,double* ep,
double* e2p, double*mp,double*m2p,FILE*f){
int k=0;
for (double T=T0;T>0.0226;T+=-0.02265){
	*(ep+k)=*(ep+k)/(float)Nprom/Nsites; //mean value + normalization 
	*(e2p+k)=*(e2p+k)/(float)Nprom/Nsites/Nsites; 
	*(mp+k)=*(mp+k)/(float)Nprom/Nsites;  
	*(m2p+k)=*(m2p+k)/(float)Nprom/Nsites/Nsites;
	fprintf(f,"%f %f %f %f %f\n",T,*(mp+k),*(m2p+k),*(ep+k),*(e2p+k) ) ;
	k+=1;
	}
}



void iniciar(int*p, int*q,int N1,int N2,double J, double J2,double*e ,double*m){
	state0(N1,N2,p); /* estado inicial subred A*/ 
	state0(N1,N2,q);/* estado inicial subred B*/ 
	//printf ("Energy0/site= %f\n",Energy/(double)(N1*N2)/2.);
	*e=calc_Energy(p,q,N1,N2,J,J2); // initial energy
	*m=calc_mag(p,q,N1, N2); // initial total mag 
}
//----------------------- fin iniciar-



void export_data(int* p, int* q, int N1, int N2, double T,FILE* d){

/* 
if(T<Tc ){fprintf(d,"%d ",0);}//label
else{fprintf(d,"%d ",1);};
*/
fprintf(d,"%f ",T);//Temperature
for(int i=0;i<N1;i++){ 
	int j=mod(-2*i,N2); 
	for(int k=0;k<N2;k++){ // N2 times
		fprintf(d,"%d %d ",*(p+mod(j,N2)+i*N2),*(q+ mod(j,N2) +i*N2) );
		j=mod(j+1,N2);
	}
}
fprintf(d,"\n");

}//end export_data ---------------------------------------------------


void AnnealingT(double* e, double* m, int Nsites, int Neq,
double T0, int N1, int N2,int* p, int* q, double J,double J2,FILE* d,
double*ep,double*e2p,double*mp,double*m2p,int exportflag){
int ii=0;
	for (double T=T0;T>0.0226;T+=-0.02265){
		// equilibrium estimation for each temperature:
		for (int t=0;t<Nsites*Neq;t++){montecarlo_step(N1,N2,J,J2,T,p,q,e,m);}
		//printf("T=%f, e=%f, m=%f\n",T,*e/Nsites,*m/Nsites);
		if(abs(*m/Nsites)<1.1){
			if (exportflag==1){export_data(p,q,N1,N2,T,d);};
			*(ep+ii)=*(ep+ii)+*e;
			*(e2p+ii)=*(e2p+ii)+pow(*e,2);
			*(mp+ii)=*(mp+ii)+abs(*m);  
			*(m2p+ii)=*(m2p+ii)+pow(*m,2);
			}
		else{
			printf("ERROR");
			}
		ii+=1;
	}
}

//--------------------------------------------

//--------------------------------------------



void plot_termalizacion(int N1,int N2,double J,double J2,double T,int* p,int* q,double* e,double* m, int tau){
FILE *f = fopen("termalizacion.txt", "w"); // open file to write magnetizations vs time.
if (f == NULL){
    printf("Error opening file!\n");
    exit(1);
}
for (int time=0;time<=tau;time++){
	fprintf(f,"%d %f %f\n",time,*m/(double)(N1*N2)/2.,*e/(double)(N1*N2)/2.);
	for (int j1=1;j1<=2*N1*N2;j1++){montecarlo_step(N1,N2,J,J2,T,p,q,e,m);}
}

fclose(f);
}//end plot_termalizacion

//--------------------------------------------------------------------------------

//montecarlo step

void montecarlo_step(int N1, int N2, double J,double J2,double T,int *p,int *q,double* e,double* m){
int rand_i=rand() % N1;//random cell
int rand_j=rand() % N2;
if(abs(rand_i)>=N1 || abs(rand_j)>=N2){printf("ERROR-RAND:i=%d, j=%d \n",rand_i,rand_j);};
double Delta_E=0.; //Delta_E came from eq 3.10  Barkema
double beta=1./T;
int sumNN=0;
int sumNNN=0;
double epsilonE=1.E-5;


double r1=rand_double();//random sublattice
//printf("%f\n",r1);
if (r1<0.5){ // spin Bij
	sumNN+= *(p +mod(rand_i,N1)*N2 + mod(rand_j,N2) );
	sumNN+= *(p + mod(rand_j+1,N2)+mod(rand_i,N1)*N2);
	sumNN+= *(p + mod(rand_j-1,N2)+mod(rand_i+1,N1)*N2 );

	sumNNN+=*(q +mod(rand_i-1,N1)*N2 + mod(rand_j+2,N2) );
	sumNNN+=*(q +mod(rand_i,N1)*N2   + mod(rand_j+1,N2) );
	sumNNN+=*(q +mod(rand_i+1,N1)*N2 + mod(rand_j-1,N2) );
	sumNNN+=*(q +mod(rand_i+1,N1)*N2 + mod(rand_j-2,N2) );
	sumNNN+=*(q +mod(rand_i,N1)*N2   + mod(rand_j-1,N2) );
	sumNNN+=*(q +mod(rand_i-1,N1)*N2 + mod(rand_j+1,N2) );

	Delta_E+=2.*  *(q + rand_j+rand_i*N2)  * (J* sumNN+J2* sumNNN);
	//printf("DeltaE1= %f\n",Delta_E);
	if(Delta_E<=epsilonE){ 
		*(q + rand_j+rand_i*N2) = - *(q + rand_j+rand_i*N2);
		*e=*e + Delta_E;
		*m=*m -2.* *(q + rand_j+rand_i*N2) ; // staggered !!!!!!
	}// if the energy decreases or stays the same, flip.
	else{
		double r2=rand_double();
		if(r2<exp(-beta*Delta_E) ){
			*(q + rand_j+rand_i*N2) = -*(q + rand_j+rand_i*N2);
			*e=*e +Delta_E;
			*m=*m -2.*  *(q + rand_j+rand_i*N2) ; // staggered !!!!!!
		}
	}
}
else{ // i.e., r>0.5 // spin Aij
	sumNN+= *(q + mod(rand_j,N2)+mod(rand_i,N1)*N2);
	sumNN+= *(q + mod(rand_j+1,N2)+mod(rand_i-1,N1)*N2);
	sumNN+= *(q + mod(rand_j-1,N2)+mod(rand_i,N1)*N2 );

	sumNNN+=*(p +mod(rand_i-1,N1)*N2 + mod(rand_j+2,N2) );
	sumNNN+=*(p +mod(rand_i,N1)*N2   + mod(rand_j+1,N2) );
	sumNNN+=*(p +mod(rand_i+1,N1)*N2 + mod(rand_j-1,N2) );
	sumNNN+=*(p +mod(rand_i+1,N1)*N2 + mod(rand_j-2,N2) );
	sumNNN+=*(p +mod(rand_i,N1)*N2   + mod(rand_j-1,N2) );
	sumNNN+=*(p +mod(rand_i-1,N1)*N2 + mod(rand_j+1,N2) );


	//Delta_E+=2.*J* *(p + rand_j+rand_i*N2) *sumNN;
	Delta_E+=2.* *(p + rand_j+rand_i*N2)  * (J* sumNN+J2 * sumNNN);
	//printf("DeltaE2= %f\n",Delta_E);
	if(Delta_E<=epsilonE){ 
		*(p + rand_j+rand_i*N2) = -*(p + rand_j+rand_i*N2);
		*e=*e + Delta_E;
		*m=*m+ 2.*  *(p + rand_j+rand_i*N2);
}// if the energy decreases or stays the same, flip.
else{
	double r2=rand_double();
	if(r2<exp(-beta*Delta_E) ){
	*(p + rand_j+rand_i*N2) = -*(p + rand_j+rand_i*N2);
	*e=*e + Delta_E;
	*m=*m + 2.*  *(p + rand_j+rand_i*N2);
	}
}
}
//printf("%f\n",exp(-beta*Delta_E));
//printf("Delta_E= %f\n",Delta_E);
//printf("E= %f\n",*e);
return ;
}//end montecarlo_step




//calculates the energy of the spin configuration: 
double calc_Energy(int *p,int *q,int N1,int N2,double J,double J2){ 
//printf("  \n");
double EnergyP=0.;
for (int i=0 ;i<N1 ;i++ ){  
	for (int j=0 ;j<N2 ;j++ ){
		EnergyP+=(-1)*J*( *(q + j+i*N2  ) * *(p + mod(j,N2)+mod(i,N1)*N2)  );
		EnergyP+=(-1)*J*( *(q + j+i*N2  ) * *(p + mod(j+1,N2)+mod(i,N1)*N2) );
		EnergyP+=(-1)*J*( *(q + j+i*N2  ) * *(p +mod(j-1,N2)+mod(i+1,N1)*N2) );


		EnergyP+=(-1)* J2 *( *(p +j+i*N2  ) * *(p +mod(i-1,N1)*N2 + mod(j+2,N2)  ));
		EnergyP+=(-1)* J2 *( *(p +j+i*N2  ) * *(p +mod(i,N1)*N2 + mod(j+1,N2)  ));
		EnergyP+=(-1)* J2 *( *(p +j+i*N2  ) * *(p +mod(i+1,N1)*N2 + mod(j-1,N2)  ));

		EnergyP+=(-1)* J2 *( *(q +j+i*N2  ) * *(q +mod(i-1,N1)*N2 + mod(j+2,N2)  ));
		EnergyP+=(-1)* J2 *( *(q +j+i*N2  ) * *(q +mod(i,N1)*N2 + mod(j+1,N2)  ));
		EnergyP+=(-1)* J2 *( *(q +j+i*N2  ) * *(q +mod(i+1,N1)*N2 + mod(j-1,N2)  ));
	}
}
return EnergyP;
}//end of  calc_Energy()

//calculates the magnetization of the spin configuration:
double calc_mag(int *p,int *q,int N1,int N2){ 
double magg=0.;
for (int i=0 ;i<N1 ;i++ ){ 
for (int j=0 ;j<N2 ;j++ ){
magg+=*(p + j+i*N2);//sublattice A
magg+=-*(q + j+i*N2);//sublattice B       staggered! 
}
}
return magg;
}//end of calc_Energy()




/* random initial state*/
void state0(int N1, int N2, int *pp){

for ( int j=0; j< N1*N2;j++){
double r=rand_double();
if (r<0.5) {
*(pp+j)=-1; }
else {
*(pp+j)=1;} 

//printf("%i",*(p+j));
}
} /* end state0 */




/* */ 
double rand_double() {
   return rand()/(double)RAND_MAX;
/*fuente: http://www.cs.utsa.edu/~wagner/CS2073/random/random.html */
}



int mod (int a, int b) 
{
   if (b < 0) //you can check for b == 0 separately and do what you want
     return mod(a, -b);   
   int ret = a % b;
   if (ret < 0)
     ret+=b;
   return ret;
}
/* */


