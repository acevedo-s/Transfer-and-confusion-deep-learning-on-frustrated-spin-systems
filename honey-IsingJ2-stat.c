#include <stdio.h> // para imprimir en consola
#include <stdlib.h>
#include <time.h>
#include <math.h> // for the exp function 
#include <string.h> // for concatenate chars

/*Function declarations-- */
double rand_double(); /* 0.0 to 1.0, unif. dist. */
double rand_double2(double x, double y); /* a to b */
int rand_int(int a, int b); /* a to b, inclusive */
int mod (int a, int b);
void state0(int N1, int N2, int* p);/* esta funcion inicializa el sistema aleatoriamente*/
double calc_Energy(int *p,int *q,int N1,int N2,double J, double J2);
double calc_mag(int *p,int *q,int N1,int N2);
void montecarlo_step(int N1, int N2 ,double J,double J2, double T,int *p,int *q,double* e,double*m);
//la funcion montecarlo step realiza *un* intento de flip. 
void plot_termalizacion(int N1,int N2,double J,double J2,double T,int* p,int* q,double* e,double* m,int tau);
void Annealing(double* e, double* m, int Nsites, int Neq,
double T0, int N1, int N2,int* p, int* q, double J,double J2,FILE* d,double Tc);
void export_data(int* p, int* q, int N1, int N2, double T, double Tc,FILE* d);
void iniciar(int*p, int*q,int N1,int N2,double J, double J2,double*e,double* m);
void AnnealingT(double* e, double* m, int Nsites, int Neq,
double T0, int N1, int N2,int* p, int* q, double J,double J2,FILE* d,double Tc,
double*ep,double*e2p,double*mp,double*m2p);

/* -------------------------- Main -------------------------- */
int main() {

int N2=5;
int N1=2*N2;
int Nprom=400;
int Neq=250;// numero de montecarlosteps / sites para que termalice el sistema.
//int tau=10000;// tiempo medido en montecarlo steps / site 
int lattA[N1][N2]; //defino arrays sin iniciar 
int lattB[N1][N2]; //sin iniciar 
for (int i=0;i<N1;i++){//inicializo lattice a cero.
	for (int j=0;j<N2;j++){
		lattA[i][j]=0;
		lattB[i][j]=0;
	}
}
double mag=0.; //magnetización total del sistema 
double Energy=0.; //Energia total del sistema  

int *p =&lattA[0][0]; /* p apunta al primer elemento de lattA */
int *q =&lattB[0][0];/* q apunta al primer elemento de lattB */
double *m=&mag;
double *e=&Energy;


double J=-1.; // NN coupling, positivo es ferro.
double T0=4.53;// temperatura inicial
//double dT=-0.02265; // paso de temperaturas
int NT=200; //Numero de temperaturas
double eprom[NT];
double e2prom[NT];
double mprom[NT];
double m2prom[NT];
//contadores 
for ( int j=0; j< NT;j++){
	eprom[j]=0;
	e2prom[j]=0;
	mprom[j]=0;
	m2prom[j]=0;
} 


double* ep=eprom;
double* e2p=e2prom;
double* mp=mprom;
double* m2p=m2prom;
//T0=.5;

int Nsites=2*N1*N2; // numero de sitios en toda la red.

//chequeo exportacion 
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
//chequeo termalizacion
double J2=0.;
double T=.5;
iniciar(p,q,N1,N2,J,J2,e,m);
plot_termalizacion( N1, N2, J, J2, T,p,q,e,m,tau);
*/


double J2=-0.0;///NNN coupling, positivo es ferro.
double Tc=1.519;//teorica, ferro y anti ferro para la honeycombe-NN


/*
double J2=-0.1;///NNN coupling, positivo es ferro.
double Tc=1.;//teorica, ferro y anti ferro para la honeycombe-NN
*/

/*
double J2=-0.2;///NNN coupling, positivo es ferro.
double Tc=0.42; // la saco del Cv, depende de J2.
*/

/*
double J2=-0.3;///NNN coupling, positivo es ferro.
double Tc=0.31; // la saco del Cv, depende de J2.
*/

/*
double J2=-0.4;///NNN coupling, positivo es ferro.
double Tc=0.45; // la saco del Cv, depende de J2.
*/

/*
double J2=-1.0;///NNN coupling, positivo es ferro.
double Tc=0.9; // la saco del Cv, depende de J2.
*/



char filenameA [100];
sprintf (filenameA,"/home/s-l/Desktop/arrays_%d_%.1f.txt",N1,-J2);
FILE *d = fopen(filenameA, "w"); // open file for export arrays

char filenameT [100];
sprintf (filenameT,"dataAF/termodinamica_%d_%.1f.txt",N1,-J2);
FILE *f = fopen(filenameT, "w"); // open file to escribir las cantidades termodinamicas


for(int ii=1;ii<=Nprom;ii++){// Nprom annealings
	iniciar(p,q,N1,N2,J,J2,e,m);
	AnnealingT(e,m,Nsites,Neq,T0,N1,N2,p,q,J,J2,d,Tc,ep,e2p,mp,m2p);
	printf("%d, m=%f, e=%f",ii,*m/(double)Nsites,*e/(double)Nsites);
}





//exporto valores medios de observables termodinamicos 
int k=0;
for (double T=T0;T>0.0226;T+=-0.02265){
	//printf("%f, %f, %f\n",T,*(mp+k)/Nprom/Nsites, *(ep+k)/Nprom/Nsites);
	*(ep+k)=*(ep+k)/(float)Nprom/Nsites; //tomo valor medio + normalizacion
	*(e2p+k)=*(e2p+k)/(float)Nprom/Nsites/Nsites; //tomo valor medio  + normalizacion
	//*(cvpoint+k)= Nsites * (  *(e2p+k) - pow(*(ep+k),2)   )/T/T; //Barkema
	*(mp+k)=*(mp+k)/(float)Nprom/Nsites; //tomo valor medio 
	*(m2p+k)=*(m2p+k)/(float)Nprom/Nsites/Nsites; //tomo valor medio 
	//*(chipoint+k)=Nsites*   (  *(m2p+k) - pow(*(mp+k),2)  )/T; //Barkema
	fprintf(f,"%f %f %f %f %f\n",
		T,*(mp+k),*(m2p+k),*(ep+k),*(e2p+k) ) ;
//*(chipoint+k) *(cvpoint+k)
k+=1;
}





printf("\n");
printf("J2=%.1f\n",-J2);



fclose(d); //cierro el file donde escribo las configuraciones
fclose(f);//cierro el file donde escribo la termodinamica



return 0;
} /* ---------------------- End Main ---------------------- */


void iniciar(int*p, int*q,int N1,int N2,double J, double J2,double*e ,double*m){
	srand((long)time(NULL)); /* initialize RNG rand() */
	state0(N1,N2,p); /* estado inicial subred A*/ 
	state0(N1,N2,q);/* estado inicial subred B*/ 
	//printf ("Energy0/site= %f\n",Energy/(double)(N1*N2)/2.);
	*e=calc_Energy(p,q,N1,N2,J,J2); //Energia de la configuracion inicial
	*m=calc_mag(p,q,N1, N2); //mag de la configuracion inicial
}
//----------------------- fin iniciar-



void export_data(int* p, int* q, int N1, int N2, double T, double Tc,FILE* d){

if(T<Tc ){fprintf(d,"%d ",0);}//label
else{fprintf(d,"%d ",1);};
fprintf(d,"%f ",T);//Temperature

for(int i=0;i<N1;i++){ 
	int j=mod(-2*i,N2); 
	for(int k=0;k<N2;k++){ // repito N2 veces 
		fprintf(d,"%d %d ",*(p+mod(j,N2)+i*N2),*(q+ mod(j,N2) +i*N2) );
		j=mod(j+1,N2);
	}
}
fprintf(d,"\n");

}//end export_data ---------------------------------------------------


void AnnealingT(double* e, double* m, int Nsites, int Neq,
double T0, int N1, int N2,int* p, int* q, double J,double J2,FILE* d,double Tc,
double*ep,double*e2p,double*mp,double*m2p){
int ii=0;
	for (double T=T0;T>0.0226;T+=-0.02265){
		//estimacion del equilibrio a la dada temperatura:
		for (int t=0;t<Nsites*Neq;t++){montecarlo_step(N1,N2,J,J2,T,p,q,e,m);}
		//printf("T=%f, e=%f, m=%f\n",T,*e/Nsites,*m/Nsites);
		if(abs(*m/Nsites)<1.1){
			export_data(p,q,N1,N2,T,Tc,d);
			*(ep+ii)+=*e;
			*(e2p+ii)+=pow(*e,2);
			*(mp+ii)+=abs(*m); // me quedo con la positiva para promediar 
			*(m2p+ii)+=pow(*m,2);
			}
		else{
			printf("E");
			//double EE=calc_Energy(p,q,N1,N2,J,J2); //Energia de la configuracion
			//double MM=calc_mag(p,q,N1, N2); //mag de la configuracion
			//printf("e=%f, m=%f",EE,MM);
			}
		ii+=1;
	}
}

//--------------------------------------------

void Annealing(double* e, double* m, int Nsites, int Neq,
double T0, int N1, int N2,int* p, int* q, double J,double J2,FILE* d,double Tc){

	for (double T=T0;T>0.0226;T+=-0.02265){
		//estimacion del equilibrio a la dada temperatura:
		for (int t=0;t<Nsites*Neq;t++){montecarlo_step(N1,N2,J,J2,T,p,q,e,m);}
		if(abs((float)*m/Nsites)<1.1){
			//export_data(p,q,N1,N2,T,Tc,d);
			}
		else{printf("ERROR");}
	}
}

//--------------------------------------------



void plot_termalizacion(int N1,int N2,double J,double J2,double T,int* p,int* q,double* e,double* m, int tau){
FILE *f = fopen("termalizacion.txt", "w"); // open file to write magnetizations vs time.
if (f == NULL){
    printf("Error opening file!\n");
    exit(1);
}

// estoy midiendo el tiempo en montecarlosteps
for (int time=0;time<=tau;time++){
	fprintf(f,"%d %f %f\n",time,*m/(double)(N1*N2)/2.,*e/(double)(N1*N2)/2.);
	for (int j1=1;j1<=2*N1*N2;j1++){montecarlo_step(N1,N2,J,J2,T,p,q,e,m);}
}

fclose(f);
}//end plot_termalizacion

//--------------------------------------------------------------------------------

//montecarlo step

void montecarlo_step(int N1, int N2, double J,double J2,double T,int *p,int *q,double* e,double* m){
int rand_i=rand_int(0,N1-1);//sorteo una celda unidad
int rand_j=rand_int(0,N2-1);//sorteo una celda unidad
//printf("%d, %d\n",rand_i, rand_j);
double Delta_E=0.; //Delta_E sale de la ecuacion 3.10 del Barkema
double beta=1./T;
int sumNN=0;
int sumNNN=0;
double epsilonE=1.E-5;


double r1=rand_double();//sorteo la subred
//printf("%f\n",r1);
if (r1<0.5){ //miro el espin Bij
	sumNN+= *(p +mod(rand_i,N1)*N2 + mod(rand_j,N2) );
	sumNN+= *(p + mod(rand_j+1,N2)+mod(rand_i,N1)*N2);
	sumNN+= *(p + mod(rand_j-1,N2)+mod(rand_i+1,N1)*N2 );

	sumNNN+=*(q +mod(rand_i-1,N1)*N2 + mod(rand_j+2,N2) );
	sumNNN+=*(q +mod(rand_i,N1)*N2   + mod(rand_j+1,N2) );
	sumNNN+=*(q +mod(rand_i+1,N1)*N2 + mod(rand_j-1,N2) );
	sumNNN+=*(q +mod(rand_i+1,N1)*N2 + mod(rand_j-2,N2) );
	sumNNN+=*(q +mod(rand_i,N1)*N2   + mod(rand_j-1,N2) );
	sumNNN+=*(q +mod(rand_i-1,N1)*N2 + mod(rand_j+1,N2) );

	Delta_E+=2.*  *(q + rand_j+rand_i*N2)  *(J* sumNN+J2* sumNNN);
	//printf("DeltaE1= %f\n",Delta_E);
	if(Delta_E<=epsilonE){ 
		*(q + rand_j+rand_i*N2) = - *(q + rand_j+rand_i*N2);
		*e+=Delta_E;
		*m+=-2.* *(q + rand_j+rand_i*N2) ; // staggered !!!!!!
	}// si baja la energia, flip.
	else{
		double r2=rand_double();
		if(r2<exp(-beta*Delta_E) ){
			*(q + rand_j+rand_i*N2) = -*(q + rand_j+rand_i*N2);
			*e+=Delta_E;
			*m+=-2.*  *(q + rand_j+rand_i*N2) ; // staggered !!!!!!
		}
	}//cierra el else
}//cierra el if de arriba
else{ // o sea r>0.5 // miro el espin Aij
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
	Delta_E+=2.* *(p + rand_j+rand_i*N2)  * (J* sumNN+J2* sumNNN);
	//printf("DeltaE2= %f\n",Delta_E);
	if(Delta_E<=epsilonE){ 
		*(p + rand_j+rand_i*N2) = -*(p + rand_j+rand_i*N2);
		*e+=Delta_E;
		*m+=2.*  *(p + rand_j+rand_i*N2);
}// si baja la energia, flip.
else{
	double r2=rand_double();
	if(r2<exp(-beta*Delta_E) ){
	*(p + rand_j+rand_i*N2) = -*(p + rand_j+rand_i*N2);
	*e+=Delta_E;
	*m+=2.*  *(p + rand_j+rand_i*N2);
	}
}
}
//printf("%f\n",exp(-beta*Delta_E));
//printf("Delta_E= %f\n",Delta_E);
//printf("E= %f\n",*e);
}//end montecarlo_step






//Calculo la energia de la configuracion 
double calc_Energy(int *p,int *q,int N1,int N2,double J,double J2){ 
printf("  \n");
double EnergyP=0;
for (int i=0 ;i<N1 ;i++ ){ //recorro todas las celdas 
for (int j=0 ;j<N2 ;j++ ){
EnergyP+=(-1)*J*( *(q + j+i*N2  ) * *(p + mod(j,N2)+mod(i,N1)*N2)  );
EnergyP+=(-1)*J*( *(q + j+i*N2  ) * *(p + mod(j+1,N2)+mod(i,N1)*N2) );
EnergyP+=(-1)*J*( *(q + j+i*N2  ) * *(p +mod(j-1,N2)+mod(i+1,N1)*N2) );


EnergyP+=(-1)*J2*( *(p +j+i*N2  ) * *(p +mod(i-1,N1)*N2 + mod(j+2,N2)  ));
EnergyP+=(-1)*J2*( *(p +j+i*N2  ) * *(p +mod(i,N1)*N2 + mod(j+1,N2)  ));
EnergyP+=(-1)*J2*( *(p +j+i*N2  ) * *(p +mod(i+1,N1)*N2 + mod(j-1,N2)  ));

EnergyP+=(-1)*J2*( *(q +j+i*N2  ) * *(q +mod(i-1,N1)*N2 + mod(j+2,N2)  ));
EnergyP+=(-1)*J2*( *(q +j+i*N2  ) * *(q +mod(i,N1)*N2 + mod(j+1,N2)  ));
EnergyP+=(-1)*J2*( *(q +j+i*N2  ) * *(q +mod(i+1,N1)*N2 + mod(j-1,N2)  ));
}
}

return EnergyP;

}//fin calc_Energy()

//Calculo la magnetizacion de la configuracion inicial
double calc_mag(int *p,int *q,int N1,int N2){ 
double magg=0.;
for (int i=0 ;i<N1 ;i++ ){ 
for (int j=0 ;j<N2 ;j++ ){
magg+=*(p + j+i*N2);//sublattice A
magg+=-*(q + j+i*N2);//sublattice B       staggered! 
}
}

return magg;

}//fin calc_Energy()




/* esta funcion inicializa el sistema aleatoriamente */
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
/* */

/* */
double rand_double2(double x, double y) {
   return (y - x)*rand_double() + x;
}
/* */

/* */
int rand_int(int a, int b) {
   return (int)((b - a + 1)*rand_double()) + a;
}
/* */

/*Le robe esta rutina a Ines que la robo de otro lado, funcion modulo para la indexacion*/
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


