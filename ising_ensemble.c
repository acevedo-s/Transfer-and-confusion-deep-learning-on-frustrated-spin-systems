#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "pcg_basic.h" //(for this code to work you need to compile it together with the file "pcg_basic.c")

void initialize(int *sentinel, int size){
	int filler(void);
	srand(time(0));
	for(int count=0;count<size;count++){
		*sentinel= filler();
		sentinel++;
	}
}

int filler(void){
	double randnum=rand()/(double)RAND_MAX;
	if (randnum>=0.5)
		return 1;
	else
		return -1;
}



int mod (int a, int b) //not the C-defined mod function
{
   if (b < 0) //you can check for b == 0 separately
     return mod(a, -b);   
   int ret = a % b;
   if (ret < 0)
     ret+=b;
   return ret;
}



/*Ising Model*/
int main(){
	int L=30; //lattice with LxL sites
	int N=L*L;
	int thermalization_steps=900*N;
	double J1=1; //H= J sum(S_i*S_j) //first-neighbor coupling
	double J2=0.46; //second-neighbor coupling

	/* we declare the lattice */
	int lat[L][L];

	pcg32_random_t rngptr1, rngptr2, myrng; //for the random number generator
	

	int i,j; //this will be the indexes used to identify the site (i,j) that will be modified
	FILE *pointer_lat; //pointer used to save the configurations in a text file


	/* --- MONTECARLO --- */
	pointer_lat= fopen("0.46.txt", "w");
	int numdata=200; // 200 temperature values
	double deltaE;
	for (int zzz=1; zzz<401; zzz++){ //400 independent initial configurations
		printf("%d\n",zzz);
		initialize(&lat[0][0],N);
		//int label=1; //comment the label-related lines if you do not know the Tc value
		for(int xx=numdata; xx>=1; xx--){
			
			//double T = xx * 4.53 / numdata ; //for the AFM square lattice
			double T = xx * 4. / numdata ;
			//double T = xx * 2.4 / numdata ; //for the AFM triangular lattice (to have a balanced data set)
			double beta = 1./T; //k=1

			//if (T <= 3.640957){
			//label=0;} //Tc=3.640957 for FM on TRIANGULAR lat
			//if (T <= 2.269185){
			//label=0;} //Tc=2.269185 for FM/AFM on SQUARE lat
			//if (T <= 1.2){
			//label=0;} //Tc=1.2 for AFM on TRIANGULAR lat (obtained from another Montecarlo)

			//labels are used by the neural networks to classify spin configurations
			
			/* --- WE LET THE LATTICE THERMALIZE --- */
			for(int steps=0; steps<thermalization_steps; steps++){
				
					/* choose site (i,j) */
					int i = pcg32_boundedrand_r(&rngptr1, L);
					int j = pcg32_boundedrand_r(&rngptr2, L);
					double rand_value = ldexp(pcg32_random_r(&myrng), -32);
					/* METROPOLIS */
					deltaE= -2*J1*lat[i][j]*(lat[i][mod(j-1,L)]+lat[i][mod(j+1,L)]+lat[mod(i-1,L)][j]+lat[mod(i+1,L)][j]);
					//deltaE= deltaE - 2*J1*lat[i][j]*(lat[mod(i-1,L)][mod(j-1,L)]+lat[mod(i+1,L)][mod(j+1,L)]); //uncomment for TRIANGULAR
					deltaE= deltaE - 2*J2*lat[i][j]*(lat[mod(i-1,L)][mod(j-1,L)]+lat[mod(i+1,L)][mod(j+1,L)]); //uncomment for SQUARE (1ST- AND 2ND-NEIGHBOR INT)
					deltaE= deltaE - 2*J2*lat[i][j]*(lat[mod(i-1,L)][mod(j+1,L)]+lat[mod(i+1,L)][mod(j-1,L)]); //uncomment for SQUARE (1ST-NEIGHBOR INT)
					if (deltaE<0 || rand_value<exp(-beta*deltaE))
						lat[i][j]= -lat[i][j];
					///---------------
			}

			//fprintf(pointer_lat, "%d ", label); //comment the label-related lines if you do not know the Tc value
			fprintf(pointer_lat, "%f ", T);
			for (int ii=0;ii<L;ii++){
					for (int jj=0;jj<L;jj++){
						fprintf(pointer_lat, "%d ", lat[ii][jj]);
					}
			}
			fprintf(pointer_lat, "\n");
		}
	}
	
	
	fclose(pointer_lat);
	return 0;
}