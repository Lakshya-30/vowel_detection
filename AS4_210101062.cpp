// AS3_210101062.cpp : Defines the entry point for the console application.

#include "stdafx.h"
#include <iostream>
#include <cstring>
#include <limits>
#include "stdlib.h"
#include "math.h"
#include "ctype.h"
using namespace std;

#define framesize 320	//samples in one frame
#define P 12 
#define F 5				//number of stable frames
#define PI 3.142857142857

double steadyFrames[F][framesize], R[F][P+1], A[F][P+1], C[F][P+1], avgCi[25][P+1], allCi[100][F][P+1], refCi[F][P+1], tokhuraDist[5];
double samples[100000] = {0}, energy[10000] = {0};
int energySize, sampleSize = 0;
const double limit = 5000.0;
double dcShift=0, nFactor=1;
int steadyStart=0, steadyEnd=0;
double tokhuraWeights[]={1.0, 3.0, 7.0, 13.0, 19.0, 22.0, 25.0, 33.0, 42.0, 50.0, 56.0, 61.0};			//tokhura's distance
int u = 0, v = 0, files = 0;		//index helper to store the values of ci to allCi
char vowels[5] = {'a', 'e', 'i', 'o', 'u'};			//vowels
int totalCorrect = 0, individualCorrect= 0;			//variable to find accuracy
double accuracy[5];

//function to compute the DC offset of the speech signal
void getDCShift(){
	for(int i=0; i<sampleSize; i++)
	{
		dcShift += samples[i];
	}
	dcShift /= sampleSize;
}

//function to compute the normalization factor of the speech signal
void getNFactor(){
	int max = INT_MIN;
	for(int i=0; i<sampleSize; i++)
	{
		if(abs(samples[i] )> max){
			max = abs(samples[i]);
		}
		nFactor = (double)limit/max;
	}
}

//function for applying hamming window to all the stable frames
void hammingWindow(){
	for(int i=0; i<F; i++){
		for(int j=0; j<framesize; j++){
			steadyFrames[i][j] *= 0.54-0.46*cos(2*PI*steadyFrames[i][j]/framesize-1);
		}
	}
}

//function to find the steady frames
void findSteadyFrames(){
	int highestEnergyIndex = 0, energySize = 0;
	double en = 0, maxEnergySoFar = 0;

	//calculating and marking the highest energy frame
	for(int i=0; i<sampleSize; i++){
		//if the n is a multiple of N then we have to store energy to the array
		if(i%framesize == 0 && i!=0){
			en /= framesize;				//taking average
			if(maxEnergySoFar < en){
				maxEnergySoFar = en;
				highestEnergyIndex = energySize;
			}
			energy[energySize++] = en;	//store the energy value
			en = 0;						//resetting en to 0
		}
		en += samples[i] * samples[i];
	}
	steadyStart = highestEnergyIndex > 2 ? (highestEnergyIndex-2)*framesize : 0;										//marking start
	steadyEnd = highestEnergyIndex <energySize - 3 ? (highestEnergyIndex+3)*framesize : energySize*framesize;			//marking end
	
	int f=0, j=0;
	for(int i = steadyStart; i<steadyEnd; i++){			//storing the frames into steadyFrames array
		steadyFrames[f][j++] = samples[i];
		if(j % framesize==0){
			f++;
			j=0;
		}
	}
}

//function to apply the raised Sin Window to the Ci's of each frame
void raisedSinWindow(){
	long double sum=0;
	for(int f = 0; f<F; f++){
		for(int m=1;m<=P;m++){
			sum = (P/2)*sin((PI*m)/P);
			C[f][m]*=sum;
		}
	}
}

//function to compute the Ri values
void computeRis(){
	for(int f = 0; f<F; f++){
		for(int m =0; m<=P; m++){
			R[f][m] = 0;
			for(int k=0; k<framesize-m; k++){
				R[f][m] += steadyFrames[f][k]*steadyFrames[f][k+m];
			}
		}
	}
}

//function to compute the Ai values using Levinson Durbin algorithm
void computeAis(){
	double Alpha[13][13],E[13],K[13];
	double sum=0;

	for(int f = 0; f<F; f++){
		if(R[f][0] == 0){
			cout<<"Divide by 0 error\n";
			break;
		}
		E[0] = R[f][0];			    //step 1
		for(int i=1;i<=P;i++){
			sum=0;
			for(int j=1;j<=i-1;j++){
				sum += Alpha[i-1][j]*R[f][i-j];	
			}
			
			K[i]=(R[f][i]-sum)/E[i-1];			//step 2
			Alpha[i][i]=K[i];					//step 3
			for(int j=1;j<=i-1;j++){
				Alpha[i][j]=Alpha[i-1][j] - K[i]*Alpha[i-1][i-j];	//step 4
			}
			E[i]=(1-(K[i]*K[i]))*E[i-1];		//step 5
		}				

		for(int i=1;i<=P;i++){
			A[f][i]= Alpha[P][i];
		}
	}
}

//function to store the values of Ci for 5 frames in the allCi array
void storeCito3D(){
	for(int f=0; f<F; f++){
		for(v=0;v<P;v++){
			allCi[files][f][v+1]=C[f][v+1];
		}
	}
	files++;
}

//function to store the avg ci values to file
void computeAvgCi(){
	FILE *ref;
	char filename[30];
	int index = 0;
	
	for(int v=0; v<5; v++){			//loop over vowel to save file/vowel 
		sprintf(filename, "AvgCi/reference_ci_%c.txt", vowels[v]);
		ref = fopen(filename, "w");
		for(int f=0; f<F; f++){			//looping over frames
			for(int p=0; p<P; p++){		//looping over p
				double sum = 0;
				//summing up all the file wise values and taking average
				for(int file=v*20; file<(v+1)*20; file++){
					sum += allCi[file][f][p+1];
				}
				sum /= 20.0;
				avgCi[index][p+1] = sum;
				fprintf(ref, "%lf ", sum);
			}
			index++;
			fprintf(ref, "\n");
		}
		cout<<filename<<" has been created\n";
		fclose(ref);
	}
}

//This function calulate the cepstral coeff Ci's
void calculateCis(){
	double sum=0;
	for(int f = 0; f<F; f++){
		C[f][0]=log(R[f][0]*R[f][0]);

		for(int m=1;m<=P;m++){
			sum=0;
			for(int k=1;k<m;k++){
				sum += (k*C[f][k]*A[f][m-k])/(m*1.0);
			}
			C[f][m]=A[f][m]+sum;
		}
	}
	
	raisedSinWindow();		//applying raised sin window on cis
	storeCito3D();			//storing the ci values to 3D matrix
}

//function to execute training
void training(){
	char filenames[5][30] = {"vowels/210101062_a_01.txt", "vowels/210101062_e_01.txt", "vowels/210101062_i_01.txt", "vowels/210101062_o_01.txt", "vowels/210101062_u_01.txt"};

	for(int i=0; i<5; i++){			//loop over vowels
		for(int j = 0; j<20; j++){	//loop over files
			//read file
			FILE *f1 = NULL;    // Read the data (x) from a text file (a user input)
			double x;
			int err = fopen_s(&f1, filenames[i], "r");
			if(err != NULL)     //error handling for file opening
			{
				printf("\nCannot open the file\n");
				system("pause");
				exit(1);
			}   
			sampleSize = 0;
			while( !feof(f1) )
			{
				if( fscanf_s( f1, "%lf", &x) == 1)
				{
					samples[ sampleSize ] = x;
					sampleSize++;
				}
				else{
					printf("Error in line %d\n", sampleSize);
				}
			}

			cout<<"Training "<<filenames[i]<<endl;
			getDCShift();
			getNFactor();
			for (int i = 0; i <= sampleSize; i++)					//DC offset correction and normalization of speech signal
			{
				samples[i] = (samples[i]-dcShift)*nFactor;
			}

			findSteadyFrames();
			hammingWindow();
			computeRis();
			computeAis();
			calculateCis();
			fclose( f1 );

			if(filenames[i][20] == '9'){	//next file
				filenames[i][19]++;
				filenames[i][20] = '0';
			}
			else filenames[i][20]++;					
		}
	}
}

//fucntion to calculate the Tokhura distance using dump Ci values of training set
double tokhuraDistance(FILE *ip){
	char line[3000];
	int f = 0;
	//read the reference file and store in restoredCi
	while(!feof(ip) && f<F){
		fgets(line, sizeof(line), ip);
		sscanf(line, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf", &refCi[f][1], &refCi[f][2], &refCi[f][3], &refCi[f][4], &refCi[f][5], &refCi[f][6], &refCi[f][7], &refCi[f][8], &refCi[f][9], &refCi[f][10], &refCi[f][11], &refCi[f][12]);
		f++;
	}

	double finalDist = 0;
	for(int i=0; i<F; i++){
		double dist = 0;
		for(int p=1; p<=P; p++){
			double d = (C[i][p]- refCi[i][p]);
			dist += tokhuraWeights[p-1]*d*d;
		}
		finalDist += dist/(P*1.0);	//taking average
	}
	return finalDist/(F*1.0);	//taking average
}

//function to make prediction based on the Tokhura distance
char calculateTokhura(){
	char filename[30];
	FILE *ip;
	double minDist = DBL_MAX;
	char predVowel;

	for(int i=0; i<5; i++){
		sprintf(filename, "AvgCi/reference_ci_%c.txt", vowels[i]);
		ip = fopen(filename, "r");
		if(ip == NULL) printf("Error in opening file %s\n", filename);
		double distance = tokhuraDistance(ip);			//calculating distance from vowels[i]'s reference file
		tokhuraDist[i] = distance;						//storing distance to print
		if(minDist > distance){							//checking whether it is minimum than the minimum distance we found?
			minDist = distance;				//update incase new distance is minimum
			predVowel = vowels[i];		//update the vowel to be returned
		}
	}
	return predVowel;
}

//testing files
void testing(){
	char filename[5][30] = {"vowels/210101062_a_21.txt", "vowels/210101062_e_21.txt", "vowels/210101062_i_21.txt", "vowels/210101062_o_21.txt", "vowels/210101062_u_21.txt"};
	files = 0;		//making files variable 0 to store the ci values in 3D array

	for(int v = 0; v<5; v++){
		individualCorrect = 0;
		for(int file=0; file<10; file++){
			//reading test file
			FILE *f1 = NULL;    // Read the data (x) from a text file (a user input)
			double x;
			int err = fopen_s(&f1, filename[v], "r");
			if(err != NULL)     //error handling for file opening
			{
				printf("\nCannot open the file\n");
				system("pause");
				exit(1);
			}
			sampleSize = 0;
			while( !feof(f1) )
			{
				if( fscanf_s( f1, "%lf", &x) == 1)
				{
					samples[ sampleSize ] = x;
					sampleSize++;
				}
				else{
					printf( "Error in line %d\n", sampleSize  );
				}
			}
			
			cout<<"Testing: "<<filename[v]<<endl;

			getDCShift();
			getNFactor();
			cout<<"DC shift offset: "<<dcShift<<endl;
			cout<<"Normalization factor: "<<nFactor<<endl;
			for (int i = 0; i <= sampleSize; i++)					//DC offset correction and normalization of speech signal
			{
				samples[i] = (samples[i]-dcShift)*nFactor;
			}

			findSteadyFrames();
			//cout<<"Steady frame start sample: "<<steadyStart<<endl;
			//cout<<"Steady frame end sample: "<<steadyEnd<<endl;

			hammingWindow();
			computeRis();
			computeAis();
			calculateCis();

			char prediction = calculateTokhura();						//tokhura calculation will return the predicted character
			printf("Predicted vowel: %c\n", prediction);
			printf("Tokhura weights distance from (a, e, i, o, u) = (%lf, %lf, %lf, %lf, %lf)\n\n", tokhuraDist[0], tokhuraDist[1], tokhuraDist[2], tokhuraDist[3], tokhuraDist[4]);
			if(prediction == vowels[v]) totalCorrect++, individualCorrect++;		//if prediction is correct then increasing counts
			fclose( f1 );
			if(filename[v][20] == '9'){	//next file
				filename[v][19]++;
				filename[v][20] = '0';
			}
			else filename[v][20]++;
		}
		accuracy[v] = (individualCorrect/10.0)*100;		//calculating the vowel recognition accuracy/vowel
	}
	//printing accuracy
	printf("---------------------------------------------------------------\n");
	printf("Accuracy of vowel a: %.2lf %% \n", accuracy[0]);
	printf("Accuracy of vowel e: %.2lf %% \n", accuracy[1]);
	printf("Accuracy of vowel i: %.2lf %% \n", accuracy[2]);
	printf("Accuracy of vowel o: %.2lf %% \n", accuracy[3]);
	printf("Accuracy of vowel u: %.2lf %% \n", accuracy[4]);
	printf("Overall Accuracy of the system: %.2lf %% \n", (totalCorrect/50.0)*100);
}

int main()
{
	training();
	cout<<endl;
	computeAvgCi();
	cout<<endl;
	testing();
	return 0;
}