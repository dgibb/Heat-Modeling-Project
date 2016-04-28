//David Gibb and Daniel Mathieu
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

	FILE *powerTraceFile; // pointer to all input files
	FILE *outputFile;
	FILE *ambientFile;
	FILE *paraFile;
	double cap[4]; 
	double rest[4][4];
	double power[4];
	double temp[4];
	
	double time;	
	double ambient;
	double t;
	void rk(double n);
	double alphaFunction(double temp);
	double f(double t, double y, int x);
	void age(double ti);

int main (int argc, char *argv[]){

	
	assert(argc>= 3);// all files have been inputed

	
	if(argc== 5){ // checks if there was an optional ambient temperature file
	ambientFile= fopen( argv[3], "r");
	assert(ambientFile != NULL);
	fscanf( ambientFile, "%f", &ambient); // scans in ambient temperature
	printf( "Ambient =  %d\n", ambient); // check delete later
	paraFile = fopen(argv[1],"r");
	powerTraceFile = fopen(argv[2],"r");
	outputFile = fopen(argv[4],"wt");
	printf("loaded four  files\n");} // check delete later

	else { // if no ambient load all files 
	ambient = 300;// ambient is 300 K
	paraFile = fopen(argv[1],"r"); /// open all files 
	powerTraceFile = fopen(argv[2],"r");
	outputFile = fopen(argv[3],"wt");
	} 

	assert(paraFile != NULL); // checks for the themal parameter file
	assert( powerTraceFile != NULL); // checks for the power trace file
	assert( outputFile != NULL); // checks for the output file
	int i;
	int j;
	for( i=0; i<4; i++){ // created capacitor array
	fscanf(paraFile, "%lf",&cap[i]);
	}
	
	for( i=0; i<4; i++){// create resistor matrix
	for( j=0; j<4;j++){
	fscanf(paraFile,"%lf", &rest[i][j]);}}
		
	for( i=0; i<4; i++){ // set the inital condition
	temp[i] = ambient;}
	//time =0;
	 t =0;
	
	while ( fscanf(powerTraceFile,"%f\n", &time)== 1) {// calls the RK for the entire file
	for(  i=0; i<4; i++){// gets the power
	fscanf(powerTraceFile, "%lf",&power[i]);
	printf(" power %lf", power[i]);}
		printf( " time %lf\n", time);
	rk( time );
	} 
	printf("done");
	fclose(powerTraceFile);  // close all files
	fclose(outputFile);
	fclose(ambientFile);
	fclose(paraFile);
	return 0;
		
} 
void rk(double n){
    double k[4][4];  //values of k[x][core#]
	
 	  double h=0.005; 
  
	int x;
	int i, j;

 	while (t<=n){    //while t<n  
	
       for (x=0; x<4; x++){        //for each core find k1

		k[0][x]=h*f(t,temp[x],x);           
        }
	for (x=0; x<4; x++){        //for each core find k2
	
		k[1][x]=h*f(t+h/2,temp[x]+k[0][x]/2,x);           
	}
	for (x=0; x<4; x++){        //for each core find k3
		k[1][x]=h*f(t+h/2,temp[x]+k[1][x]/2,x);           
	}
	for (x=0; x<4; x++){        //for each core find k4
		k[3][x]=h*f(t+h,temp[x]+k[2][x],x);           
	}
	for (x=0; x<4; x++){        //for each core find y t+h
	
		temp[x]=temp[x]+(k[0][x]+2*k[1][x]+2*k[2][x]+k[3][x])/6;//assumes initial conditions already in a[x}
	}
	t=t+h;
	    
    }
	age(t);// call age function      
}
double alphaFunction(double temp){ // age acceleration factor
double Ea = .8;
double k = .00008617;

double alpha = exp(-Ea/(k*temp));
return alpha;

}

void age(double ti){ // age of the core
	double tempAmbient;
	double output[4];
	double tempDevice[4];
	int i;
	tempAmbient = alphaFunction(ambient);
	fprintf(outputFile, "%f ", time);
	for(i = 0; i < 4; i++){
		tempDevice[i] = alphaFunction(temp[i]);
		output[i] = (tempDevice[i]/tempAmbient);
		fprintf(outputFile, "%f ", temp[i]);
		fprintf(outputFile, "%lf ", output[i]);
		fprintf(outputFile, "\n");
	}
}

double f( double t, double y, int x){
	double var = temp[x];
	temp[x] =y;
	int i;
	int j;	
	double change[4];	
	for(  i=0; i<4; i++){
	for(j =0; j<4; j++){
	if( j!=i){
	change[i] = -(((temp[i]-temp[j])/(rest[i][j]*cap[i]))+(power[i]/cap[i]));}}}
	temp[x]=var;
	return(change[0] + change[1] + change[2] + change[3]); 
	} 


