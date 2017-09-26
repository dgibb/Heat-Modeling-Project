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

    if(argc== 5){       // checks if there was an optional ambient temperature file

        ambientFile= fopen( argv[3], "r");
        assert(ambientFile != NULL);
        fscanf( ambientFile, "%lf", &ambient); // scans in ambient temperature
        printf( "Ambient =  %d\n", ambient); // check delete later
        paraFile = fopen(argv[1],"r");
        powerTraceFile = fopen(argv[2],"r");
        outputFile = fopen(argv[4],"wt");
        printf("loaded four  files\n");

    } else {

        ambient = 300;// ambient is 300 K
        paraFile = fopen(argv[1],"r"); /// open all files
        powerTraceFile = fopen(argv[2],"r");
        outputFile = fopen(argv[3],"wt");
    }

    assert(paraFile != NULL); // checks for the themal parameter file
    assert( powerTraceFile != NULL); // checks for the power trace file
    assert( outputFile != NULL); // checks for the output file

    int i, j;

    printf("creating capacitor array size [4]\n");
    for( i=0; i<4; i++){ // created capacitor array
        fscanf(paraFile, "%lf", &cap[i]);
        printf("[%lf], ", cap[i]);
    }
    printf(" \n");
    printf(" \n");

    printf("creating resistor matrix size [4][4]\n");
    for( i=0; i<4; i++){// create resistor matrix
        for( j=0; j<4;j++){
            fscanf(paraFile,"%lf", &rest[i][j]);
            printf("[%lf], ", rest[i][j]);
        }
        printf(" \n");
    }
    printf(" \n");
    printf(" \n");

    printf("setting initial temperature...\n");
    for( i=0; i<4; i++){ // set the inital condition
        temp[i] = ambient;
        printf("[%lf], ", temp[i]);
    }
    printf(" \n");
    printf(" \n");

    printf("setting initial time...\n");
    t =0; //time =0;
    printf(" \n");

    int count = 0;

    while ( fscanf(powerTraceFile,"%lf\n", &time)!=EOF) {// calls the RK for the entire file
    printf("Running rk at time = %lf \n", t);
    printf("power at each core =  ");
        for(  i=0; i<4; i++) {  // gets the power
            fscanf(powerTraceFile, "%lf",&power[i]);
            printf("[%lf], ", power[i]);
        }
        printf(" \n");
        rk( time );
        printf(" \n");
    count++;
    }

    printf("done \n");
    fclose(powerTraceFile);  // close all files
    fclose(outputFile);
    fclose(paraFile);
    return 0;
}

void rk(double n){

  double k[4][4];  //values of k[x][core#]
  double h=0.005;
    int x, i, j;

     while (t<=n){

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

        printf("calculating temp\n");
        for (x=0; x<4; x++){        //for each core find y t+h
            temp[x]=temp[x]+(k[0][x]+2*k[1][x]+2*k[2][x]+k[3][x])/6;//assumes initial conditions already in a[x]
            printf("[%lf]", temp[x]);
        }
        printf("\n");

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

void age(double ti){ // age of each core

    double tempAmbient;
    double output[4];
    double tempDevice[4];
    int i;
    tempAmbient = alphaFunction(ambient);
    fprintf(outputFile, "Time = %lf\n", time);

    for(i = 0; i < 4; i++){
        tempDevice[i] = alphaFunction(temp[i]);
        output[i] = (tempDevice[i]/tempAmbient);
        fprintf(outputFile, "core %d: temp = %lf, output = %lf\n", i, temp[i], output[i]);
    }
}

double f( double t, double y, int x){
    double var = temp[x];
    temp[x]=y;
    int i, j;
    double change[4];
    for(  i=0; i<4; i++){
        for(j =0; j<4; j++){
            if( j!=i){
                change[i] = (((temp[i]-temp[j])/(rest[i][j]*cap[i]))+(power[i]/cap[i]));
            }
        }
    }
    temp[x]=var;
    return(change[0] + change[1] + change[2] + change[3]);
}
