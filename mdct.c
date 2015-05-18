#include <stdio.h>
#include <stdlib.h>
#include <values.h>
#include <math.h>


void printShArray(short* a, int N);
void printFlArray(long double* a, int N);


void twoTriangles32768(short* *arr, int* N) {
    int i;
    short* a;

    *N = 4 * MAXSHORT + 3;

    a = calloc(sizeof(short), *N);
    *arr=a;
    for (i = 0; i < MAXSHORT; i++) {
        a[i]            = i;
        a[MAXSHORT + i] = MAXSHORT - i;
    }

    for (i = 0; i <= MAXSHORT+1; i++) {
        a[2*MAXSHORT+i]     = -i;
        a[3*MAXSHORT+i+1]   = MINSHORT + i;
    }
    printf("%d", *N);
    //printShArray(*arr, *N);
}
void twoTriangles(short* *arr, int n, int* N) {
    int i;
    short* a;

    *N = 4 * n;

    a = calloc(sizeof(short), *N);
    *arr=a;
    for (i = 0; i < n; i++) {
        a[i]            = i;
        a[n + i] = n - i;
    }

    for (i = 0; i < n; i++) {
        a[2*n+i]     = -i;
        a[3*n+i]   = -n + i;
    }
    //printf("%d", *N);
    //printShArray(*arr, *N);
}

void zeroOne(short* *arr, short step, short K, short amplitude, int* N){
    int i, j;
    short* a;
    *N = (int)step * K *2;
    a = calloc(sizeof(short), *N);
    *arr=a;
    for (i=0;i<K;i++)
        for (j=0;j<step;j++) {
            a[2 * i * step + j]     = 0;
            a[(2 * i + 1) * step + j] = amplitude;
        }
}
void x64(short* *arr, short x, short K, int* N){
    int i, j;
    *N = x * K;
    short* a;
    a = calloc(sizeof(short), *N);
    *arr=a;
    for (i=0;i<K;i++)
        for (j=0;j<x;j++)
            a[i * x + j] = 64*j;
}
void rand128K(short* a, short K){
    int i, j, N = 128 * K;
    a = calloc(sizeof(short), N);

    for (i=0;i<N;i++)
        a[i] = (short) ((1 - (2 * (rand() % 2))) * (rand() % MAXSHORT));//RANDOM NUMBER BETWEEN -32767 nad +32767
}

void kSin(short* a, short K){
    int i;
    a = calloc(sizeof(short), K);

    for (i=0;i<K;i++)
        a[i]= round(256.* sin(i*M_PI/256));
}


void kCos(short* * arr, int K, int* N){
    int i;
    short* a;
    *N=K;
    a = calloc(sizeof(short), *N);
    *arr=a;

    for (i=0;i<K;i++)
        a[i]= round(256.* cos(i*M_PI/(K/4.)));
}

void writeInputDataToFile(){

     /*   K=atoi(argv[0]);
        filename= argv[1];
        wavFile = fopen(fileName, "r");
            if (wavFile ==  NULL){

        }
        fread(N, sizeof(int), 1, wavFile);
        x=calloc(n, sizeof(int));

        for (i=0;i<N;i++)
            fread(x[i], sizeof(int), 1, wavFile);

        fclose(wavFile);

    */
}

void loadData(FILE* inputFile, char* address, short* *arr){
 //   K=atoi(argv[0]);
 //   filename= argv[1];
 //   inputFile = fopen(fileName, "r");
 //   if (wavFile ==  NULL){

//    }
//    fread(N, sizeof(int), 1, wavFile);
//    x=calloc(n, sizeof(int));

//    for (i=0;i<N;i++)
//        fread(x[i], sizeof(int), 1, wavFile);

//    fclose(wavFile);
}
void extendArrayForMdct(short* *arr, short N){
    short* a=calloc(2*N, sizeof(short));
    short i;
    for (i=0;i<N/2; i++)
        a[i]=a[3*N/2+i]=0;
    for (i=0;i<N;i++)
        a[i+N/2]=(*arr)[i];
    *arr=a;
}

void mdct(short* x, int N, long double** XPointer, short K) {
    int k, n;
    long double* X;
    X = calloc(K, sizeof(long double));

   // printShArray(x, 2*N);

    *XPointer= X;
    for (k = 0; k < K; k++) {
        X[k] = 0;
//        printf("");
        for (n = 0; n < 2 * N; n++)
            X[k] += (long double)x[n] * cos(M_PI * (n + .5 * (1.+N) ) * (k + .5)/ N );
    }
//    printFlArray(X, K);
    //getchar();
}
void imdct(long double* X, int K, long double** xPointer, int N) {
    int k, n;
    long double* x;
    x = calloc(2*N, sizeof(long double));

    //printFlArray(X, K);


    *xPointer= x;
    for (n = 0; n < 2*N; n++) {
        x[n] = 0;
        for (k = 0; k < K; k++) {
            x[n] += X[k] * cos(M_PI * (n + .5 * (1. + N)) * (k + .5)/ N);
     //       printf("%d - argument-%f, cosinus- %f, x[%d]=%f\n", k, (n + .5 * (1. + N)) * (k + .5) / N, cos(M_PI * (n + .5 * (1. + N)) * (k + .5) / N), n, x[n]);
        }
        x[n]/=N;
       // printf("x[%d]=%f\n", n, x[n]);
       // getchar();
    }
  //  printFlArray(x, 2*N);
  //  getchar();
}

void printShArray(short* a, int N){
    int i;
    printf("Array of Shorts\n%d", a[0]);
    for (i=1;i<N; i++){
        printf(", %d", a[i]);
        if (i%5000==0) {printf("more...");getchar();}
    }
    getchar();
    printf("\n");
}

void printFlArray(long double* a, int N){
    int i;
    printf("array of long doubles\n%Lf, \n", a[0]);
    for (i=1;i<N; i++){
        printf("%Lf, \n", a[i]);
      //  if (i%5==0) printf("\n");
        if (i%1000==0) {printf("more...");getchar();}
    }
    printf("\n");
}
void transformResults(long double* a, int N){
    int i;
    char svakitreci=0;
    printf("Nonzero values\n");
    for (i=0;i<N; i++)
        if (fabs(a[i]>0.00001)) {
            printf("a[%d]=%Lf\t", i, a[i]);
            svakitreci++;
            if (svakitreci%3==0) printf("\n");
        }
    printf("\n");
}
/*    printf("\nResults of the transform:\n");
    printf("%f", X[0]);
    for (i = 1; i < K; i++)
        printf(", %f", X[i]);
    printf("\n");*/


int main(int argc, char* argv[]) {
    short *x;
    long double *X, *X2, *X3, *invX, *invX2, *invX3;
    int i, k, n, N, K=1024;
    FILE* wavFile=NULL;
    char* fileName=NULL;
    char* fileFormat, *fileSize, *soundFormat, *fmt, *chunkSize, *PCMHeader, *numChannels, *sampleRate, *byteRate,*blockSize, *bitRate, *dataText, *dataSize;


    //Loading data from a file
    /*
        K=atoi(argv[0]);
        filename= argv[1];
        wavFile = fopen(fileName, "r");
            if (wavFile ==  NULL){

        }
        fread(N, sizeof(int), 1, wavFile);
        x=calloc(n, sizeof(int));

        for (i=0;i<N;i++)
            fread(x[i], sizeof(int), 1, wavFile);

        fclose(wavFile);

    */


  /*  x64(&x, 4, 64, &N);
    twoTriangles(&x, 64, &N);
//    kCos(&x, 256, &N);
 //   zeroOne(&x, /*128*/16, /*512*/8, 256, &N);
     //   printf("%d\n", );*/
/*
    printf("Orioginal array:\n");
    printShArray(x,N);
*/

    extendArrayForMdct(&x, N);
    printShArray(x,2* N);

    //N=N/2;

    //printf("%d\n", N);
    mdct(x, N/2, &X, N/256);
    mdct(x+N/2, N/2, &X2, N/256);
    mdct(x+N, N/2, &X3, N/256);


    imdct(X, N/256, &invX, N/2 );
    imdct(X2, N/256, &invX2, N/2 );
    imdct(X3, N/256, &invX3, N/2 );

    printf("0");
    for (i=1;i<N;i++)
        printf(", %d", i);
    printf("\n");

    for (i=0;i<N/2;i++)
        invX[i]=invX[i+N/2]+invX2[i];
    for (i=N/2;i<N;i++)
        invX[i]=invX2[i]+invX3[i-N/2];

    printf("MDCT+IMDCT result:\n");
    printFlArray(invX, N);

    printf("Difference:\n");
    for (i=0;i<N;i++) {
        int counter = 0;
        if (fabs(invX[i] - (long double) x[N / 2 + i] > 0.0)) {
            counter++;
            printf("difference[%d]=%Lf  ", i, invX[i] - x[N / 2 + i]);
            if (counter%3==0) printf("\n");
        }
    }
    printf("\n");
    /*printf("\nResults of the transform:\n");

    transformResults(X, K);*/



/*    printf("Enter file name:");
    scanf("%s", fileName);

    wavFile = fopen(fileName, "r");
    if (wavFile ==  NULL){

    }

// LOADING HEADER - FILE INFO

    fread(fileFormat, sizeof(char), 4, wavFile);
    if (fileFormat!= "RIFF"){

    }

    fread(fileSize, sizeof(char), 4, wavFile);
    fread(soundFormat, sizeof(char), 4, wavFile);
    if (soundFormat!=  "WAVE"){

    }

    fread(fmt, sizeof(char), 4, wavFile);
    if (fmt !=  "fmt "){

    }

    fread(chunkSize, sizeof(char), 4, wavFile);
    fread(PCMHeader, sizeof(char), 2, wavFile);
    fread(numChannels, sizeof(char), 2, wavFile);
    if ((numChannels !="1") && (numChannels!= "2")){

    }

    fread(sampleRate, sizeof(char), 4, wavFile);
    fread(byteRate, sizeof(char), 4, wavFile);
    fread(blockSize, sizeof(char), 2, wavFile);
    fread(bitRate, sizeof(char), 2, wavFile);
    fread(dataText, sizeof(char), 4, wavFile);
    if (dataText !=  "data"){

    }
    fread(dataSize, sizeof(char), 4, wavFile);

*/
// HEADER LOADED


  /*  printf("Enter N number of discrete samples:");
    scanf("%d", &N);
    printf("Enter samples in time domain:");

    printf("Enter K (K<=N number of frequences:");
    scanf("%d", &K);
    printf("Enter samples in time domain:");
*/

  /*  x = calloc(2 * N, sizeof(long double));
    for (i = 0; i <= 2 * N - 1; i++)
        scanf("%f", &x[i]);
*/

}