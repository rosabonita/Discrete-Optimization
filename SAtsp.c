/* Taken from Numerical Recipes p. 444--Simulated Annealing for TSP
All added comments refer to the line directly above unless noted otherwise.
When you compile, type: "gcc SAtsp.c rngs.c rvgs.c nrutil.c -lm"  */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "nrutil.h"
#include "rngs.h"
#include "rvgs.h"
#define ALEN(a,b,c,d) sqrt(((b)-(a))*((b)-(a))+((d)-(c))*((d)-(c)))
#define TEMPFACTOR 0.90

int metrop(float de, float t);
float revcst(float x[], float y[], int iorder[], int N, int n[]);
void reverse(int iorder[], int N, int n[]);
/*reverse move*/
float trncst(float x[], float y[], int iorder[], int N, int n[]);
void trnspt(int iorder[], int N, int n[]);
/*transport move*/

void anneal(float x[], float y[], int iorder[], int N){
/* This algorithm finds the shortest round-trip path to N cities
   whose coordinates are in the arrays x[1..N],y[1..N]. The
   array iorder[1..N] specifies the order in which the cities are
   visited.  On input, the   elements of iorder may be set to any
   permutation of the numbers 1 to N.  This routine will return
   the best alternative path it can find.  */
    int ans, L, sizefactor, initprob, alpha, minpercent, changes, i1, i2, tcount, fcount, trials, FLIMIT;
    int counter = 0;
    int i, j, k, nsucc, nn, idec;
    static int n[7];
    float path, de, t;
    t = 0.1;     /*temperature */
    tcount = 1;  /* Counter for temperatures*/
    trials = 0;  /* Counter for temperatures*/
    fcount = 0;  /* Counter for frozen*/
    changes = 0; /* Counter for moves*/
    FLIMIT = 5;  /* Limit for frozen*/
    sizefactor = 100;
    L = sizefactor * N; /* max number of paths at any temp */
    initprob = 0.40;    /* initial probability */
    minpercent = 0.01;  /* Breakout minpercent */
    path = 0.0;         /*initial path length*/
    alpha = 0.05;       /*alpha cutoff percentage */

    printf("Initial order:\n");
    /*the following for loop calculates the distance of the TSP */
    for (i=1; i <= N; i++)
        printf("city %d \n", iorder[i]);


    for (i=1; i<N; i++) {
        /* calculate initial path length */
        i1 = iorder[i];
        /*Determine position in the permuted solution space.*/
        i2 = iorder[i+1];
        path += ALEN(x[i1], x[i2], y[i1], y[i2]);
    }
    i1 = iorder[N];
    i2 = iorder[1];
    path += ALEN(x[i1], x[i2], y[i1], y[i2]);
    printf("path length %12.6f \n", path);
    /*print path length*/

    do {
        changes = trials = 0;
        do {
            trials = trials + 1;
            /* try up to 100 temperature steps */
            nsucc = 0;
            for (k = 1; k <= L; k++) {
                do {
                    n[1] = 1 + (int) (N * Random());
                    /*n[1] is the start of the segment*/
                    n[2] = 1 + (int) ((N - 1) * Random());
                    /*n[2} is the end of the segment*/
                    if (n[2] >= n[1]) ++n[2];
                    nn = 1 + ((n[1] - n[2] + N - 1) % N);
                    /*nn is the number of cities not in the segment*/
                }
                while (nn < 3);
                idec = Equilikely(0, 1);
                /* decide between segment reversal and transport */
                if (idec == 0) {
                    /* do a transport */
                    n[3] = n[2] + (int) (abs(nn - 2) * Random()) + 1;
                    n[3] = 1 + ((n[3] - 1) % N);
                    /* transport to a location NOT on the path */
                    de = trncst(x, y, iorder, N, n);
                    /* find cost */
                    if (de <= 0.0) {
                        /* a downhill move */
                        changes = changes + 1;
                        counter = 0;
                        /* reset counter to 0 */
                    }
                    else {
                        /*an uphill move */
                        if (0.5 >= exp(-de / t)) {
                            /*allow if true */
                            changes = changes + 1;
                            counter = 0;
                        }
                    }
                    ans = metrop(de, t);
                    if (ans) {
                        ++nsucc;
                        path += de;
                        trnspt(iorder, N, n);
                        /* transport */
                    }
                }
                else {
                    /* do a path reversal */
                    de = revcst(x, y, iorder, N, n);
                    /* find cost, i.e de = change in objective function */
                    if (de <= 0.0) {
                        /* a downhill move */
                        changes = changes + 1;
                        counter = 0;
                        /* reset counter to 0 */
                    }
                    else {
                        /*an uphill move */
                        if (0.5 >= exp(-de / t)) {
                            /*allow if true */
                            changes = changes + 1;
                            counter = 0;
                        }
                    }
                    ans = metrop(de, t);
                    if (ans) {
                        ++nsucc;
                        path += de;
                        reverse(iorder, N, n);
                        /* reversal */
                    }
                }
                if ((changes / trials) <= minpercent) {
                    fcount = fcount + 1;
                    printf("Fcount: %i\n", fcount);
                }
                /*break out of a higher temperature once a number of successes is reached
            (i.e. we finish early if we have enough successes */
            }
            printf("\n %s %10.6f :%s %12.6f \n", "T =", t, "    Path Length =", path);
            printf("Successful Moves: %6d\n", nsucc);
            printf("T Count: %d", tcount);
            t *= TEMPFACTOR;
            fcount = 0;
            tcount = tcount + 1;
            if (nsucc == 0) return;
        } while ((trials < L) && (changes < (alpha*N)));
    } while (fcount < FLIMIT);
}


void reverse(int iorder[], int N, int n[]) {
    int nn, j, k;
    int m;
    int itmp;
    nn = (1 + ((n[2]-n[1]+N) % N))/2;
    for (j=1; j<=nn; j++) {
        k = 1 + ((n[1]+j-2) % N);
        m = 1 + ((n[2] - j + N) % N);
        itmp = iorder[k];
        iorder[k] = iorder[m];
        iorder[m] = itmp;
    }
}

void trnspt(int iorder[], int N, int n[]) {
    int m1, m2, m3, nn, j, jj, *jorder;
    jorder = ivector(1, N);
    m1 = 1 + ((n[2] - n[1] + N) % N);
    m2 = 1 + ((n[5] - n[4] + N) % N);
    m3 = 1 + ((n[3] - n[6] + N) % N);
    nn = 1;
    for (j=1; j <= m1; j++) {
        jj = 1 + ((j + n[1] - 2) % N);
        jorder[nn++] = iorder[jj];
    }
    for (j=1; j <= m2; j++) {
        jj = 1 + ((j + n[4] - 2) % N);
        jorder[nn++] = iorder[jj];
    }
    for (j=1; j <= m3; j++) {
        jj = 1 + ((j + n[6] - 2) % N);
        jorder[nn++] = iorder[jj];
    }
    for (j=1; j <= N; j++)
        iorder[j] = jorder[j];
    free_ivector(jorder, 1, N);
}


float trncst(float x[], float y[], int iorder[], int N, int n[]){
    float xx[7], yy[7], de;
    int j, ii;
    n[4] = 1 + (n[3] % N);
    n[5] = 1 + ((n[1] + N - 2) % N);
    n[6] = 1 + (n[2] % N);
    for (j = 1; j <= 6; j++){
        ii = iorder[n[j]];
        xx[j] = x[ii];
        yy[j] = y[ii];
    }
    de = -ALEN(xx[2], xx[6], yy[2], yy[6]);
    de -= ALEN(xx[1], xx[5], yy[1], yy[5]);
    de -= ALEN(xx[3], xx[4], yy[3], yy[4]);
    de += ALEN(xx[1], xx[3], yy[1], yy[3]);
    de += ALEN(xx[2], xx[4], yy[2], yy[4]);
    de += ALEN(xx[5], xx[6], yy[5], yy[6]);
    return de;
}

float revcst(float x[], float y[], int iorder[], int N, int n[]) {
    float xx[5], yy[5], de;
    int j, ii;
    n[3] = 1 + ((n[1] + N - 2) % N);
    n[4] = 1 + (n[2] % N);
    for (j = 1; j <= 4; j++){
        ii = iorder[n[j]];
        xx[j] = x[ii];
        yy[j] = y[ii];
    }
    de = -ALEN(xx[1], xx[3], yy[1], yy[3]);
    de -= ALEN(xx[2], xx[4], yy[2], yy[4]);
    de += ALEN(xx[1], xx[4], yy[1], yy[4]);
    de += ALEN(xx[2], xx[3], yy[2], yy[3]);
    return de;
}

int metrop(float de, float t) {
    return de < 0.0 || Random() < exp(-de/t);
}

#define CITIES 194

int main (void) {
    int i, j;
    int N = CITIES;
    int iorder[CITIES+1], node[CITIES+1];
    float x[CITIES+1], y[CITIES+1];
    FILE *fp;

    long SEED= 13213.45;
    /*keeping the same seed value allows you to duplicate an experiment.
    The seed value can be anything. */
    PlantSeeds(SEED);

    for (i=1; i <= N; i++)
        iorder[i] = i;
    /*coordinates for original cities, to read in another file,
    use an 'io' (input output) routine to read in a data file*/

    for(i=1; i<= N; i++)
        iorder[i] = N -i +1;

    fp = fopen("qatar2.dat", "r");

    if(fp == NULL) {
        printf("\nError: Unable to open	file.\n");
        return;
    }

    for(j = 1; j<= N; j++){
        if(fscanf(fp, "%i%f%f", &node[j],&x[j], &y[j]) != 3){
            printf("\nError: Incorrect input.\n");
            return;
        }
    }
    anneal(x,y,iorder,N);

    printf("\nOrder after anneal:\n");
    for (i=1; i <= N ; i++)
        printf("city %d \n", iorder[i]);
}

