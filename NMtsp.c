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
#define TFACTR 0.9

int metrop(float de, float t);

float revcst(float x[], float y[], int iorder[], int ncity, int n[]);
void reverse(int iorder[], int ncity, int n[]);
/*reverse move*/
float trncst(float x[], float y[], int iorder[], int ncity, int n[]);
void trnspt(int iorder[], int ncity, int n[]);
/*transport move*/

void anneal(float x[], float y[], int iorder[], int ncity){
/* This algorithm finds the shortest round-trip path to ncity cities
   whose coordinates are in the arrays x[1..ncity],y[1..ncity]. The
   array iorder[1..ncity] specifies the order in which the cities are
   visited.  On input, the   elements of iorder may be set to any
   permutation of the numbers 1 to ncity.  This routine will return
   the best alternative path it can find.  */
    int ans, nover, nlimit, i1, i2, count;
    int counter = 0;
    int i, j, k, nsucc, nn, idec;
    static int n[7];
    float path, de, t;
    t = 0.1;
    count = 1;
    /*temperature */
    nover = 100 * ncity;
    /* max number of paths at any temp */
    nlimit = 10 * ncity;
    /* max number of successful path changes */
    path = 0.0;
    /*initial path length*/
    printf("Initial order:\n");
    /*the following for loop calculates the distance of the TSP */
    for (i=1; i <= ncity; i++)
        printf("city %d \n", iorder[i]);
    for (i=1; i<ncity; i++) {
    /* calculate initial path length */
        i1 = iorder[i];
        /*Determine position in the permuted solution space.*/
        i2 = iorder[i+1];
        path += ALEN(x[i1], x[i2], y[i1], y[i2]);
    }
    i1 = iorder[ncity];
    i2 = iorder[1];
    path += ALEN(x[i1], x[i2], y[i1], y[i2]);
    printf("path length %12.6f \n", path);
    /*print path length*/
    for (j=1; j <= 100; j++,count++) {
    /* try up to 100 temperature steps */
        nsucc = 0;
        for (k = 1; k <= nover; k++) {
            do {
                n[1] = 1 + (int) (ncity * Random());
                /*n[1] is the start of the segment*/
                n[2] = 1 + (int) ((ncity - 1)*Random());
                /*n[2} is the end of the segment*/
                if (n[2] >= n[1]) ++n[2];
                nn = 1 + ((n[1] - n[2] + ncity - 1) % ncity);
                /*nn is the number of cities not in the segment*/
            }
            while (nn < 3);
            idec = Equilikely(0,1);
            /* decide between segment reversal and transport */
            if (idec == 0) {
            /* do a transport */
                n[3] = n[2] + (int) (abs(nn-2) * Random()) + 1;
                n[3] = 1 + (( n[3] - 1) % ncity);
                /* transport to a location NOT on the path */
                de = trncst (x,y,iorder,ncity,n);
                /* find cost */
                if (de < 0.0)
                /* a downhill move */
                    counter = 0;
                    /* reset counter to 0 */
                ans = metrop(de,t);
                if (ans) {
                    ++nsucc;
                    path += de;
                    trnspt(iorder,ncity,n);
                    /* transport */
                }
            }
            else {
            /* do a path reversal */
                de = revcst(x,y,iorder,ncity,n);
                /* find cost, i.e de = change in objective function */
                if (de < 0.0)
                /* a downhill move!!! */
                    counter = 0;
                    /* reset counter to 0 */
                ans = metrop(de,t);
                if (ans) {
                    ++nsucc;
                    path += de;
                    reverse(iorder,ncity,n); /* reversal */
                }
            }
            if ( nsucc >= nlimit ) break;
            /*break out of a higher temperature once a number of successes is reached
            (i.e. we finish early if we have enough successes */
        }
        printf("\n %s %10.6f :%s %12.6f \n","T =",t, "    Path Length =",path);
        printf("Successful Moves: %6d\n",nsucc);
        printf("Count: %d", count);
        t *= TFACTR;
        if (nsucc == 0) return;
    }
}

void reverse(int iorder[], int ncity, int n[]) {
    int nn, j, k;
    int m;
    int itmp;
    nn = (1 + ((n[2]-n[1]+ncity) % ncity))/2;
    for (j=1; j<=nn; j++) {
        k = 1 + ((n[1]+j-2) % ncity);
        m = 1 + ((n[2] - j + ncity) % ncity);
        itmp = iorder[k];
        iorder[k] = iorder[m];
        iorder[m] = itmp;
    }
}

void trnspt(int iorder[], int ncity, int n[]) {
    int m1, m2, m3, nn, j, jj, *jorder;
    jorder = ivector(1, ncity);
    m1 = 1 + ((n[2] - n[1] + ncity) % ncity);
    m2 = 1 + ((n[5] - n[4] + ncity) % ncity);
    m3 = 1 + ((n[3] - n[6] + ncity) % ncity);
    nn = 1;
    for (j=1; j <= m1; j++) {
        jj = 1 + ((j + n[1] - 2) % ncity);
        jorder[nn++] = iorder[jj];
    }
    for (j=1; j <= m2; j++) {
        jj = 1 + ((j + n[4] - 2) % ncity);
        jorder[nn++] = iorder[jj];
    }
    for (j=1; j <= m3; j++) {
        jj = 1 + ((j + n[6] - 2) % ncity);
        jorder[nn++] = iorder[jj];
    }
    for (j=1; j <= ncity; j++)
        iorder[j] = jorder[j];
    free_ivector(jorder, 1, ncity);
}


float trncst(float x[], float y[], int iorder[], int ncity, int n[]){
    float xx[7], yy[7], de;
    int j, ii;
    n[4] = 1 + (n[3] % ncity);
    n[5] = 1 + ((n[1] + ncity - 2) % ncity);
    n[6] = 1 + (n[2] % ncity);
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

float revcst(float x[], float y[], int iorder[], int ncity, int n[]) {
    float xx[5], yy[5], de;
    int j, ii;
    n[3] = 1 + ((n[1] + ncity - 2) % ncity);
    n[4] = 1 + (n[2] % ncity);
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

#define CITIES 10

int main (void) {
    int i;
    int ncity = CITIES;
    int iorder[CITIES+1];
    float x[CITIES+1], y[CITIES+1];
    long SEED= 31764.0;
    /*keeping the same seed value allows you to duplicate an experiment.
    The seed value can be anything. */
    PlantSeeds(SEED);

    for (i=1; i <= ncity; i++)
        iorder[i] = i;
        /*coordinates for original cities, to read in another file,
        use an 'io' (input output) routine to read in a data file*/
    x[1]=1;
    y[1]=3;
    x[2]=2;
    y[2]=2;
    x[3]=5;
    y[3]=9;
    x[4]=8;
    y[4]=8;
    x[5]=2;
    y[5]=7;
    x[6]=7;
    y[6]=1;
    x[7]=5;
    y[7]=5;
    x[8]=6;
    y[8]=4;
    x[9]=3;
    y[9]=4;
    x[10]=1;
    y[10]=5;
    anneal(x,y,iorder,ncity);

    printf("\nOrder after anneal:\n");
    for (i=1; i <= ncity ; i++)
        printf("city %d \n", iorder[i]);
}

