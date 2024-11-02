/*
This file is supplemental material to the paper
"Attempting the impossible: enumerating extremal submodular functions for n=6"
by E P Csirmaz and L Csirmaz.

Copyright 2024 E P Csirmaz and L Csirmaz

This program is free software: you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by the Free
Software Foundation, either version 3 of the License, or (at your option)
any later version.

This program is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
more details.

You should have received a copy of the GNU General Public License along with
this program.  If not, see <https://www.gnu.org/licenses/>.
*/

/* read a file of an iteration result; combine all positive/negative
   rays and check if they are extremal rays of the final cone.
   usage: extract <iteration> <axno> <input> <output>
*/

#define T_FACTOR    static int

#if RANK==3
#  include "ax3.c"
#  include "sym3.c"
#elif RANK==4
#  include "ax4.c"
#  include "sym4.c"
#elif RANK==5
#  include "ax5.c"
#  include "sym5.c"
# elif RANK==6
#  include "ax6.c"
#  include "sym6.c"
# else
#  error RANK should be defined to 3,4,5,6
#endif

#define    DIM    VARS

#define stringify(x)    #x
#define mkstringof(x)   stringify(x)


/* =================================================================== */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static void usage(void)
{
    printf(
"usage: extr" mkstringof(RANK) " [flags] <n> <newax> <input> <output> [<skipfile>]\n"
"extracting valid extremal rays from an incomplete iteration\n"
"  --ray    save the ray, not the orbit\n"
"  --v      verbose\n"
"  --x      do not check extremity of the combined ray\n"
"  n        iteration count (also in input file)\n"
"  newax    axiom number to be probed\n"
"  input    last generation of extremal rays\n"
"  output   file to save the result\n"
"  skipfile stop processing if exists\n");
    exit(1);
}


static int verbose;          /* verbose setting */
static int checkrank;        /* 0: check; 1: no */
static int noorbit;          /* 1 if no minimal version */
static int iteration;        /* iteration we are in */
static int axiomno[AXIOMS];  /* axioms handled so far */
static int newax;            /* the new axiom to be added */
static FILE *ifile;          /* input file */
static FILE *ofile;          /* output file */
static const char *skipfile; /* the skipfile */


static void handle_params(int argc, const char*argv[])
{
    const char *w, *input,*output;
    verbose = 0; 
    checkrank = 0; 
    noorbit = 0;
    while(argc > 1 && argv[1][0] == '-') {
       if(strcmp(argv[1],"--ray") == 0) {
           argc--; argv++; noorbit=1;
       } else if(strcmp(argv[1],"--v") == 0) {
           argc--; argv++; verbose=1;
       } else if(strcmp(argv[1],"--x") == 0) {
           argc--; argv++; checkrank=1;
       } else {
           printf("Unknown flag %s\n",argv[1]); 
           usage();
       }
    }
    if(argc != 5 && argc != 6) usage();
    skipfile = argc == 6 ? argv[5] : "skip";
    w = argv[1]; 
    iteration = 0;
    while('0'<= *w && *w <= '9'){ iteration = iteration * 10 + (*w - '0'); w++; }
    if(*w || iteration<DIM || iteration >= AXIOMS) {
       printf("iteration number %d out of range\n", iteration); 
       exit(1);
    }
    newax = 0; 
    w = argv[2];
    while('0' <= *w && *w <= '9'){ newax = newax * 10 + (*w - '0'); w++; }
    if(*w || argv[2][0] == 0 || newax >= AXIOMS) {
       printf("wrong new axiom number %d\n", newax); 
       exit(1);
    }
    input = argv[3];
    output = argv[4];
    if(!input || !*input || (ifile = fopen(input, "r")) == NULL) {
        printf("input file %s not found\n", input); 
        exit(1);
    }
    if(!output || !*output || (ofile = fopen(output, "a")) == NULL) {
        printf("Cannot append to output file %s\n", output); 
        exit(1);
    }
}


/* ==========================================================================*/

#define EPS    5e-9

static inline int is_zero(double x) {
     return (-EPS < (x) && (x) < EPS);
}


static inline int not_zero(double x) {
    return ((x) <= -EPS || (x) >= EPS);
}


static int iinner(int axno,const int *ray)
{
    int v=0; 
    const int *ax=axioms[axno];
    for(int i = 0; i < DIM; i++) v += ax[i] * ray[i];
    return v;
}


/* =================================================================== */
static int posrays=0, negrays=0, zerorays=0;    // number of rays
static int nextray[DIM];             // reading the  next ray
static int newrays=0;                // rays in the output


/* read the input file first; compute the number of positive, negative rays
*  relative to the next axiom */
static void first_input_scan(void)
{
    int w;
    // iteration number
    w=0;
    if(fscanf(ifile, "%d", &w) != 1 || w != iteration) {
        printf("wrong iteration number %d in input file\n", w); 
        exit(1);
    }
    // iteration many axiom numbers handled so far
    for(int i=0; i < iteration; i++) {
        w=0;
        if(fscanf(ifile, "%d", &w) != 1 || w < 0 || w >= AXIOMS) {
            printf("axiom list in input file corrupted (i=%d, w=%d)\n", i, w); 
            exit(1);
        }
        axiomno[i] = w;
    }
    while(1) {  // read rays
        for(int i = 0; i < DIM; i++) {
            if(fscanf(ifile, "%d", &w) != 1) {
                if(i == 0 && feof(ifile)) { // last item read
                    if(posrays + negrays + zerorays < DIM) {
                        printf("Not enough data in input file\n"); exit(1);
                    }
                    if(posrays <= 0 || negrays <= 0) {
                        printf("consistency error: rays %d/%d (positive/negative)\n",
                            posrays, negrays); 
                        exit(1);
                    }
                    rewind(ifile); // prepare for re-reading
                    return; 
                }
                printf("input file is corrupt\n"); 
                exit(1);
            }
            nextray[i] = w;
            if(w < -5000 || w > 5000) {
                printf("large value in input file: %d\n", w); 
                exit(1);
            }
        }
        w = iinner(newax, nextray);
        if(w < 0) { negrays++; }
        else if(w > 0) { posrays++; }
        else { zerorays++; }
    }
}


/* ===================================================================== */
static void save_nextray(int ridx); // forward declaration
/* ===================================================================== */


/* second scan; output file is open */
#include "bitmap.h"
static int *pos_store, *neg_store; // storing positive and negative rays

#define STSIZE    (DIM+1)

static void second_input_scan(void)
{
    int w;
    int pidx, nidx, zidx; 
    int rayidx;
    /* allocate memory for positive and negative rays */
    neg_store = (int *)malloc(sizeof(int) * STSIZE * (negrays + posrays));
    if(neg_store == NULL) {
        printf("Cannot allocate memory for rays %d / %d\n", posrays, negrays);
        exit(1);
    }
    pos_store = neg_store + (STSIZE * negrays);
    if(init_bitmaps(iteration + 1, posrays + negrays + zerorays) == 0) {
        printf("Cannot allocate memory for bitmaps: %d / %d / %d\n", posrays, negrays, zerorays);
        exit(1);
    }
    /* go over the input file again */
    if(fscanf(ifile, "%d", &w) != 1 || w != iteration) {
        printf("input file consistency error #1\n"); 
        exit(1);
    }
    for(int i = 0; i < iteration; i++) if(fscanf(ifile, "%d", &w) != 1 || axiomno[i] != w) {
        printf("input file consistency error #2\n"); 
        exit(1);
    }
    pidx=0; 
    nidx=0; 
    zidx=0;
    while(1) {  // read rays
        for(int i = 0; i < DIM; i++) {
            if(fscanf(ifile, "%d", &w) != 1) {
                if(i == 0 && feof(ifile)) { // sanity check
                    if(pidx != posrays || nidx != negrays || zidx != zerorays) {
                        printf("second scan inconsistency\n"); 
                        exit(1);
                    }
                    fclose(ifile);
                    return;
                }
                printf("input file consistency error #3\n"); 
                exit(1);
            }
            nextray[i] = w;
        }
        w = iinner(newax, nextray);
        if(w < 0) { // negative
            memcpy(neg_store + STSIZE * nidx, nextray, sizeof(int) * DIM);
            neg_store[STSIZE * nidx + DIM] = w;
            rayidx = nidx;
            nidx++;
        } else if(w > 0) { // positive
            memcpy(pos_store + STSIZE * pidx, nextray, sizeof(int) * DIM);
            pos_store[STSIZE * pidx + DIM] = w;
            rayidx = pidx + negrays;
            pidx++;
        } else { // zero
            rayidx = zidx + negrays + posrays;
            zidx++;
        }
        if(w >= 0) save_nextray(-1); // do not check
        // check if the ray is adjacent to the first axioms
        for(int ax = 0; ax < iteration; ax++) {
             w = iinner(axiomno[ax], nextray);
             if(w < 0){ printf("ray is on the negative side of axiom %d\n", axiomno[ax]); exit(1); }
             if(w == 0) // it is adjacent
                 add_bitmap(ax, rayidx);
        }
    }
}


/* =============================================================================== */
/* combine a negative and a positive ray */

static int gcd(int a, int b) 
{
    if(a == 1) return 1;
    if(a == 0) return b < 0 ? -b : b;
    if(a < 0) { a = -a; } 
    if(b < 0) { b = -b; }
again:
    b %= a; if(b == 0) return a; if(b == 1) return 1;
    a %= b; if(a == 0) return b; if(a == 1) return 1;
    goto again;
}


static void simplify_nextray(void)
{
    int v;
    v = nextray[0]; 
    if(v < 0) v = -v;
    for(int i = 1; v != 1 && i < DIM; i++) v = gcd(v, nextray[i]);
    if(v <= 0){ printf("nextray is all zero (%d)\n", v); exit(1); }
    if(v == 1) return;
    for(int i = 0; i < DIM; i++) nextray[i] /= v;
}


static void combine_rays(int posidx,int negidx)
{
    const int *r1 = pos_store + STSIZE * posidx;
    const int *r2 = neg_store + STSIZE * negidx;
    int r1d = r1[DIM], r2d = -r2[DIM];
    if(r1d <= 0 || r2d <= 0){ printf("wrong sign\n"); exit(1); }
    for(int i = 0; i < DIM; i++) {
        // if it has a negative coeff, then not extremal
        if((nextray[i] = r2d * r1[i] + r1d * r2[i]) < 0) return;
    }
    simplify_nextray();
    save_nextray(checkrank ? -2 : negrays+posidx);
}


/* ===================================================================== */
/* timer */
#include <sys/time.h>
static struct timeval prev_time;

static float elapsed_time(int set)
{
    struct timeval current_time;
    if(set){ gettimeofday(&prev_time, NULL); return 0.0; }
    gettimeofday(&current_time, NULL);
    return (current_time.tv_sec - prev_time.tv_sec) + 
        (current_time.tv_usec - prev_time.tv_usec) / 1000.0 / 1000.0;
}


static void printtime(float t)
{
    int d, h, m, s;
    if(t < 59.994){ printf("%.2f", t); return; }
    s = (int)(t + 0.49);
    m = s / 60; 
    s = s % 60;
    if(m < 60) { printf("%d:%02d", m, s); return; }
    h = m / 60; 
    m = m % 60;
    if(h < 24){ printf("%d:%02d:%02d", h, m, s); return; }
    d = h / 24; 
    h = h % 24;
    printf("%dd %d:%02d:%02d", d, h, m, s);
}



/* =============================================================== */
/* statistics */
static int negidx; 
static size_t combchkyes;

static void statistics(void)
{
    static float lastreport = 0.0; 
    float elapsed;
    elapsed = elapsed_time(0);
    if(elapsed - lastreport < 30.0) return;
    lastreport = elapsed;
    printf("Elapsed "); 
    printtime(elapsed);
    printf(", rays %d, checked %4.2f%% (%zu), idx=%d\n",
       newrays,
       100.0 * ((float)combchkyes) / ((float)posrays) / ((float)negrays),
       combchkyes, negidx);
}


#include <sys/stat.h>
/* return 1 if file "skip" exists */
static int check_skipfile(void)
{
    struct stat sb;
    return !skipfile ? 0 : lstat(skipfile, &sb) < 0 ? 0 : 1;
}


static void find_adjacent_rays(void)
{
    int pairs_checked = 0;
    combchkyes = 0;
    for(negidx = 0; negidx < negrays; negidx++) {
        if(check_skipfile()) {
            statistics();
            printf("Program run aborted by skipfile\n");
            exit(1);
        }
        if(pairs_checked >= 5000000) {
            pairs_checked = 0; 
            if(verbose) statistics();
        }
        prepare_adjacency(negidx, DIM - 2);
        for(int pidx = 0; pidx < posrays; pidx++) {
            combchkyes++;
            pairs_checked++;
            if(!fast_raycheck(negrays+pidx)) 
                combine_rays(pidx, negidx);
        }
    }
}


/* =============================================================================== */
/* the ray to be handled is in nextray[] */
/* find the minimal permutation of nextray[] */


static int minperm[DIM];    // the minimal permutation found so far 
static int dual[DIM];        // the dual of intray[]
static int nextperm[DIM];    // next permutation is generated here

static void fill_dual(void)    // fill dual[]
{   
    for(int i = 0; i < DIM; i++) {
        dual[i] = 0;
        for(int j = 0; j < DIM; j++){ dual[i] += nextray[j] * dualmatrix[i][j]; }
    }
}


static int cmp_as_string(int a1,int b1)
{   
    if(a1 == b1) return 0; 
    if(a1 < b1) return +1; 
    return -1; 
}


/* if nextperm[] is lexicographycally smaller than mimperm, replace */
static void lexminperm(void)
{
    int smaller = 0;
    for(int i = 0; smaller == 0 && i < DIM; i++) {
        smaller = cmp_as_string(nextperm[i], minperm[i]);
    }
    if(smaller > 0) {
       for(int i = 0; i < DIM; i++) minperm[i] = nextperm[i];
    }
}


/* copy the final result to minperm[] */
static void check_minperm(void)
{   
    // fill the minimal value
    for(int i = 0; i < DIM; i++){ minperm[i] = nextray[i]; }
    if(noorbit) return;
    fill_dual(); // and compute the dual
    for(int pm = 0; pm < VPERMNO; pm++) {
       for(int i = 0; i < DIM; i++) nextperm[i] = nextray[VPERMS[pm][i]];
       lexminperm();
       for(int i = 0; i < DIM; i++) nextperm[i] = dual[VPERMS[pm][i]];
       lexminperm();
    }
}


/* ===================================================================== */
// how many axioms are satisfied by minperm[]
static int axn_given=0;            // weight of the ray
static int handled_axioms[AXIOMS];    // axioms which are zero

static int find_axrank(void)
{  
    axn_given = 0;
    for(int axn = 0; axn < AXIOMS; axn++) {
        int w = iinner(axn, nextray);
        if(w < 0) return 1;
        if(w == 0) { handled_axioms[axn_given] = axn; axn_given++; }
   }
   if(axn_given < DIM - 1) { printf("axn=%d < DIM-1\n", axn_given); exit(1); }
   return 0;
}


#ifdef ALGEBRAIC_TEST
// axiom numbers this ray is on are in handled_axioms[]
// the rank must be exactly DIM-1
// axn_given is at least DIM-1
static double M[AXIOMS][DIM];

static int check_extremity(void)
{
    int rank=0;
    // fill the matrix M
    for(int j = 0; j < axn_given; j++)
        for(int i = 0; i < DIM; i++) M[j][i] = axioms[handled_axioms[j]][i];
    // go over all rows
    for(int j = 0; j < axn_given; j++) { 
        double *Mj, maxM; 
        int pivot;
        // find the maximal absolut value in row j
        Mj = &(M[j][0]); 
        maxM = Mj[0]; 
        if(maxM < 0.0) { maxM = -maxM; } 
        pivot=0;
        for(int i = 1; i < DIM; i++){ 
            double w = Mj[i];
            if(w < 0.0) {
                if(maxM + w < 0.0) { maxM = -w; pivot = i; }
            } else if(maxM < w) { maxM = w; pivot = i; }
        }
        if(maxM < EPS) continue; // all zero line
        // the maximal value is in column pivot
        rank++;
        maxM = 1.0 / Mj[pivot];
        for(int i = 0; i < DIM; i++) Mj[i] *= maxM;
        Mj[pivot]=1.0;
        // subtract row j from all subsequent rows making column pivot zero
        for(int jj = j + 1; jj < axn_given; jj++) {
            double w = M[jj][pivot];
            if(not_zero(w)) {  // subtract w*Mj[i]
                for(int i = 0; i < DIM; i++) { M[jj][i] -= w*Mj[i]; }
            }
            M[jj][pivot] = 0.0;
        }
    }
    // now we must have rank==DIM-1
    if(rank == DIM - 1) return 0;
    return 1;
}
#endif
/* =============================================================================== */


// do not touch nextray[]
static void save_nextray(int ridx)
{
    if(find_axrank()) return; // not all axioms are >=0 
    if(ridx >= 0 && !are_adjacent_rays(ridx)) return; // not extremal
    check_minperm();    // the result is in minperm[]
    // and print it
    newrays++;
    fprintf(ofile, "%d: %d", axn_given, minperm[0]);
    for(int i = 1; i < DIM; i++) fprintf(ofile, ",%d", minperm[i]);
    fprintf(ofile, "\n");
}


int main(int argc, const char*argv[])
{
    handle_params(argc, argv);
    if(check_skipfile()){ printf("Skipfile %s exists, aborting\n", skipfile); return 1; }
    elapsed_time(1);        // set time
    first_input_scan();     // find positive/negative rays
    printf("Extracting rays, iteration=%d, newax=%d, pos/zero/neg=%d/%d/%d\n",
        iteration, newax, posrays, zerorays, negrays);
    second_input_scan();    // closes ifile
    find_adjacent_rays();
    fclose(ofile);
    printf("Done, rays=%d, idx=%d, time=", newrays, negidx); 
    printtime(elapsed_time(0)); 
    printf("\n");
    return 0;
}

/* EOF */

