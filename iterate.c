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

/* This program act as a pipe while iterating the DD method.
*  it reads a DD pair, the next inequality to be added, and outputs the
*  next DD pair.
*
*  Usage: iterate <n> <axno> <input> <output> <stopfile>
*    n      iteration count, inequalities handled so far, also in the input file
*    axno   inequality number to be added next
*    input  input file, see description below
*    output output file
*    stopfile abort running if the file exists
*
*  input file format:
*    only decimal numbers separated by white spaces
*    <iteration count>
*    <axiom numbers> iteration many axioms handled so far 
*    <ray1> <ray2> ... <rayn>
*    each ray is a list of DIM many (possibly signed) integers
*  We assume that no coefficient has absolute value above 1700,
*  so integer arithmetic can be used. This condition is checked.
*
*  The created output file can be used for input for the next iteration.
*
*  The input file is read twice. First the syntax is checked; each
*  ray should be on the >=0 size of the first ITERATION axioms. The
*  number of positive, negative and zero rays are determined.
*  Ray indices are sorted as positive, negative, zero; zero
*  rays are not stored, only passed to the output file (along with
*  the positive rays). Next, for each positive/negative ray pair the
*  adjacency test is performed; if they are adjacent, their linear
*  combination is created and reduced by dividing by the gcd of all
*  coeffs.
*/

#define T_FACTOR    static int
#if RANK==3
#  include "ax3.c"
#elif RANK==4
#  include "ax4.c"
#elif RANK==5
#  include "ax5.c"
#elif RANK==6
#  include "ax6.c"
#else
#  error RANK should be defined to 3,4,5,6
#endif

#define DIM    VARS

#define stringify(x)    #x
#define mkstringof(x)    stringify(x)

/* ============================================================================ */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static int verbose = 0;        // print additional info
static int iteration = 0;      // iteration we are in
static int axiomno[AXIOMS];    // axioms handled so far
static FILE *ifile, *ofile;    // input and output file
static const char *skipfile;   // the skipfile


static void usage(void) {
    printf(
"usage: iter" mkstringof(RANK) " [-v] <n> <newax< <input> <output> [<skipfile>]\n"
"cutting the given cone with the specified axiom for rank=" mkstringof(RANK) "\n"
"  -v       (verbose) print additional info\n"
"  n        iteration count, also in the input file\n"
"  newax    axiom number to be added next\n"
"  input    input file containing the last generation of extremal rays\n"
"  output   output file for the next generation of extremal rays\n"
"  skipfile abort the program if exists; default is \"skip\"\n"
"The output file should not exist.\n");
    exit(1);
}


static void handle_params(int argc, const char *argv[]) {
    const char *w, *input, *output; 
    int n;
    
    while(argc > 1 && argv[1][0] == '-') {
        if(argv[1][1] == 'v') { verbose = 1; }
        else usage();
        argc--; 
        argv++;
    }
    if(argc != 5 && argc != 6) usage();
    skipfile = argc == 6 ? argv[5] : "skip";
    w = argv[1]; 
    iteration = 0;
    while(*w && '0' <= *w && *w <= '9') { iteration = iteration * 10 + (*w - '0'); w++; }
    if(*w || iteration < DIM || iteration >= AXIOMS) {
        printf("Iteration number %d is out of range\n", iteration); 
        exit(1);
    }
    n = 0; 
    w = argv[2];
    while(*w && '0' <= *w && *w <= '9'){ n = n * 10 + (*w - '0'); w++; }
    if(*w || n >= AXIOMS){ printf("wrong axiom number %d\n", n); exit(1); }
    // 0 .. iteration-1 : previous iterations; these axoms are added
    axiomno[iteration] = n;
    input = argv[3]; 
    output = argv[4];
    if(!input || !*input || (ifile = fopen(input, "r")) == NULL) {
        printf("input file \"%s\" not found\n", input);
        exit(1);
    }
    if(!output || !*output || (ofile = fopen(output, "r")) != NULL) {
        printf("output file \"%s\" exists, cannot continue\n", output);
        exit(1);
    }
    ofile = fopen(output, "w");
    if(ofile == NULL) {
        printf("Cannot create oputput file \"%s\"\n", output);
        exit(1);
    }
}


/* =================================================================== */
static int posrays = 0, negrays = 0, zerorays = 0; // number of rays
static int nextray[DIM];            // reading the next ray
static int newrays = 0;                // rays in the output


/* integer dot product */
static int intdot(const int*a, const int*b) {
    int sum = 0;
    for(int i = 0; i < DIM; i++) sum += a[i] * b[i];
    return sum;
}


/* read the input file first; compute the number of positive, negative rays
*  relative to the next axiom */
static void first_input_scan(void) {
    int w;
    // iteration number
    w = 0;
    if(fscanf(ifile, "%d", &w) != 1 || w != iteration) {
        printf("wrong iteration number %d in input file\n", w); 
        exit(1);
    }
    // iteration many axiom numbers handled so far
    for(int i = 0; i < iteration; i++) {
         w = 0;
         if(fscanf(ifile, "%d", &w) != 1 || w < 0 || w >= AXIOMS) {
             printf("axiom list in input file corrupted (i=%d, w=%d)\n", i, w); 
             exit(1);
         }
         axiomno[i] = w;
    }
    // axiomno[0 .. iteration] must be all different
    for(int i = 0; i <= iteration; i++) {
        w = axiomno[i];
        for(int j = 0; j < i; j++)
            if(axiomno[j] == w) {
                printf("axiomno %d is repeated (%d,%d)\n", w, i, j); 
                exit(1);
            }
    }    
    while(1) {  // read rays
        for(int i = 0; i < DIM; i++) {
            if(fscanf(ifile, "%d", &w) != 1) {
                if(i == 0 && feof(ifile)) { // last item read
                    if(posrays + negrays + zerorays < DIM) {
                        printf("Not enough data in input file\n"); 
                        exit(1);
                    }
                    if(posrays <= 0 || negrays <= 0) {
                        printf("consistency error: rays %d/%d (positive/negative)\n",
                            posrays, negrays); 
                        exit(1);
                    }
                    rewind(ifile); // prepare for re-reading
                    return; 
                }
                printf("input file is corrupt\n"); exit(1);
            }
            nextray[i] = w;
            if(w < -1700 || w > 1700) {
                printf("large value in input file: %d\n", w); exit(1);
            }
        }
        w = intdot(axioms[axiomno[iteration]], nextray);
        if(w < 0) { negrays++; }
        else if(w > 0) { posrays++; }
        else { zerorays++; }
    }
}


/* ===================================================================== */
static void save_nextray(void) {
    for(int i = 0; i < DIM; i++) fprintf(ofile, "%d%c", nextray[i], i == DIM - 1 ? '\n' : ' ');
    newrays++;
}


/* ===================================================================== */
/* second scan; output file is open */

#include "bitmap.h"
static int *pos_store, *neg_store; // storing positive and negative rays

#define STSIZE    (DIM+1)    // the ray coordinate plus distance from next axiom


static void second_input_scan(void) {
    int w; 
    size_t total_memory;
    int pidx, nidx, zidx; // index of positive, negative and zero rays
    int rayidx;           // index of the actual ray
    /* allocate memory for positive and negative rays */
    neg_store = (int *)malloc(sizeof(int) * STSIZE * (negrays + posrays));
    if(neg_store == NULL) {
       printf("Cannot allocate memory for rays %d/%d\n", posrays, negrays);
       exit(1);
    }
    pos_store = neg_store + (STSIZE * negrays);
    total_memory = init_bitmaps(iteration + 1, posrays + negrays + zerorays);
    if(total_memory == 0) {
        printf("Cannot allocate memory for bitmaps: pos=%d, neg=%d, zero=%d\n",
            posrays, negrays, zerorays);
        exit(1);
    }
    printf("Iteration: %d, newax=%d, pos/zero/neg=%d/%d/%d\n",
        iteration, axiomno[iteration], posrays, zerorays, negrays);
    /* header of the next iteration */
    fprintf(ofile, "%d\n", iteration + 1);
    for(int i = 0; i <= iteration; i++) fprintf(ofile, " %d", axiomno[i]);
    fprintf(ofile, "\n");
    /* go over the input file again */
    if(fscanf(ifile, "%d", &w) != 1 || w != iteration) {
        printf("input file consistency error #1\n"); exit(1);
    }
    for(int i = 0; i < iteration; i++)
        if(fscanf(ifile, "%d", &w) != 1 || w != axiomno[i]) {
            printf("input file consistency error #2\n"); exit(1);
        }
    pidx = 0; 
    nidx = 0; 
    zidx = 0;
    while(1) {  // read rays
        for(int i = 0; i < DIM; i++) {
            if(fscanf(ifile, "%d", &w) != 1) {
                if(i == 0 && feof(ifile)) { // sanity check
                    if(pidx != posrays || nidx != negrays || zidx != zerorays) {
                        printf("second scan inconsistency\n"); exit(1);
                    }
                    fclose(ifile);
                    return;
                }
                printf("input file consistency error #3\n"); exit(1);
            }
            nextray[i]=w;
        }
        w = intdot(axioms[axiomno[iteration]], nextray);
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
        if(w >= 0) { // it is also a ray in the next iteration
            save_nextray();
        }
        // check if the ray is adjacent to the first axioms
        for(int ax = 0; ax < iteration; ax++) {
            w = intdot(axioms[axiomno[ax]], nextray);
            if(w < 0){ printf("ray is on the negative side of axiom %d\n", axiomno[ax]); exit(1); }
            if(w == 0) // it is adjacent
                add_bitmap(ax, rayidx);
        }
    }
}


/* ---------------------------------------------------------------------------- */

/* fast iterative gcd algorithm */
static int gcd(int a, int b) {
    if(a == 1) return 1;
    if(a == 0) return b < 0 ? -b : b;
    if(a < 0) { a = -a; } if(b < 0) { b = -b; }
again:
    b %= a; if(b == 0) return a; if(b == 1) return 1;
    a %= b; if(a == 0) return b; if(a == 1) return 1;
    goto again;
}


static void simplify_nextray(void) {
    int v;

    v = nextray[0]; 
    if(v < 0) v = -v;
    for(int i = 1; v != 1 && i < DIM; i++) v = gcd(v, nextray[i]);
    if(v <= 0) { printf("nextray is all zero (%d)\n", v); exit(1); }
    if(v == 1) return;
    for(int i = 0; i < DIM; i++) nextray[i] /= v;
}


static void combine_rays(int posidx, int negidx) {
    const int *r1 = pos_store + STSIZE * posidx;
    const int *r2 = neg_store + STSIZE * negidx;
    int r1d = r1[DIM], r2d = -r2[DIM];
    if(r1d <= 0 || r2d <= 0) { printf("wrong sign\n"); exit(1); }
    for(int i = 0; i < DIM; i++) {
        nextray[i] = r2d * r1[i] + r1d * r2[i];
    }
#if TEST /* save some time not checking */
    int w = intdot(axioms[axiomno[iteration]], nextray);
    if(w != 0) { printf("wrong combine rays: %d %d\n", posidx, negidx); exit(1); }
#endif
    simplify_nextray();
    save_nextray();
}


/* ===================================================================== */
/* timer */
#include <sys/time.h>
static struct timeval prev_time;

static float elapsed_time(int set) {
    struct timeval current_time;
    
    if(set) { gettimeofday(&prev_time, NULL); return 0.0; }
    gettimeofday(&current_time, NULL);
    return (current_time.tv_sec - prev_time.tv_sec) + 
        (current_time.tv_usec - prev_time.tv_usec) / 1000.0 / 1000.0;
}


/* ===================================================================== */
/* find adjacent pos/neg ray pairs */
static size_t fastchkno = 0, combchkyes = 0, combchkno = 0;

static void printtime(float t) {
    int d, h, m, s;
    
    if(t < 59.994) { printf("%.2f", t); return; }
    s = (int)(t + 0.49);
    m = s / 60; 
    s = s % 60;
    if(m < 60) { printf("%d:%02d", m, s); return; }
    h = m / 60; 
    m = m % 60;
    if(h < 24) { printf("%d:%02d:%02d", h, m, s); return; }
    d = h / 24; 
    h = h % 24;
    printf("%dd %d:%02d:%02d", d, h, m, s);
}


static void statistics(void) {
    static float lastreport = 0.0; 
    float elapsed;
    
    if(verbose == 0) return;
    elapsed = elapsed_time(0);
    if(elapsed - lastreport < 30.0) return;
    lastreport = elapsed;
    printf("Elapsed "); 
    printtime(elapsed);
    printf(", rays %d, checked %4.2f%%, fast/no/yes=%zu/%zu/%zu\n",
       newrays,
       100.0 * ((float)(fastchkno + combchkyes + combchkno)) / ((float)posrays) / ((float)negrays),
       fastchkno, combchkno, combchkyes);
}


#include <sys/stat.h>
/* return 1 if file "skip" exists */

static int check_skipfile(void) {
    struct stat sb;
    return !skipfile ? 0 : lstat(skipfile, &sb) < 0 ? 0 : 1;
}


static void find_rays_negpos(void) { 
    int pairs_checked = 0;
    
    for(int nidx = 0; nidx < negrays; nidx++) {
        prepare_adjacency(nidx, DIM - 2);
        for(int pidx = 0; pidx < posrays; pidx++) {
            if(pairs_checked >= 5000000) {
                pairs_checked = 0;
                statistics();
                if(check_skipfile()) {
                    printf("Program run aborted by skipfile %s, exiting\n", skipfile); exit(1);
                }
            }
            if(fast_raycheck(negrays + pidx)) {
                fastchkno++;
            } else if(are_adjacent_rays(negrays + pidx)) {
                combchkyes++; 
                pairs_checked++;
                // rays nidx, pidx are adjacent, generate the corresponding ray
                combine_rays(pidx, nidx);
            } else {
                combchkno++; 
                pairs_checked++; 
            }
        }
    }
}


static void find_rays_posneg(void) { 
    int pairs_checked = 0;
    
    for(int pidx = 0; pidx < posrays; pidx++) {
        prepare_adjacency(negrays + pidx, DIM - 2);
        for(int nidx = 0; nidx < negrays; nidx++) {
            if(pairs_checked >= 5000000) {
                pairs_checked = 0;
                statistics();
                if(check_skipfile()) {
                    printf("Program run aborted by skipfile %s, exiting\n",skipfile); exit(1);
                }
            }
            if(fast_raycheck(nidx)){ 
                fastchkno++;
            } else if(are_adjacent_rays(nidx)){
                combchkyes++; 
                pairs_checked++;
                combine_rays(pidx, nidx);
            } else {
                combchkno++; 
                pairs_checked++;
            }
        }
    }
}


static void find_adjacent_rays(void) {   
    if(negrays < posrays) find_rays_negpos();
    else find_rays_posneg();
}


/* ===================================================================== */

int main(int argc, const char *argv[]) {
    float firsttime, secondtime, gentime;
    
    handle_params(argc, argv);
    if(check_skipfile()) {  printf("Skipfile %s exists, aborting\n", skipfile); return 1; }
    elapsed_time(1);
    first_input_scan(); // determine number of positive/negative rays
    firsttime = elapsed_time(0); 
    elapsed_time(1);
    second_input_scan();
    secondtime = elapsed_time(0); 
    elapsed_time(1);
    find_adjacent_rays();
    gentime = elapsed_time(0);
    fclose(ofile);
    printf("Rays=%d, fast/no/yes=%zu/%zu/%zu, time=",
        newrays, fastchkno, combchkno, combchkyes);
    printtime(firsttime + secondtime + gentime);
    printf(" (");
    printtime(firsttime + secondtime);
    printf(")\n");
    return 0;
}


/* EOF */
