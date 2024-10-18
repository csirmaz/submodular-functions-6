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

/* Iterate the DD method to find the neighbors of an extremal ray.
*  Only a single axiom is fixed, see Section IV-B of the paper.
*  
*  Usage: neig <n> <axno> <input> <output> <stopfile>
*    n      iteration count, inequalities handled so far
*    axno   inequality number to be added next
*    input  input file, see description below
*    output output file, can be used for the next iteration
*    stopfile abort running if the file exists
*
*  input file format:
*    only defimal numbers separated by white space
*    1 <fixed ineq>   this is 1, followed by the number of an inequality
*    <iteration count>
*    <ineq number>    inequality numbers handled so far
*    <ray1> <ray2> ... <rayn>
*    each ray is a list of DIM many (possibly signed) integers
* 
*  We assume that no coefficient has absolute value above 1700,
*  so integer arithmetic can be used. This condition is checked.
*
*  This is a modified version of "iterate.c".
*
*/

#define T_FACTOR	static int
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

#define DIM	VARS

#define stringify(x)	#x
#define mkstringof(x)	stringify(x)

/* ============================================================================ */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static int verbose=0; 		// print additional info
static int iteration=0;		// iteration we are in
#define fixed	1		// axioms fixed
static int axiomno[AXIOMS];	// axioms handled and fixed
static FILE *ifile, *ofile;	// input and output file
static const char* skipfile;	// the skipfile

static void usage(void){
    printf(
"usage: neig" mkstringof(RANK) " [-v|-q] <n> <newax> <input> <output> [<skipfile>]\n"
"iterating partial extremal rays algorithms for rank=" mkstringof(RANK) "\n"
" -v        verbose (printing additional information)\n"
" -q        quiet, no printing\n"
"  n        iteration count; number of axioms handled before\n"
"  newax    next axiom number to be added next\n"
"  input    input file containing the last generation of extremal rays\n"
"  output   output file for the next generation of extremal rays\n"
"  skipfile abort the program if exists; default is \"skip\"\n"
"The output file should not exist.\n");
    exit(1);
}

static void handle_params(int argc,const char*argv[])
{const char *w, *input,*output; int n;
    if(argc>1 && argv[1][0]=='-' && argv[1][1]=='v'){
        verbose=1; argc--; argv++;
    }
    if(argc>1 && argv[1][0]=='-' && argv[1][1]=='q'){
        verbose=2; argc--; argv++;
    }
    if(argc!=5 && argc!=6) usage();
    skipfile=argc==6 ? argv[5] : "skip";
    w=argv[1]; iteration=0;
    while(*w && '0'<=*w && *w<='9'){ iteration =iteration*10+(*w-'0'); w++; }
    if(*w || iteration<2 || iteration+fixed>=AXIOMS-1){ printf("wrong arg <n>\n"); exit(1); }
    n=0; w=argv[2];
    while(*w && '0'<=*w && *w<='9'){ n = n*10+(*w-'0'); w++; }
    if(*w || n>=AXIOMS){ printf("wrong arg <newax>\n"); exit(1); }
    input=argv[3]; output=argv[4];
    // 0 .. iteration-1 : previous iteration;
    // iteration+1 .. iteration+fixed : fixed axioms
    axiomno[iteration]=n;
    if(!input || !*input || (ifile=fopen(input,"r"))==NULL){
         printf("input file \"%s\" not found\n",input);
         exit(1);
    }
    if(!output || !*output || (ofile=fopen(output,"r"))!=NULL){
         printf("output file \"%s\" exists, cannot continue\n",output);
         exit(1);
    }
    ofile=fopen(output,"w");
    if(ofile==NULL){
         printf("Cannot create oputput file \"%s\"\n",output);
         exit(1);
    }
}

/* =================================================================== */
static int posrays=0,negrays=0,zerorays=0;	// number of rays
static int nextray[DIM];			// reading the next ray
static int newrays=0;				// rays in the output
#define MINS	(DIM-2-fixed)			// minimal intersection size

static int intdot(const int*a,const int*b)
{int sum=0;
    for(int i=0;i<DIM;i++) sum+=a[i]*b[i];
    return sum;
}

/* read the input file first; compute the number of positive, negative rays
*  relative to the next axiom */
static void first_input_scan(void)
{int w,i;
   // number of fixed axioms
   w=0;if(fscanf(ifile,"%d",&w)!=1 || w!=fixed){
      printf("error in number of fixed axioms %d\n",w); exit(1);
   }
   // fixed axioms
   for(i=0;i<fixed;i++){
       w=0;if(fscanf(ifile,"%d",&w)!=1 || w<0 || w>=AXIOMS){
           printf("error in the %d-th fixed axiom: %d\n",i,w); exit(1);
       }
       axiomno[iteration+1+i]=w;
   }
   // handled axioms so far
   w=0;if(fscanf(ifile,"%d",&w)!=1 || w!=iteration){
       printf("iteration number %d error in input file\n",w); exit(1);
   }
   for(i=0;i<iteration;i++){
       w=0;if(fscanf(ifile,"%d",&w)!=1 || w<0 || w>=AXIOMS){
         printf("axiom list in input file is corrupted\n"); exit(1);
       }
       axiomno[i]=w;
   }
   // axiomno[0 .. iteration+fixed]  must all be different
   for(i=1;i<=iteration+fixed;i++){
       w=axiomno[i];
       for(int j=0;j<i;j++) if(axiomno[j]==w){
         printf("axiomno %d is repeated (%d,%d)\n",w,i,j); exit(1);
       }
   }
   while(1){ // read rays ...
      for(i=0;i<DIM;i++){
        if(fscanf(ifile,"%d",&w)!=1){
            if(i==0 && feof(ifile)){ // last item read
                if(posrays+negrays+zerorays<DIM-fixed){
                    printf("Not enough data in input file\n"); exit(1);
                }
                if(posrays<=0 || negrays<=0){
                    printf("consistency error: pos/zero/neg rays %d/%d/%d\n",
                         posrays,zerorays,negrays); exit(1);
                }
                rewind(ifile); // prepare for re-reading
                return; 
            }
            printf("input file is corrupt\n"); exit(1);
        }
        nextray[i]=w;
        if(w<-1700 || w>1700){
             printf("Large value in input file: %d\n",w); exit(1);
        }
      }
      w=intdot(axioms[axiomno[iteration]],nextray);
      if(w<0){negrays++;}
      else if(w>0){ posrays++; }
      else {zerorays++; }
   }
}

/* ===================================================================== */
/* second scan; output file is open */

#include "bitmap.h"
static int *pos_store, *neg_store; // storing positive and negative rays

#define STSIZE	(DIM+1)	// the ray coordinate plus distance from next axiom

static void second_input_scan(void)
{int w,i; size_t total_memory;
 int pidx,nidx,zidx; // index of positive, negative and zero rays
 int rayidx;         // index of the actual ray
    /* allocate memory for positive and negative rays */
    neg_store=(int *)malloc(sizeof(int)*STSIZE*(negrays+posrays));
    if(neg_store==NULL){
       printf("Cannot allocate memory for rays %d/%d\n",posrays,negrays);
       exit(1);
    }
    pos_store=neg_store+(STSIZE*negrays);
    total_memory=init_bitmaps(iteration+1,posrays+negrays+zerorays);
    if(total_memory==0){
       printf("Cannot allocate memory for bitmaps: pos=%d, neg=%d, zero=%d\n",
          posrays,negrays,zerorays);
       exit(1);
    }
    if(verbose!=2){
       printf("Iteration: %d, newax=%d,  pos/zero/neg=%d/%d/%d\n",
           iteration,axiomno[iteration],posrays,zerorays,negrays); 
    }
    /* header of the next iteration */
    fprintf(ofile,"%d",fixed);
    for(i=1;i<=fixed;i++) fprintf(ofile," %d",axiomno[iteration+i]);
    fprintf(ofile,"\n%d",iteration+1);
    for(i=0;i<=iteration;i++) fprintf(ofile," %d",axiomno[i]);
    /* go over the input file again */
    if(fscanf(ifile,"%d",&w)!=1 || w!=fixed){
       printf("input file consistency error #1\n"); exit(1);
    }
    for(i=0;i<fixed;i++) if(fscanf(ifile,"%d",&w)!=1 || w!=axiomno[iteration+1+i]){
       printf("input file consistency error #2\n"); exit(1);
    }
    if(fscanf(ifile,"%d",&w)!=1 || w!=iteration){
       printf("input file consistency error #3\n"); exit(1);
    }
    for(i=0;i<iteration;i++) if(fscanf(ifile,"%d",&w)!=1 || w!=axiomno[i]){
       printf("input file consistency error #4\n"); exit(1);
    }
    pidx=0; nidx=0; zidx=0;
    while(1){ // read rays
        for(i=0;i<DIM;i++){
          if(fscanf(ifile,"%d",&w)!=1){
              if(i==0 && feof(ifile)){ // sanity check
                   if(pidx!=posrays || nidx!=negrays || zidx!=zerorays){
                       printf("second scan inconsistency, new rays=%d\n",newrays); exit(1);
                   }
                   fprintf(ofile,"\n");
                   fclose(ifile);
                   return;
              }
              printf("input file consistency error #5\n"); exit(1);
          }
          nextray[i]=w;
        }
        w=intdot(axioms[axiomno[iteration]],nextray);
        if(w<0){ // negative
            memcpy(neg_store+STSIZE*nidx,nextray,sizeof(int)*DIM);
            neg_store[STSIZE*nidx+DIM]=w;
            rayidx=nidx;
            nidx++;
        } else if(w>0){ // positive
            memcpy(pos_store+STSIZE*pidx,nextray,sizeof(int)*DIM);
            pos_store[STSIZE*pidx+DIM]=w;
            rayidx=pidx+negrays;
            pidx++;
        } else { // zero
            rayidx=zidx+negrays+posrays;
            zidx++;
        }
        if(w>=0){ // it is also a ray in the next iteration
            newrays++;
            fprintf(ofile,"\n%d",nextray[0]);
            for(i=1;i<DIM;i++) fprintf(ofile," %d",nextray[i]);
        }
        // check if the ray is adjacent to the first axioms
        for(int ax=0;ax<iteration;ax++){
            w=intdot(axioms[axiomno[ax]],nextray);
            if(w<0){printf("next ray is on negative side of axiom %d\n",axiomno[ax]); exit(1); }
            if(w==0) // it is adjacent
                add_bitmap(ax,rayidx);
        }
        // check if the ray is adjacent to the fixed axioms
        for(int ax=1;ax<=fixed;ax++){
            w=intdot(axioms[axiomno[iteration+ax]],nextray);
            if(w!=0){printf("next ray is not on the fixed axiom %d\n",axiomno[iteration+ax]); exit(1);}
        }
    }
}

/* ---------------------------------------------------------------------------- */

/* fast iterative gcd algorithm */
static int gcd(int a,int b){
    if(a==1) return 1;
    if(a==0) return b<0?-b:b;
    if(a<0){a=-a;} if(b<0){b=-b;}
again:
    b %= a; if(b==0) return a; if(b==1) return 1;
    a %= b; if(a==0) return b; if(a==1) return 1;
    goto again;
}

static void simplify_nextray(void)
{int v;
    v=nextray[0]; if(v<0) v=-v;
    for(int i=1;v!=1 && i<DIM;i++) v=gcd(v,nextray[i]);
    if(v<=0){printf("nextray is all zero (%d)\n",v); exit(1);}
    if(v==1) return;
    for(int i=0;i<DIM;i++) nextray[i] /= v; 
}

static void combine_rays(int posidx,int negidx)
{const int *r1=pos_store+STSIZE*posidx;
 const int *r2=neg_store+STSIZE*negidx;
    int r1d=r1[DIM],r2d=-r2[DIM];
    if(r1d<=0 || r2d<=0){ printf("wrong sign\n"); exit(1); }
    for(int i=0;i<DIM;i++){
        nextray[i]=r2d*r1[i]+r1d*r2[i];
    }
#if TEST /* save time by not checking */
    int w=intdot(axioms[axiomno[iteration]],nextray);
    if(w!=0){printf("wrong combine rays: %d %d\n",posidx,negidx); exit(1);}
#endif
    simplify_nextray();
    newrays++;
    for(int i=0;i<DIM;i++) fprintf(ofile,"%d%c",nextray[i],i==DIM-1?'\n':' ');
}

/* ===================================================================== */
/* timer */
#include <sys/time.h>
static struct timeval prev_time;
static float elapsed_time(int set)
{struct timeval current_time;
   if(set){ gettimeofday(&prev_time,NULL); return 0.0; }
   gettimeofday(&current_time,NULL);
   return (current_time.tv_sec-prev_time.tv_sec)+
          (current_time.tv_usec-prev_time.tv_usec)/1000.0/1000.0;
}

/* ===================================================================== */
/* find adjacent pos/neg ray pairs */
static size_t fastchkno=0, combchkyes=0, combchkno=0;

static void printtime(float t)
{int d,h,m,s;
    if(t<59.994){ printf("%.2f",t); return; }
    s=(int)(t+0.49);
    m=s/60; s = s %60;
    if(m<60){ printf("%d:%02d",m,s); return;}
    h=m/60; m= m % 60;
    if(h<24){ printf("%d:%02d:%02d",h,m,s); return; }
    d=h/24; h = h % 24;
    printf("%dd %d:%02d:%02d",d,h,m,s);
}

static void statistics(void){
static float next_report=29.5;
    if(verbose!=1) return;
    float elapsed=elapsed_time(0);
    if(elapsed<next_report) return;
    next_report=elapsed+29.5;
    printf("Elapsed: "); printtime(elapsed);
    printf(", rays=%d, checked %4.2f%%, fast/no/yes=%zu/%zu/%zu\n",
       newrays,
       100.0*((float)(fastchkno+combchkyes+combchkno))/((float)posrays)/((float)negrays),
       fastchkno,combchkno,combchkyes);
}

#include <sys/stat.h>
/* return 1 if file "skip" exists */
static int check_skipfile(void){
struct stat sb;
   return !skipfile ? 0 : lstat(skipfile,&sb)< 0 ? 0 : 1;
}

static void find_rays_negpos(void)
{int pairs_checked=0;
    for(int nidx=0;nidx<negrays;nidx++){
        prepare_adjacency(nidx,MINS);
        for(int pidx=0;pidx<posrays;pidx++){
            if(pairs_checked>=5000000){
                pairs_checked=0;
                statistics();
                if(check_skipfile()){
                    printf("Program run aborted by skipfile %s, exiting\n",skipfile);
                    exit(1);
                }
            }
            if(fast_raycheck(negrays+pidx)){
                fastchkno++;
            } else if(are_adjacent_rays(negrays+pidx)){
                combchkyes++; pairs_checked++;
                // rays nidx, pidx are adjacent, generate the corresponding ray
                combine_rays(pidx,nidx);
            } else {
                combchkno++; pairs_checked++; 
            }
        }
    }
}

static void find_rays_posneg(void)
{int pairs_checked=0;
    for(int pidx=0;pidx<posrays;pidx++){
        prepare_adjacency(negrays+pidx,MINS);
        for(int nidx=0;nidx<negrays;nidx++){
            if(pairs_checked>=5000000){
                pairs_checked=0;
                statistics();
                if(check_skipfile()){
                    printf("Program run aborted by skipfile %s, exiting\n",skipfile);
                    exit(1);
                }
            }
            if(fast_raycheck(nidx)){
                fastchkno++;
            } else if(are_adjacent_rays(nidx)){
                combchkyes++; pairs_checked++;
                // rays nidx, pidx are adjacent, generate the corresponding ray
                combine_rays(pidx,nidx);
            } else {
                combchkno++; pairs_checked++; 
            }
        }
    }
}

static void find_adjacent_rays(void)
{   if(negrays<posrays) find_rays_negpos();
    else find_rays_posneg();
}

/* ===================================================================== */

int main(int argc, const char*argv[])
{float firsttime,secondtime,gentime;
    handle_params(argc,argv);
    if(check_skipfile()){ printf("Skipfile %s exists, aborting\n",skipfile); return 1; }
    elapsed_time(1);
    first_input_scan(); // determine number of positive/negative rays
    firsttime=elapsed_time(0); elapsed_time(1);
    second_input_scan();
    secondtime=elapsed_time(0); elapsed_time(1);
    find_adjacent_rays();
    gentime=elapsed_time(0);
    fclose(ofile);
    if(verbose!=2){
      printf("Rays=%d, fast/no/yes=%zu/%zu/%zu, time=",
         newrays,fastchkno,combchkno,combchkyes);
      printtime(firsttime+secondtime+gentime);
      printf(" ("); printtime(firsttime+secondtime); printf(")\n");
    }
    return 0;
}


/* EOF */
