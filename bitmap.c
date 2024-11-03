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

/* bitmap declarations and functions, Version 2.0 */

/* The bitmap algorithm works as follows. For an axiom a, B(a) is the set
of rays which satisfy a: B(a) = {u: u*a=0}. B(a) is stored as a bitmap,
1 bit for each extremal ray.

Given extremal rays u and v, store { address(B(a)): a \in  Z(u) \cap Z(v) }
u and v is adjacent iff the intersection of the bitmaps in this set
contains u and v only.

For efficiency the bitmaps are intersected in 64 bit chunks, stopping as
soon as the intersection is empty. Before taking the intersection (and
checking whether it is empty), bits corresponding to u and v are cleared
(if they fall into the corresponding 64 bit chunk).

Using Zolotykh's observation, the following variant has been implemented.
Fix the ray u, it is going to be compared to many other rays v.
The first step is to check if |Z(u)\cap Z(v)| >= d-2. If not, (u,v) are 
not adjacent (this is fast checking). If yes, create the bitmap
        C(w)=0 if w=u, or if |Z(u) \cap Z(w)| < d-2,
where w runs over all rays. By the observation, if Z(w) extends 
Z(u) \cap Z(v), then C(w)=1. Consequently u and v are adjacent iff the 
intersection of C(w) and the bitmaps B(a) for a \in Z(u) \cap Z(v) contains
v only (as u has been cleared in the bitmap C ).

The bitmap C is generated once for each extremal ray u, and then it
is used for subsequent tests where u is one element of the pair.
Also, checking |Z(u) \cap Z(v)| >= d-2 can be done by looking up the
corresponding bit from C(w).
*/

#include "bitmap.h"
#include <stdint.h>    /* types uint32_t, uint64_t */
#include <stdlib.h>    /* malloc, free */
#include <strings.h>   /* bzero */

/* bitmaps are unsigned words of 32 or 64 bits; have type BITMAP_t
*   declare a bitmap array containing "total" many bits:
*       BITMAP_t bmap[bmarray_size(total)];
*   W2~within a BITMAP_t word  bits are numbered starting from the least
*   significant position (i.e., 1)
*/

#ifdef DEBUG
#  include <stdio.h>
#  define ASSERT(expr)  ((void)((expr)||(assert_fails(#expr,__LINE__),1)))
  static void assert_fails(const char* text, int line) {
    printf("Line %d: assert \"%s\" failed\n", line, text);
  }
#else
#  define ASSERT(expr)    /* empty */
#endif

/* decide the bitmap size. We check here UINT64_MAX to be defined; 
   it can be changed to anything else */

#ifdef UINT64_MAX        /* using 64 bit bitmaps */
  #define BITMAP_64
  typedef uint64_t BITMAP_t;
  #define packsizelog   3    /* sizeof(BITMAP_t) == 1 << packsizelof */
#else
  #warning Using 32 bit bitmaps
  #define BITMAP_32
  typedef uint32_t BITMAP_t;
  #define packsizelog   2
#endif

#define packsize    (1<<packsizelog)    /* 4 or 8 bytes */
#define packshift   (packsizelog+3)     /* in bits, 5 or 6 */
#define bitsperword (1<<packshift)      /* 32 or 64 */
#define packmask    ((1<<packshift)-1)  /* 31 or 63 */
#define BITMAP1     ((BITMAP_t)1u)
#define BITMAP0     ((BITMAP_t)0u)

/* Basic macros handling a bitmap array with more than 32 (or 64)
*  bitcount:
*   clear / set bit "cnt" in a bitmap array:
*       clear_bit(b,cnt); set_bit(b,cnt);
*   extract bit "cnt" returning zero or one:
*       extract_bit(b,cnt)
*   check if bit "cnt" is set this way:
*       if(extract_bit(b.cnt))
*   find the index of a bitmap array containing bit number "cnt":
*       b[(cnt)>>packshift]
*   number of BITMAP_t words containing cnt many bits
*/

/* bool extract bit(BITMAP-t bm[], int cnt)
   void clear_bit  (BITMAP_t bm[], int cnt)
   void set_bit    (BITMAP_t bm[], int cnt) */
#define extract_bit(bm, cnt)  \
    (((bm)[(cnt)>>packshift]>>((cnt)&packmask))&1)

#define clear_bit(bm, cnt)    \
    (bm)[(cnt)>>packshift] &= ~(BITMAP1<<((cnt)&packmask))

#define set_bit(bm, cnt)      \
    (bm)[(cnt)>>packshift] |= (BITMAP1<<((cnt)&packmask))

/* int bmarray_size(int cnt)
*   number of BITMAP_t words to contain cnt bits
*/
#define bmarray_size(cnt)     \
    (((cnt)+packmask)>>packshift)

/* int get bitcount(BITMAP_t v)
*     count the number of bits in a BITMAP word
*/
#if __x86_64__ || __amd64__    /* assembler code on these architectures */

#warning Using x86 assembler code
static inline int get_bitcount(BITMAP_t v) {
    register BITMAP_t res;
    asm ("popcnt %[w], %[t]"
         :[t] "=rm" (res)
         :[w] "rm"  (v));
    return (int)res;
}

#else /* not x86 architecture */
static int get_bitcount(BITMAP_t v) {
//    these exceptional cases do not seem to help
//    if(v==BITMAP0){ return 0; }
//    if(v==~BITMAP0){ return 1<<packshift; }
#   ifdef BITMAP_32     /* 32 bit bitmap */
    v = (v & 0x55555555u) + ((v >> 1) & 0x55555555u);
    v = (v & 0x33333333u) + ((v >> 2) & 0x33333333u);
    v = (v & 0x0F0F0F0Fu) + ((v >> 4) & 0x0F0F0F0Fu);
    v = (v & 0x00FF00FFu) + ((v >> 8) & 0x00FF00FFu);
    return (int)((v & 0x0000FFFF) + (v >> 16));
#   else                /* 64 bit bitmap */
    v = (v & 0x5555555555555555ul) + ((v >> 1) & 0x5555555555555555ul);
    v = (v & 0x3333333333333333ul) + ((v >> 2) & 0x3333333333333333ul);
    v = (v & 0x0F0F0F0F0F0F0F0Ful) + ((v >> 4) & 0x0F0F0F0F0F0F0F0Ful);
    v = (v & 0x00FF00FF00FF00FFul) + ((v >> 8) & 0x00FF00FF00FF00FFul);
    v = (v & 0x0000FFFF0000FFFFul) + ((v >> 16) & 0x0000FFFF0000FFFFul);
    return (int)((v & 0xFFFF) + (v >> 32));
#   endif        /* 64 bit bitmap */
}

#endif /* x86 assembler code */


/* ======================================================================= */
/* Adjacency bitmaps of axioms are in
*      BITMAP_t* AxAdj[AXIOMS];
*  each AxAdj(i) is a BITMAP_t array of length bmarray_size(RAYS).
*  Adjacency bitmaps of rays are in
*      BITMAP_t* RayAdj[RAYS];
*  each RayAdj(j) is a BITMAP_t array of length bmarray_size(AXIOMS).
*/

static int AXIOMS,RAYS;
static int AxAdjSize;          // bmarray_size(AXIOMS)
static int RayAdjSize;         // bmarray_size(RAYS)
static BITMAP_t *AxAdjStore;   // AXIOMS bitmap arrays, each for RAYS bits
static BITMAP_t *RayAdjStore;  // RAYS bitmap arrays, each for AXIOMS bits
static BITMAP_t *BigRays;      // bitmap storing rays with large intersection
static BITMAP_t **joint;       // AXIOM bitmaps for intersection

#define AxAdj(ax)    \
    (AxAdjStore + (ax * RayAdjSize))
#define RayAdj(ray)  \
    (RayAdjStore + (ray * AxAdjSize))

/* allocate memory for the bitmaps
*    return 0 on memory error, otherwise the allocated memory size
*/

size_t init_bitmaps(int axiomno, int rayno) {
    AXIOMS = axiomno; 
    RAYS=rayno;
    AxAdjSize = bmarray_size(AXIOMS);
    RayAdjSize = bmarray_size(RAYS);
    if(!(AxAdjStore = (BITMAP_t*)malloc(sizeof(BITMAP_t) * AXIOMS * RayAdjSize)))
        return 0;
    if(!(RayAdjStore = (BITMAP_t*)malloc(sizeof(BITMAP_t) * RAYS * AxAdjSize)))
        return 0;
    if(!(BigRays = (BITMAP_t*)malloc(sizeof(BITMAP_t) * RayAdjSize)))
        return 0;
    if(!(joint = (BITMAP_t**)malloc(sizeof(BITMAP_t*) * AXIOMS)))
        return 0;
    bzero(AxAdjStore, sizeof(BITMAP_t) * AXIOMS * RayAdjSize);
    bzero(RayAdjStore, sizeof(BITMAP_t) * RAYS * AxAdjSize);
    return sizeof(BITMAP_t) * (AXIOMS + 1) * RayAdjSize + sizeof(BITMAP_t) * RAYS * AxAdjSize
        + sizeof(BITMAP_t*) * AXIOMS;
}


void close_bitmaps(void) {
    free(AxAdjStore); 
    free(RayAdjStore); 
    free(joint);
}


void add_bitmap(int axiom, int ray) {
    ASSERT(0 <= axiom && axiom < AXIOMS);
    ASSERT(0 <= ray && ray < RAYS);
    set_bit(AxAdj(axiom), ray);
    set_bit(RayAdj(ray), axiom);
}


/* Given ray r1, create the ray bitmap BigRays as
*     B(r2)=1 if r1!=r2 and intersect_size(r1,r2)>=dim-2
*/
static BITMAP_t *RayAdjR1;

void prepare_adjacency(int r1, int DIMminus2) {
    register BITMAP_t *L2; 
    BITMAP_t *g; 
    BITMAP_t onebit;
    
    bzero(BigRays, sizeof(BITMAP_t) * RayAdjSize);
    g = BigRays; 
    L2 = RayAdj(0); 
    RayAdjR1 = RayAdj(r1); 
    onebit = BITMAP1;
    for(int r2 = 0; r2 < RAYS; r2++) {
        if(__builtin_expect(r2 != r1, 1)) {
            int total = 0; 
            register BITMAP_t *L1 = RayAdjR1;
            for(int i = 0; i < AxAdjSize; i++, L1++, L2++)
                total += get_bitcount((*L1) & (*L2));
            if(total >= DIMminus2) *g |= onebit;
        } else {
           L2 += AxAdjSize;
        }
        if((onebit <<= 1) == 0) { onebit = BITMAP1; g++; }
    }
}


int fast_raycheck(int r2) { 
    return extract_bit(BigRays, r2) == 0;
}


/* Rays r1 and r2 are adjacent, it there are at least dim-2 
*  axioms containing both of them, moreover intersection the
*  bitmaps of these axioms, only r1 and r2 should remain.
*  We store the bitmap pointers of the joint axioms in the
*  array joint[].
*/

int are_adjacent_rays(int r2) {
    int nax = 0;             // number of axioms
    int idx2;                // bitmap word counter
    register BITMAP_t *L1, *L2;
    
    L1 = RayAdjR1; 
    L2 = RayAdj(r2);
    for(int i = 0; i < AXIOMS; i += bitsperword, L1++, L2++) {
        int idx = i;
        BITMAP_t w = (*L1) & (*L2);
        while(w) {
            while((w & 0xF) == 0) { idx += 4; w >>= 4; }
            if(w & 1) { joint[nax] = AxAdj(idx); nax++; }
            idx++; 
            w >>= 1;
        }
    }
    // prepare the intersection of BigRays and the first axiom
    idx2 = r2 >> packshift;
    // go over all bitmap words, abort as soon as it is empty
    L1 = BigRays; 
    L2 = joint[0];
    for(int i = 0; i < RayAdjSize; i++, L1++, L2++) {
        register BITMAP_t w = (*L1) & (*L2);
        if(__builtin_expect(i == idx2, 0)) w &= ~(BITMAP1 << (r2 & packmask));
        for(int j = 1; w != 0 && j < nax; j++) { w &= (joint[j])[i]; }
        if(w != 0) return 0;
    }
    // the intersection is empty; they are adjacent
    return 1;
}


/* EOF */

