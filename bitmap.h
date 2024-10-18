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

/* bitmap.h -- bitmap routines */

/* The bitmap part of a single iteration of the DD method.
*  There are AXIOMS many axioms and RAYS many rays.
*  Bitmaps are initialized by calling 
*     init_bitmaps(axiomno,rayno,dimension)
*  which specifies the number of axioms, rays, and the dimension, allocates
*  space; returns -1 if cannot allocate the bitmaps; otherwise the number of
*  allocated bytes.
*     add_bitmap(axiom,ray)
*  establishes that "ray" satisfies "axiom". Should be called for all such
*  pairs to fill the bitmaps.
*     prepare_adjacency(r1,dim-2)
*  prepares a bulk check of ray pair adjacency.
*     fast_raycheck(r2)
*  checks if ray r2 intersect r1 in size less than dim-2; returns if yes
*     are_adjacent_rays(r2)
*  returns 1 if rays r1 and r2 are adjacent, and 0 otherwise.
*/

#ifndef BITMAP_H
#define BITMAP_H
#include <stddef.h>


size_t init_bitmaps(int axiomno, int rayno);
/* return value: 0: allocate error, otherwise total memory */
void close_bitmaps(void);
/* release all memory allocated for the bitmaps */

void add_bitmap(int axiom,int ray);
/* axiom with index "axiom" vanished at ray with index "ray".
   if DEBUG is defined, checks the range of arguments */

void prepare_adjacency(int r1,int DIMminus2);
/* prepare bulk check of rays starting with r1 */
int fast_raycheck(int r2);
/* number of axioms adjacent to both r1 and r2 is too small */

int are_adjacent_rays(int r2);
/* returns 1 it r1 and r2 are adjacent, 0 otherwise. */

#endif

/* EOF */

