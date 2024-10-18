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

#define VPERMNO     24

T_FACTOR VPERMS[VPERMNO][VARS]={ /* permutations of variables */
{0,1,2,3,4,5,6,7,8,9,10}, /* 0 */
{1,0,2,3,5,4,6,7,9,8,10}, /* 1 */
{7,0,8,1,9,2,6,3,10,4,5}, /* 2 */
{0,7,8,1,2,9,6,3,4,10,5}, /* 3 */
{7,1,9,0,8,2,6,3,10,5,4}, /* 4 */
{3,7,10,0,4,8,6,1,5,9,2}, /* 5 */
{1,7,9,0,2,8,6,3,5,10,4}, /* 6 */
{7,3,10,0,8,4,6,1,9,5,2}, /* 7 */
{3,0,4,7,10,8,6,1,5,2,9}, /* 8 */
{3,7,10,1,5,9,6,0,4,8,2}, /* 9 */
{1,3,5,7,9,10,6,0,2,4,8}, /* 10 */
{0,3,4,7,8,10,6,1,2,5,9}, /* 11 */
{7,3,10,1,9,5,6,0,8,4,2}, /* 12 */
{3,1,5,7,10,9,6,0,4,2,8}, /* 13 */
{1,7,9,3,5,10,6,0,2,8,4}, /* 14 */
{1,3,5,0,2,4,6,7,9,10,8}, /* 15 */
{7,1,9,3,10,5,6,0,8,2,4}, /* 16 */
{3,1,5,0,4,2,6,7,10,9,8}, /* 17 */
{0,7,8,3,4,10,6,1,2,9,5}, /* 18 */
{0,3,4,1,2,5,6,7,8,10,9}, /* 19 */
{0,1,2,7,8,9,6,3,4,5,10}, /* 20 */
{7,0,8,3,10,4,6,1,9,2,5}, /* 21 */
{3,0,4,1,5,2,6,7,10,8,9}, /* 22 */
{1,0,2,7,9,8,6,3,5,4,10}  /* 23 */
};

T_FACTOR dualmatrix[VARS][VARS]={ /* computing the dual ray */
/* 0 */ {1,0,0,0,0,0,0,0,0,0,0},
/* 1 */ {0,1,0,0,0,0,0,0,0,0,0},
/* 2 */ {1,1,0,0,0,0,-1,0,0,0,1},
/* 3 */ {0,0,0,1,0,0,0,0,0,0,0},
/* 4 */ {1,0,0,1,0,0,-1,0,0,1,0},
/* 5 */ {0,1,0,1,0,0,-1,0,1,0,0},
/* 6 */ {1,1,0,1,0,0,-1,1,0,0,0},
/* 7 */ {0,0,0,0,0,0,0,1,0,0,0},
/* 8 */ {1,0,0,0,0,1,-1,1,0,0,0},
/* 9 */ {0,1,0,0,1,0,-1,1,0,0,0},
/* 10 */ {0,0,1,1,0,0,-1,1,0,0,0}
};

