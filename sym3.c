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

#define VPERMNO     6

T_FACTOR VPERMS[VPERMNO][VARS]={ /* permutations of variables */
{0,1,2,3}, /* 0 */
{1,0,2,3}, /* 1 */
{3,0,2,1}, /* 2 */
{0,3,2,1}, /* 3 */
{3,1,2,0}, /* 4 */
{1,3,2,0}  /* 5 */
};

T_FACTOR dualmatrix[VARS][VARS]={ /* computing the dual ray */
/* 0 */ {1,0,0,0},
/* 1 */ {0,1,0,0},
/* 2 */ {1,1,-1,1},
/* 3 */ {0,0,0,1}
};

