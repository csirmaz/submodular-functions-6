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
#define LABEL "SET_N=3 TIGHT=True"
#define AXIOMS 6
#define VARS 4
T_FACTOR axioms[AXIOMS][VARS] = {
{1,1,-1,0} /* #0 f({0}) + f({1}) >= f({0, 1, 2})*/,
{0,0,1,-1} /* #1 f({0, 1, 2}) >= f({2})*/,
{1,0,-1,1} /* #2 f({0}) + f({2}) >= f({0, 1, 2})*/,
{0,-1,1,0} /* #3 f({0, 1, 2}) >= f({1})*/,
{0,1,-1,1} /* #4 f({1}) + f({2}) >= f({0, 1, 2})*/,
{-1,0,1,0} /* #5 f({0, 1, 2}) >= f({0})*/
};
#ifdef READABLE
char* human_readable_axioms[AXIOMS] = {
"#0 f({0}) + f({1}) >= f({0, 1, 2})",
"#1 f({0, 1, 2}) >= f({2})",
"#2 f({0}) + f({2}) >= f({0, 1, 2})",
"#3 f({0, 1, 2}) >= f({1})",
"#4 f({1}) + f({2}) >= f({0, 1, 2})",
"#5 f({0, 1, 2}) >= f({0})"
};
#endif
