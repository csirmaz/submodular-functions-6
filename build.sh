#!/bin/bash

# This file is supplemental material to the paper
# "Attempting the impossible: enumerating extremal submodular functions for n=6"
# by E P Csirmaz and L Csirmaz.
#
# Copyright 2024 E P Csirmaz and L Csirmaz
#
# This program is free software: you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the Free
# Software Foundation, either version 3 of the License, or (at your option)
# any later version.
#
# This program is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
# FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
# more details.
#
# You should have received a copy of the GNU General Public License along with
# this program.  If not, see <https://www.gnu.org/licenses/>.

mkdir bin
for n in 4 5 6; do
gcc -O3 -o bin/iter$n -D RANK=$n iterate.c bitmap.c
gcc -O3 -o bin/neig$n -D RANK=$n neig.c bitmap.c
gcc -O3 -o bin/newr$n -D RANK=$n newrays.c
gcc -O3 -o bin/extr$n -D RANK=$n extract.c bitmap.c
done

