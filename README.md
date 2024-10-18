
# Enumerating extremal submodular functions

This repository contains
supplementary material for the paper:

E. P. Csirmaz and L. Csirmaz:
Attempting the impossible: enumerating extremal submodular functions for n=6.

## Build

To build the necessary binaries, ensure `gcc` is installed, and run `build.sh`.

## Examples

To compute all extremal rays for n=4 and output the results to the
file `/tmp/ex4-24.txt` with recursive order, use
`perl base.pl 4 --purge --rec /tmp/ex4`.

To generate all adjacent extremal rays of the Vamos-ray defined in `vamos.txt`
and output the results to the
file `/tmp/v-all.txt`, use
`perl runjob.pl 4 --ray vamos.txt /tmp/v`.

## Files

`ax<N>.c` -  p-standardized submodular inequalities in C structure for base set
sizes N=3,4,5,6. Human-readable versions of the inequalities are included.

`sym<N>.c` - Conversion matrices for computing permutations and the dual, as
C structures.

`bitmap.c`, `bitmap.h` - Implementation of bitmap-manipulating routines as required by
the DD method. It uses assembler instructions to count the bits in a 64-bit word.

`iterate.c` - Implementation of a single iteration of the DD method as a pipe.
It reads the description of the problem, the next inequality, and then
produces the result of the iteration. Should be compiled for each
base size separately. Used by the controller program `base.pl`.
To compile, use
`gcc -O3 -o bin/iter<N> -D RANK=<N> iterate.c bitmap.c`

`neig.c` - A version of `iterate.c` for computing the neighbors of a ray. Used by
the controller program `runjob.pl`. Should be compiled for each base
size separately. To compile, use
`gcc -O3 -o bin/neig<N> -D RANK=<N> neig.c bitmap.c`

`newrays.c` - Generates the extremal rays from the neighbor enumeration. It can
produce a canonical representative from the corresponding orbit
rather than the ray itself. Used by the controller program `runjob.pl`.
Sould be compiled for each base size separately. To compile, use
`gcc -O3 -o bin/newr<N> -D RANK=<N> newrays.c`

`extract.c` - Generate extremal rays from an intermediate cone of the DD method
as produced by the controller program `base.pl`. It can produce a
canonical representative from the corresponding orbit rather than
the ray itself. Should be compiled for each base size separately.
To compile, use
`gcc -O3 -o bin/extr<N> -D RANK=<N> extract.c bitmap.c`

`base.pl` - Controller program to execute all iterations of the DD method for
the submodular cone. Arguments are the base size and the ordering
method. Generates the initial approximation, then calls `bin/iter<N>`
until the algorithm terminates. When stopped, the last result can be
used as input to `bin/extr<N>`. 

`runjob.pl` - Controller program to compute rays adjacent to a list of extremal
rays. Arguments are the base size, ordering method, and the list of
extremal rays. For each ray generates the initial approximation,
and calls `bin/neig<N>` until the DD method completes. Finally calls
`bin/newr<N>` to generate the neighbors or their canonical representatives.

`vamos.txt` - One of the six Vamos extremal rays of the 4-element submodular
cone.

`build.sh` - Script to call `gcc` to build all binaries.
