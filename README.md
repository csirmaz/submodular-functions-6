
# Enumerating extremal submodular functions

This repository contains
supplementary material for the paper:

E. P. Csirmaz and L. Csirmaz:
Attempting the impossible: enumerating extremal submodular functions for n=6.

The preprint is available on arXiv at https://arxiv.org/abs/2410.15502 .

## Results for N=6

A partial list of 360 billion extremal submodular functions for N=6 is available at
https://zenodo.org/records/13954788 .

## Quickstart

```bash
# Clone the repository
$ git clone https://github.com/csirmaz/submodular-functions-6.git
$ cd submodular-functions-6

# Compile the binaries
$ sudo apt-get install gcc
$ ./build.sh

# Example 1: Compute all extremal rays for N=4 using recursive order
$ perl base.pl 4 --purge --rec /tmp/ex4
# The results will be written to /tmp/ex4-24.txt

# Example 2: Generate all extremal rays adjacent to the Vamos-ray defined
# in vamos.txt
$ perl runjob.pl 4 --ray vamos.txt /tmp/v
# The results will be written to /tmp/v-all.txt

# Example 3: Continuing from the known N=6 rays in Zenodo, keep searching
# for neighboring rays. Download the latest rays:
$ wget https://zenodo.org/records/13954788/files/6rays-26.txt?download=1 -O 6rays-26.txt
$ perl runjob.pl 6 6rays-1.txt /tmp/out
# The results will be written to /tmp/out-all.txt in a
# "<weight>: <coordinates>" format.
```

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
base size N separately. Used by the controller program `base.pl`.
To compile for a given N, use
`gcc -O3 -o bin/iter<N> -D RANK=<N> iterate.c bitmap.c`

`neig.c` - A version of `iterate.c` for computing the neighbors of a ray. Used by
the controller program `runjob.pl`. Should be compiled for each base
size N separately. To compile for a given N, use
`gcc -O3 -o bin/neig<N> -D RANK=<N> neig.c bitmap.c`

`newrays.c` - Generates the extremal rays from the neighbor enumeration. It can
produce a canonical representative from the corresponding orbit
rather than the ray itself. Used by the controller program `runjob.pl`.
Sould be compiled for each base size N separately. To compile for a given N, use
`gcc -O3 -o bin/newr<N> -D RANK=<N> newrays.c`

`extract.c` - Generate extremal rays from an intermediate cone of the DD method
as produced by the controller program `base.pl`. It can produce a
canonical representative from the corresponding orbit rather than
the ray itself. Should be compiled for each base size N separately.
To compile for a given N, use
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

### Raw data

`figures-data` - Raw data for the figures and charts in the paper can be found
in the folder https://github.com/csirmaz/submodular-functions-6/tree/main/figures-data .
