#!/usr/bin/perl -W

=pod

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

=cut

# Find the neighbors of a list of extremal rays

use strict;

sub usage {
    my ($msg) = @_;

    if ($msg) { print "$msg\n"; }
    print "Usage: runjob.pl <N> [options] <rayfile> <filestub>\n";
    print "  N         size of base set: 3,4,5,6\n";
    print "  --topt    use the tail-optimal axiom order (default)\n";
    print "  --rec     recursive order of axioms\n";
    print "  --dry     print list of axioms, but do not process\n";
    print "  --ray     generate neighboring rays, not orbits\n";
    print "  rayfile   rays to be probed\n";
    print "  filestub  result is appended to filesub-all.txt; processed rays to\n";
    print "            filestub-probed.txt, and the log to filestub-log.txt\n"; 
    print "Lines in the rayfile are |<label> <ray>| where label is alphanumeric, and\n";
    print "ray is a comma-separated list of integers (without space).\n";
    print "If the file \"<filestub>-skip\" exists, the ray currently processed is skipped;\n";
    print "if \"<filestub>-stop\" exists, processing is stopped after finishing\n";
    print "the current ray.\n";
    exit(0);
}


my $N     = 6;
my $VARNO = 57;
my $AXNO  = 240;

my $helper = "bin/neig6";

# helper program to execute a single step of the DD algorithm
#  $helper -v 1 <itno> <axiomno> <input> <output>
#  -v: verbose
#  itno: iteration number, also the number of axioms executed so far
#  axiomno: axiom to be added
#  input/output: input and output files

my $exgen = "bin/newr6";

# program generating the new extremal rays from the list of neighbors
# $exgen [--ray] <input> <output> <origray>
# --ray: save the ray, not the minimal representative from its orbit
# input/output: input, output files
# origray: comma separated list of original ray coefficients
#
my $axfile   = "ax6.c";
my $skipfile = "skip";    # skip file, must be this as hard-coded inot $helper
my $stopfile = "stop";    # stop after finishinig the present iteration

my $ordertype = "topt";   # default order type
my $genrays   = 0;        # generate rays, not orbits
my $dryrun    = 0;        # print axiom list only, don't do anything

my ($rayfile, $filestub);
my ($allfile, $logfile, $probedfile);


sub read_args {
    if (scalar @ARGV <= 1) { usage(); }
    $N = $ARGV[0] + 0;
    $N == 3 || $N == 4 || $N == 5 || $N == 6 || usage("first argument should be 3,4,5,6");
    shift @ARGV;
    $helper = "bin/neig$N";
    $exgen  = "bin/newr$N";
    $axfile = "ax$N.c";
    while (scalar @ARGV >= 1 && $ARGV[0] =~ /^--/) {
        if    ($ARGV[0] eq "--topt") { $ordertype = "topt"; shift @ARGV; }
        elsif ($ARGV[0] eq "--rec")  { $ordertype = "recursive"; shift @ARGV; }
        elsif ($ARGV[0] eq "--dry")  { $dryrun = 1; shift @ARGV; }
        elsif ($ARGV[0] eq "--ray")  { $genrays = 1; shift @ARGV; }
        else                         { usage("unknown option"); }
    }
    scalar @ARGV <= 2 || usage("Additional argument");
    scalar @ARGV == 2 || usage("missing argument");
    $rayfile    = $ARGV[0];
    $filestub   = $ARGV[1];
    $allfile    = "$filestub-all.txt";
    $logfile    = "$filestub-log.txt";
    $probedfile = "$filestub-probed.txt";
    $skipfile   = "$filestub-skip";
    $stopfile   = "$filestub-stop";

    if (!-e $rayfile) {
        print "ray file $rayfile not found, exiting\n";
        exit(1);
    }
}


read_args();


#=================================================================================
# reading and sorting axioms

sub bitset {    # 0, 2, 4, 5
    my $set = shift;

    my $v   = 0;
    while ($set =~ s/(\d)//) {
        $v |= 1 << $1;
    }
    return $v;
}


sub varname {
    my $x   = shift;

    my $txt = "";
    $txt .= "0" if ($x & 1);
    $txt .= "1" if ($x & 2);
    $txt .= "2" if ($x & 4);
    $txt .= "3" if ($x & 8);
    $txt .= "4" if ($x & 16);
    $txt .= "5" if ($x & 32);
    return $txt;
}


sub singlebit {
    my $x = shift;

    return $x == 1 || $x == 2 || $x == 4 || $x == 8 || $x == 16 || $x == 32;
}


sub parse_axname {
    my ($txt, $n) = @_;

    $txt =~ s/.*\#(\d+)\s+//;
    my $axno = $1;
    if ($axno != $n) { die "wrong number: needed: $n, got: $axno\n"; }
    if ($txt =~ /^f\(\{([\d,\s]+)\}\) \+ f\(\{([\d,\s]+)\}\)/) {
        my ($v1, $v2) = ($1, $2);
        $v1 = bitset($v1);
        $v2 = bitset($v2);    # A|C=v1, B|C=v2
        my $C = $v1 & $v2;
        my $A = $v1 & ~$C;
        my $B = $v2 & ~$C;
        if (!singlebit($A) || !singlebit($B)) {
            die "error in axname $txt\n";
        }
        my ($vA, $vB, $vC) = (varname($A), varname($B), varname($C));
        return $C == 0 ? "($vA,$vB)" : "($vA,$vB|$vC)";
    }
    if ($txt =~ /^f\([^\)]+\) >= f\(\{([\d,\s]+)\}\)/) {
        my $v = bitset($1);
        my $A = -1;
        my $B = -1;
        for my $b (1, 2, 4, 8, 16, 32) {
            next if (($b & $v) != 0);
            next if ($N == 5 && $b == 32);
            next if ($N == 4 && $b >= 16);
            next if ($N == 3 && $b >= 8);
            if    ($A < 0) { $A = $b; }
            elsif ($B < 0) { $B = $b; }
            else           { die "error2 in axname $txt, v=$v, A=$A, B=$B, b=$b\n"; }
        }
        ($A > 0 && $B > 0) || die "error3 in axname $txt\n";
        my ($vA, $vB, $vC) = (varname($A), varname($B), varname($v));
        return "($vA,$vB|$vC)";
    }
    die "cannot parse string [$txt]\n";
}


my @AX    = ();    # axiom coeffs as ann array of integers
my @LABEL = ();    # axiom labels in the order given


sub read_axioms {
    open(FILE, "<", $axfile) || die "Cannot open file $axfile\n";
    my $axno = 0;
    while (<FILE>) {
        if (/^\#define AXIOMS (\d+)/) { $AXNO  = $1; next; }
        if (/^\#define VARS (\d+)/)   { $VARNO = $1; next; }
        next if (!/^\{([01\-,]+)\}\s*(.*)\*/);
        my ($ax, $txt) = ($1, $2);
        my @coeffs = split(',', $ax);
        if (scalar @coeffs != $VARNO) { die "wrong coeff number\n$ax\n"; }
        push @AX,    \@coeffs;
        push @LABEL, parse_axname($txt, $axno);

        #        print parse_axname($txt,$axno),"\n";
        $axno++;
    }
    if ($axno != $AXNO) { die "axiom number=$axno\n"; }
    close(FILE);
}


read_axioms();


#===================================================================================
# put axiom order to @ORDER

my @ORDER = ();    # $ORDER[0] .. $ORDER[$AXNO-1] is the insertion order


#-------------------------------------------------------------------------------
# tail-optimal order
sub order_axlabel {
    my ($i, $j) = @_;

    return 0 if ($i == $j);
    my ($K1, $K2) = ("", "");
    if ($LABEL[$i] =~ /\|(.*)\)/) { $K1 = $1; }
    if ($LABEL[$j] =~ /\|(.*)\)/) { $K2 = $1; }
    my $O = $N == 3 ? [1, 0] : $N == 4 ? [2, 0, 1] : $N == 5 ? [3, 1, 0, 2] : [4, 2, 0, 1, 3];
    my ($o1, $o2) = ($O->[length($K1)], $O->[length($K2)]);
    if ($o1 != $o2) { return $o1 <=> $o2; }
    return $LABEL[$j] cmp $LABEL[$i];
}


sub get_axorder {    # simple basic order
    return if ($ordertype ne "topt");
    my @id = ();
    for my $i (0 .. $AXNO - 1) { $id[$i] = $i; }
    for my $i (sort { order_axlabel($a, $b) } @id) {
        push @ORDER, $i;
    }
}


get_axorder();       # basic order to $ORDER[]


#-------------------------------------------------------------------------------
# recursive axiom order
my %LISTED  = ();       # hash containing subsets encountered so far
my $LISTidx = $AXNO;    # index counting backward


# Subsets of [a,b,c,d,e]:
#  1) abcde
#  2) Subsets(abcd); Subsets(abce); Subsets(abde); Subsets(acde); Subsets(bcde)
# every subset is kept only at the first appearance.

sub Subsets {
    my ($set, $arg) = @_;

    my $str = join('', @$set);
    return if ($LISTED{"x$str"});
    $LISTED{"x$str"} = 1;
    useSubset($str, $arg);
    my $n = -1 + scalar @$set;
    if ($n < 0)  { return; }
    if ($n == 0) { Subsets([], $arg); return; }

    for my $i (0 .. $n) {
        my @subset = ();
        for my $j (0 .. $n) {
            next if ($i + $j == $n);
            push @subset, $set->[$j];
        }
        Subsets(\@subset, $arg);
    }
}


sub useSubset {
    my ($subset, $arg) = @_;

    # $arg is "$i,$j" $subset is either empty or a condition
    my $label = "($arg" . ($subset ne "" ? "|$subset)" : ")");
    my $found = -1;
    for my $i (0 .. $AXNO - 1) {
        if ($LABEL[$i] eq $label) { $found = $i; last; }
    }
    $found >= 0 || die "label $label not found\n";
    $LISTidx--;
    $ORDER[$LISTidx] = $found;
}


sub get_altorder {
    return if ($ordertype ne "recursive");

    # go over the pair ij, and then for all subsets
    $LISTidx = $AXNO;    # fill from the back
    for my $i (0 .. $N - 2) {
        for my $j ($i + 1 .. $N - 1) {
            my @rest = ();
            for my $k (0 .. $N - 1) {
                next if ($k == $i || $k == $j);
                push @rest, $k;
            }

            # now we have ij, and the rest in array @rest
            %LISTED = ();    # nothing is there
            Subsets(\@rest, "$i,$j");
        }
    }
    $LISTidx == 0 || die "LISTidx=$LISTidx not zero\n";
    %LISTED = ();            # cleanup
}


get_altorder();


# ------------------------------------------------------------------------
# for my $i(@ORDER){ print "$i: $LABEL[$i]\n"; }
## check that $ORDER is acually set

defined($ORDER[5]) || die "axiom order not defined\n";


#===================================================================================
# find a base: take axioms as given by the permutation in $ORDER[]
#

sub is_zero {
    my $x = shift;

    return -4e-9 < $x && $x < 4e-9;
}


sub not_zero {
    my $x = shift;

    return $x <= -4e-9 || $x >= 4e-9;
}


sub inner {    # inner product
    my ($ax, $ray) = @_;

    my $v = 0;
    for my $i (0 .. $VARNO - 1) { $v += $AX[$ax]->[$i] * $ray->[$i]; }
    return $v;
}


my @INBASE = ();    # define how axioms are handled


# Before calling find_base()
#  0: axioms which lie on the ray to be handled
# -1: others; since the rank of the zeroset is $VARNO-1, at least one
#         such axiom exists
# After processing by find_base()
#  1: axioms which are kept to be zero
#  0: axioms to be processed by the helper program
#  2: axioms processed so far
# -1: axioms to be ignored

my @INITIAL_RAYS = ();    # initial rays are stored here


# using $INBASE[] find a rank $VARNO subset of axioms using
#  1) the order given in $ORDER[];
#  3) as many axioms with $INBASE[] == 0 as possible     ==> 2
#  3) a single axiom with $INBASE[] == -1                ==> 1

sub find_base {
    my @cols       = (-1) x $VARNO;    # which columns are filled
    my $filled     = 0;
    my $nextlistno = 0;
    my @M;
    my $scan = 0;
    while ($filled < $VARNO) {
        if ($nextlistno == $AXNO) {
            if ($scan < 2) { $scan++; $nextlistno = 0; next; }
            die "Cannot generate base ($filled, $VARNO)\n";
        }
        my $nextax = $ORDER[$nextlistno];
        $nextlistno++;
        if    ($scan == 0) { next if ($INBASE[$nextax] != 1); }
        elsif ($scan == 1) { next if ($INBASE[$nextax] != 0); }
        else               { next if ($INBASE[$nextax] != -1); }
        my @row = ();
        for my $i (0 .. $VARNO - 1) { $row[$i] = $AX[$nextax]->[$i]; }

        for my $i (0 .. $VARNO - 1) {
            next if ($cols[$i] < 0 || is_zero($row[$i]));
            my $pivot = $row[$i];
            my $idx   = $cols[$i];
            for my $ii (0 .. $VARNO - 1) { $row[$ii] -= $M[$idx]->[$ii] * $pivot; }
            is_zero($row[$i]) || die "not zero row elem\n";
            $row[$i] = 0;
        }

        # find a non-zero element in $row if possible
        my $new = -1;
        for my $i (0 .. $VARNO - 1) {
            if (not_zero($row[$i])) { $new = $i; last; }
        }
        if ($new < 0) {
            next if ($scan > 0);
            die "Cannot insert axiom $nextax to base\n";
        }
        $cols[$new] < 0 || die "next entry is not negative\n";
        $cols[$new] = $filled;
        my $pivot = $row[$new];
        for my $i (0 .. $VARNO - 1) { $row[$i] /= $pivot; }
        $row[$new] = 1;
        push @M, \@row;
        $filled++;
        $INBASE[$nextax] += 2;
        if ($scan == 2 && $filled < $VARNO) { die "more additional axiom requested\n"; }
    }
}


# --------------------------------------------------------------------
# integrify( $ray[] )
#  given a ray of size $VARNO, make all entries integer
# this is quite tricky since the limited precision
sub lcm {    # 0< $v <1; find the smalles $d such that $d*$v is integer
    my ($v, $info) = @_;

    for my $dd (2 .. 310, 315, 327, 332, 338, 365) {    ## exceptional values 792=8*9*11
        return $dd if (is_zero(($v * $dd - int($v * $dd + 2.4e-9)) / $info));
    }
    die "Strange value to make an integer from: v=$v, d=$info\n";
}


sub gcd {                                               # x,y non-negative integers, fast Euclidean algorithm
    my ($x, $y) = @_;

    if ($x == 1) { return 1; }
    if ($x == 0) { return $y; }
    while (1) {
        $y %= $x;
        if ($y == 0) { return $x; }
        if ($y == 1) { return 1; }
        $x %= $y;
        if ($x == 0) { return $y; }
        if ($x == 1) { return 1; }
    }
}


# given a ray, make an integer multiple; also divide all coordinates
#   by the lcm at the last step. Allow negative coordinates.
sub integrify {
    my ($ray) = @_;

    my @iray  = (0) x $VARNO;
    my @signs = (0) x $VARNO;
    for my $i (0 .. $VARNO - 1) {
        next if ($ray->[$i] >= 0);
        $signs[$i] = 1;
        $ray->[$i] = -$ray->[$i];
    }
    my $d = 1;
    for my $i (0 .. $VARNO - 1) {
        my $v  = $d * $ray->[$i];
        my $iv = int((1 + 5e-10) * $v + 2e-9 * $d);
        $iray[$i] = $iv;
        next if (is_zero(($v - $iv) / $d));
        next if ($v > 10.0 && is_zero(($v - $iv) / $v));
        if ($v - $iv < 0) { die "i=$i, v=$v, iv=$iv\n"; }
        $d = $d * lcm($v - $iv, $d);
    }
    if ($d != 1) {
        for my $i (0 .. $VARNO - 1) {
            my $v  = $d * $ray->[$i];
            my $iv = int($v + 0.01);
            $iray[$i] = $iv;
            is_zero(($v - $iv) / $d) || die "integrify problem: v=$v, iv=$iv\n";
        }
    }
    $d = gcd($iray[0], $iray[1]);
    my $idx = 2;
    while ($d != 1 && $idx < $VARNO) { $d = gcd($iray[$idx], $d); $idx++; }
    if ($d > 1) {
        for my $i (0 .. $VARNO - 1) { $iray[$i] /= $d; }
    }
    for my $i (0 .. $VARNO - 1) {
        if ($iray[$i] != 0 && $signs[$i] != 0) { $iray[$i] = -$iray[$i]; }
    }
    return \@iray;
}


# store a ray in $INITIAL_RAYS[] after integrifying.
#   check that the ray is consistent with $INBASE[] (just a safety check)
#   if $INBASE[ax]==1 then ax*ray == 0
#   if $INBASE[ax]==2 then ax*ray >= 0
sub store_initial_ray {
    my ($ray) = @_;

    # check that the ray is on the positive side of all axioms in the base
    for my $axn (0 .. $AXNO - 1) {
        next if ($INBASE[$axn] <= 0);
        my $v = inner($axn, $ray);
        if ($v < -1e-10)                      { die "ray is not on axiom $axn (v=$v)\n"; }
        if ($INBASE[$axn] == 1 && $v > 1e-10) { die "ray is not on axiom $axn (v=$v)\n"; }
    }
    push @INITIAL_RAYS, integrify($ray);
}


# invert the matrix consisting of the base vectors to generate the initial
# set of rays. Solutions corresponding to axioms violations with $INBASE[ax]==1
#   are skipped
sub find_initial_rays {
    my @M       = ();
    my $idx     = 0;
    my @skipvar = (0) x $VARNO;
    for my $i (0 .. $AXNO - 1) {
        next if ($INBASE[$i] <= 0);    # keep axioms with $INBASE[ax]== 1 or 2
        $skipvar[$idx] = 1 if ($INBASE[$i] == 1);
        my @row = ();
        for my $k (0 .. $VARNO - 1) { $row[$k] = $AX[$i]->[$k]; }
        for my $k (0 .. $VARNO - 1) { $row[$VARNO + $k] = $k == $idx ? 1 : 0; }
        push @M, \@row;
        $idx++;
    }
    $idx == $VARNO || die "matrix row number is $idx\n";

    # eliminate row $row
    for my $row (0 .. $VARNO - 1) {
        my $pcol = -1;
        for my $j (0 .. $VARNO - 1) {
            if (not_zero($M[$row]->[$j])) { $pcol = $j; last; }
        }
        $pcol >= 0 || die "The matrix is degenerate\n";
        my $pivot = 1.0 / $M[$row]->[$pcol];
        for my $i (0 .. $VARNO + $VARNO - 1) {
            $M[$row]->[$i] *= $pivot;
        }
        $M[$row]->[$pcol] = 1;

        # eliminate column $pcol
        for my $j (0 .. $VARNO - 1) {
            next if ($j == $row);
            my $v = $M[$j]->[$pcol];
            for my $i (0 .. $VARNO + $VARNO - 1) { $M[$j]->[$i] -= $v * $M[$row]->[$i]; }
            is_zero($M[$j]->[$pcol]) || die "Numerical instability\n";
            $M[$j]->[$pcol] = 0;
        }
    }

    # now the columns are the rays
    @INITIAL_RAYS = ();
    my $ray = [];
    for my $i (0 .. $VARNO - 1) {
        next if ($skipvar[$i]);
        for my $row (0 .. $VARNO - 1) {
            my $idx = -1;
            for my $s (0 .. $VARNO - 1) {
                if (not_zero($M[$row]->[$s])) { $idx = $s; last; }
            }
            $idx >= 0 || die "index not found\n";
            $ray->[$idx] = $M[$row]->[$VARNO + $i];
        }
        store_initial_ray($ray);
    }
}


sub print_axiom_order {
    print "axioms to be processed in order:\n";
    for my $i (0 .. $AXNO - 1) {
        next if ($INBASE[$ORDER[$i]] != 2);
        printf "* %3d => $LABEL[$ORDER[$i]]\n", $ORDER[$i];
    }
    for my $i (0 .. $AXNO - 1) {
        next if ($INBASE[$ORDER[$i]] != 0);
        printf "  %3d => $LABEL[$ORDER[$i]]\n", $ORDER[$i];
    }
}


# create the initial ray file for the helper program
sub write_initial_rays {
    my ($file) = @_;

    open(FILE, ">", $file) || die "Cannot create file $file\n";
    my $extracol = 0;    # how many axiioms are fixed; those with $INBASE[ax]==1
    for my $i (0 .. $AXNO - 1) {
        if ($INBASE[$i] == 1) { $extracol++; }
    }
    $extracol == 1 || die "no extra axiom is found ($extracol)\n";
    print FILE "$extracol";
    for my $i (0 .. $AXNO - 1) {
        if ($INBASE[$i] == 1) { print FILE " $i"; }
    }
    print FILE "\n", $VARNO - $extracol;    # number of axioms handled so far
    scalar @INITIAL_RAYS == $VARNO - $extracol || die "wrong number of initial rays\n";
    for my $i (0 .. $AXNO - 1) {            # axioms handled so far
        next if ($INBASE[$i] != 2);
        print FILE " $i";
        $extracol++;
    }
    print FILE "\n";
    $extracol == $VARNO || die "wrong number of handled axioms ($extracol,$VARNO)\n";

    # and the rays
    for my $r (@INITIAL_RAYS) {
        print FILE "", join(' ', @$r), "\n";
    }
    close(FILE);
}


# ================================================================================
# given the extrmal ray $RAY[], select all axioms it satisfies.
# It must have rank $VARNO-1. Then find the first axiom to be added
#  which makes the rank $VARNO.
# Prepare the iteration. The initial ray set has $VARNO-1 extremal rays
#  from i-th = 1, all others = 0 (where i runs over axioms in r_A )
#  <n> <aximos handled so far (n of them)>
#  <all extremal rays found so far>
# Calling the program gives <n>, <next axiom>, <ifile>, <ofile>
#   the helper program executes a single iteration
#   old axioms: ax-1,...,ax-k; rank=varno-m

# call the helper program to determine all neighbors of the given ray
#  make some rudimentary checks
#  filenames are derived from $filestub
#  clean up after generating all extremal rays
# return the file name for the final result (to be unlinked)
#  or the empty string if the ray is to be skipped

sub generate_neighbor_file {
    my ($ray, $filestub) = @_;

    my $zerono = 1;    # number of axioms kept to zero
    for my $axn (0 .. $AXNO - 1) {
        my $v = inner($axn, $ray);
        if ($v < 0) { die "original ray is in the negative side of axiom $axn\n"; }
        elsif ($v == 0) { $INBASE[$axn] = 0; }
        else            { $INBASE[$axn] = -1; }
    }
    find_base();       # we expect one excluded axiom (rank == VARNO-1)
    find_initial_rays();
    if ($dryrun) { print_axiom_order(); return ""; }
    my $itno = $VARNO - $zerono;    # axioms handled so far
    write_initial_rays("$filestub-$itno.txt");

    # now we call the helper program for all axioms with $INBASE[ax]==0
    if (-e $skipfile) { unlink $skipfile; }
    my $total = $itno;
    for my $i (0 .. $AXNO - 1) { $total++ if ($INBASE[$i] == 0); }
    for my $i (0 .. $AXNO - 1) {
        next if ($INBASE[$ORDER[$i]] != 0);

        # check if we should quit
        if (-e $skipfile) {
            unlink $skipfile;
            for my $fn ($VARNO - $zerono .. $itno) { unlink "$filestub-$fn.txt"; }
            return "";    # indicate failure
        }

        # call the helper program
        my $nextit = $itno + 1;
        print "Remaining: ", $total - $itno, ", axiom: ", $LABEL[$ORDER[$i]], "\n";
        system($helper, "-v", $itno, $ORDER[$i], "$filestub-$itno.txt", "$filestub-$nextit.txt", $skipfile);
        if (($? >> 8) != 0) {
            if (-e $skipfile) {
                unlink $skipfile;
                for my $fn ($VARNO - $zerono .. $nextit) { unlink "$filestub-$fn.txt"; }
                return "";    # indicate failure
            }
            else {
                die "Program $helper reported an error\n";
            }
        }
        $itno++;
    }

    # remove all but the last files
    for my $fn ($VARNO - $zerono .. $itno - 1) { unlink "$filestub-$fn.txt"; }
    return "$filestub-$itno.txt";
}


#=================================================================================
# how many symmetries a ray has
# also compute the dual of the ray
sub hw {    # Hamming weight
    my $x = shift;

    if ($x < 0) { $x = -$x; }
    my $w = 0;
    while ($x) {
        if ($x & 1) { $w++; }
        $x >>= 1;
    }
    return $w;
}


sub varidxtostr {    # idx to var representation
    my $x = shift;

    if (hw($x) > $N - 1) { $x = -1 + (1 << $N); }
    my $txt  = "";
    my $base = 0;
    while ($x) {
        if ($x & 1) { $txt .= "$base"; }
        $base++;
        $x >>= 1;
    }
    return $txt;
}


sub varidxtocompl {    # 1,2,3  => complement of representation
    my $x = shift;

    if (hw($x) >= $N - 1) { $x = -1 + (1 << $N); }
    my $txt  = "";
    my $base = 0;
    $x ^= -1 + (1 << $N);
    while ($x) {
        if ($x & 1) { $txt .= "$base"; }
        $base++;
        $x >>= 1;
    }
    return $txt;
}


my @varnamearr      = ();
my @vidxarr         = ();
my @complvarnamearr = ();


sub vname {    # e.g, 24 => "034"
    return $varnamearr[shift];
}


sub vcompl {    # complement of a varname
    return $complvarnamearr[shift];
}


sub vidx {      # e.g., "034" => 24
    my $txt = shift;

    my $v   = 0;
    my $d   = 1;
    for my $i (0 .. 5) {
        $v |= $d if ($txt =~ /$i/);
        $d <<= 1;
    }
    return $vidxarr[$v];
}


sub varshw {    # order of variables when ordered by Hamming weight
    my @varw = ();
    my $was  = 0;
    for my $i (1 .. -1 + (1 << $N)) {
        if (hw($i) >= $N - 1) { next if ($was); $was = 1; $vidxarr[-1 + (1 << $N)] = scalar @varw; }
        $vidxarr[$i] = scalar @varw;
        push @varw,            hw($i);
        push @varnamearr,      varidxtostr($i);
        push @complvarnamearr, varidxtocompl($i);
    }
    for my $i (0 .. -1 + (1 << $N)) {
        if (!defined $vidxarr[$i]) { $vidxarr[$i] = $vidxarr[-1 + (1 << $N)]; }
    }
    scalar @varnamearr == $VARNO || die "variable number error\n";
}


varshw();


# determine variable permutations generated by permuting 0,1,2,3,4,5 (6!)
sub swap01 {
    my $txt = shift;

    $txt =~ s/0/A/g;
    $txt =~ s/1/0/g;
    $txt =~ s/A/1/g;
    return $txt;
}


sub shift05 {
    my $txt = shift;

    $txt =~ s/0/A/g;
    my $Nm1 = $N - 1;
    my $f   = 1;
    my $t   = 0;
    for (1 .. $Nm1) { $txt =~ s/$f/$t/g; $f++; $t++; }
    $txt =~ s/A/$Nm1/g;
    return $txt;
}


my @VPERMS  = ();    # variable permutations
my $VPERMNO = 1;     # N factorial

for my $i (2 .. $N) { $VPERMNO *= $i; }


sub addvperm_ifnew {
    my ($p) = @_;

    foreach my $vp (@VPERMS) {
        my $same = 1;
        for my $i (0 .. -1 + $VARNO) {
            if ($p->[$i] != $vp->[$i]) { $same = 0; last; }
        }
        if ($same) { return 0; }
    }
    push @VPERMS, $p;
    return 1;
}


sub generate_varperms {
    my @id = ();
    for my $i (0 .. -1 + $VARNO) { $id[$i] = $i; }
    push @VPERMS, \@id;
    my @gen1 = ();    # swap 0 and 1
    for my $i (0 .. -1 + $VARNO) { $gen1[$i] = vidx(swap01(vname($i))); }
    my @gen2 = ();    # shift 012345 => 123450
    for my $i (0 .. -1 + $VARNO) { $gen2[$i] = vidx(shift05(vname($i))); }
    for (my $vidx = 0; $vidx < scalar @VPERMS; $vidx++) {
        my @new1 = ();
        my $old  = $VPERMS[$vidx];
        for my $i (0 .. -1 + $VARNO) { $new1[$i] = $gen1[$old->[$i]]; }
        addvperm_ifnew(\@new1);
        my @new2 = ();
        for my $i (0 .. -1 + $VARNO) { $new2[$i] = $gen2[$old->[$i]]; }
        addvperm_ifnew(\@new2);
    }
    scalar @VPERMS == $VPERMNO || die "varperms generation error\n";
}


generate_varperms();


# compute the dual of a ray
#  r^*(A)= -r(N) + r(N-A)+\sum{ r(a): a \in A }
sub dual {
    my ($ray) = @_;

    my $full = "";
    for my $i (0 .. $N - 1) { $full .= "$i"; }
    my @dual = ();
    my $K    = $ray->[vidx($full)];
    for my $j (0 .. -1 + $VARNO) {
        my $cmp = vcompl($j);
        $dual[$j] = -$K;
        if ($cmp eq "") {
            for my $i (0 .. $N - 1) { $dual[$j] += $ray->[vidx("$i")]; }
        }
        elsif (length($cmp) == -1 + $VARNO) {
            $dual[$j] = $ray->[$j];
        }
        else {
            $dual[$j] += $ray->[vidx($cmp)];
            for my $i (0 .. $N - 1) {
                next if ($cmp =~ /$i/);
                $dual[$j] += $ray->[vidx("$i")];
            }
        }
    }
    return \@dual;
}


sub symmno {    # how many symmetries the ray has
    my ($ray) = @_;

    my $symm  = 0;
    my $base  = join(',', @$ray);
    my @new   = ();
    my $dual  = dual($ray);
    for my $p (@VPERMS) {
        for my $i (0 .. $VARNO - 1) { $new[$i] = $ray->[$p->[$i]]; }
        $symm++ if ($base eq join(',', @new));
        for my $i (0 .. $VARNO - 1) { $new[$i] = $dual->[$p->[$i]]; }
        $symm++ if ($base eq join(',', @new));
    }
    return $symm;
}


sub raydegree {    # number of axioms it satisfies
    my ($ray) = @_;

    my $degree = 0;
    for my $i (0 .. $AXNO - 1) {
        $degree++ if (inner($i, $ray) == 0);
    }
    return $degree;
}


#=================================================================================
sub append_to {
    my ($file, $line) = @_;

    open(FILE, ">>", $file) || die "Cannot append to $file\n";
    print FILE "$line\n";
    close(FILE);
}


# handle a ray with label $label
sub handle_ray {
    my ($label, $ray, $origray) = @_;

    my @r = split(',', $ray);
    scalar @r == $VARNO || die "wrong dimension in ray $ray, label=$label\n";
    my $logrd     = raydegree(\@r);
    my $logsym    = symmno(\@r);
    my $starttime = time();
    my $resfile   = generate_neighbor_file(\@r, $filestub);
    if (!$resfile) {
        if ($dryrun) { return "--- dryrun $label,axdim=$logrd,sym=$logsym"; }
        return "+++ skipped $label,axdim=$logrd,sym=$logsym";
    }
    my $stored;
    if   ($genrays) { $stored = `$exgen --ray $resfile $allfile $ray`; }
    else            { $stored = `$exgen $resfile $allfile $ray`; }
    print $stored;
    ($? >> 8) == 0 || die "program $exgen exited with an error\n";
    my $rno = 0;
    if ($stored =~ /(\d+) rays appended/) { $rno = $1; }
    my $logtime = time() - $starttime;
    unlink($resfile);
    append_to($probedfile, "#$label $logrd: $origray");
    append_to($logfile,    "#$label,axn=$logrd,degree=$rno,sym=$logsym,time=$logtime,ord=$ordertype");
    return "*** done $label,axdim=$logrd,degree=$rno,sym=$logsym,time=$logtime";
}


#================================================================================

my %toskip = ();    # $toskip{#label}=1 if it is to be skipped


sub read_logfile {
    open(LOG, "<", $logfile) || return;
    while (<LOG>) {
        /^#([\w]+),/ || die "wrong line in $logfile:\n$_\n";
        $toskip{$1} = 1;
    }
    close(LOG);
}


sub process_rayfile {
    my @RAYS      = ();
    my @RAYLABELS = ();
    my $total     = 0;
    read_logfile();
    open(RAYFILE, "<", $rayfile) || die "Cannot open rayfile $rayfile\n";
    while (<RAYFILE>) {
        my ($label, $ray);
        if (/^\#?(\w+)\s+\d+:\s+([\d,]+)$/) {
            $label = $1;
            $ray = $2;
        }
        elsif (/^\#?(\w+):?\s+([\d,]+)$/) {
            $label = $1;
            $ray   = $2;
        }
        elsif (m#^//#) {
            last;
        }
        else {
            die "wrong line in $rayfile\n $_\n";
        }
        next if ($toskip{$label});
        push @RAYS,      $ray;
        push @RAYLABELS, $label;
        $total++;
    }
    close(RAYFILE);
    print "Number of rays to process: $total\n";
    unlink($stopfile);
    my $newrays     = 0;
    my $skippedrays = 0;
    for my $idx (0 .. -1 + $total) {
        print "*** doing $RAYLABELS[$idx] (", $idx + 1, " of $total)\n";
        my $res = "";
        $res = handle_ray($RAYLABELS[$idx], $RAYS[$idx], $RAYS[$idx]);
        if ($res =~ /degree=(\d+),/) { $newrays += $1; }
        else                         { $skippedrays++; }
        print "$res\n";
        if (-e $stopfile) {
            unlink $stopfile;
            print "+++ Processing stopped after ", $idx + 1, " rays\n";
            print "    total=$newrays rays added to $allfile\n";
            return;
        }
    }
    if ($skippedrays == 0) {
        print "DONE, processed=$total,added=$newrays,file=$allfile\n";
    }
    else {
        $total -= $skippedrays;
        print "DONE, processed=$total,skipped=$skippedrays,added=$newrays,file=$allfile\n";
    }
}


process_rayfile();
exit 0;


__END__

