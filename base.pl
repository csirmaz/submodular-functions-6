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

# This is the controlling program for the DD method computing the
# extremal rays of the submodular functions for |X|=3,4,5,6.
# Creates the initial DD pair using linear algebra, and then calls
# the pipe program for every remaining inequality.

use strict;

sub usage {
    my ($msg) = @_;

    if ($msg) { print "$msg\n"; }
    print <<USAGE;
usage: base.pl <N> [options] <filestub>
   N         size of base set: 3,4,5,6
   --topt    tail-optimal inequality order (default)
   --rec     recursive inequality order
   --lexmin  lexmin inequality order with random shuffling of variables
   --dry     print axiom list but do not process
   --purge   delete intermediate results
   filestub  used for final and temporary files
If <filestub>-stop exists, stop processing
USAGE
    exit(0);
}


my $N         = 5;
my $VARNO     = 26;
my $AXNO      = 80;
my $helper    = "bin/iter5";
my $axfile    = "ax$N.c";
my $stopfile  = "stop";
my $ordertype = "topt";
my $purge     = 0;
my $dryrun    = 0;
my $filestub  = "";


sub read_args {
    if (scalar @ARGV <= 1) { usage(); }
    $N = $ARGV[0] + 0;
    $N == 3 || $N == 4 || $N == 5 || $N == 6 || usage("first argument should be 3,4,5,6");
    $helper = "bin/iter$N";
    $axfile = "ax$N.c";
    shift @ARGV;
    while (scalar @ARGV >= 1 && $ARGV[0] =~ /^--/) {
        if    ($ARGV[0] eq "--topt")   { $ordertype = "topt"; }
        elsif ($ARGV[0] eq "--rec")    { $ordertype = "recursive"; }
        elsif ($ARGV[0] eq "--lexmin") { $ordertype = "lexmin"; }
        elsif ($ARGV[0] eq "--dry")    { $dryrun = 1; }
        elsif ($ARGV[0] eq "--purge")  { $purge = 1; }
        else                           { usage("wrong option: $ARGV[0]"); }
        shift @ARGV;
    }
    scalar @ARGV <= 1 || usage("Additional arguments");
    scalar @ARGV == 1 || usage("Missing stub filename");
    $filestub = $ARGV[0];
    $stopfile = "$filestub-stop";
}


read_args();


# =====================================================================================
# reading and sorting axioms


sub bitset {    # 0, 2, 4, 5
    my $set = shift;

    my $v = 0;
    while ($set =~ s/(\d)//) {
        $v |= 1 << $1;
    }
    return $v;
}


sub varname {
    my $x = shift;

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
            next if ($b >= (1 << $N));
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


my @AX    = ();    # axiom coeffs as an array of integers
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
        if (scalar @coeffs != $VARNO) { die "wrong coeff number ($VARNO)\n$ax\n"; }
        push @AX,    \@coeffs;
        push @LABEL, parse_axname($txt, $axno);
        $axno++;
    }
    if ($axno != $AXNO) { die "axiom number $AXNO != $axno\n"; }
    close(FILE);
}


read_axioms();


# =================================================================
# determine the insertion order

my @ORDER = ();    # $ORDER[0] ... $ORDER[$AXNO-1] is the insertion order


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


sub get_axorder {    # tail-optimial order
    my ($force) = @_;

    return if (!$force && $ordertype ne "topt");
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
    my ($force) = @_;

    return if (!$force && $ordertype ne "recursive");

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


# -------------------------------------------------------------------------
# lexmin
#

my @VARPERM = ();    # permutations of variables


sub lexcmp {         # $ax[i] < $ax[j]
    my ($i, $j) = @_;

    for my $k (0 .. $VARNO - 1) {
        my $t = $AX[$i]->[$VARPERM[$k]] <=> $AX[$j]->[$VARPERM[$k]];
        return $t if ($t);
    }
    return 0;
}


sub get_lexminorder {
    my ($force) = @_;

    return if (!$force && $ordertype ne "lexmin");

    # generate a random order of variables
    for my $i (0 .. $VARNO - 1) { $VARPERM[$i] = $i; }
    for my $iidx (0 .. $VARNO - 1) {
        my $ii = $VARNO - 1 - $iidx;
        my $s  = int(rand($ii + 1));
        my $t  = $VARPERM[$ii];
        $VARPERM[$ii] = $VARPERM[$s];
        $VARPERM[$s]  = $t;
    }
    for my $i (0 .. $AXNO - 1) {
        $ORDER[$i] = $i;
    }

    # bubble sort
    my $doit = 1;
    while ($doit) {
        $doit = 0;
        for my $i (0 .. $AXNO - 2) {
            if (lexcmp($ORDER[$i], $ORDER[$i + 1]) < 0) {
                $doit = 1;
                my $t = $ORDER[$i];
                $ORDER[$i] = $ORDER[$i + 1];
                $ORDER[$i + 1] = $t;
            }
        }
    }
}


get_lexminorder();


# ------------------------------------------------------------------------
# for my $i(@ORDER){ print "$i: $LABEL[$i]\n"; }
## check that $ORDER is acually set

defined($ORDER[5]) || die "axiom order not defined\n";


# =============================================================================
#
# find a base: take inequalities as given by the permutation in $ORDER[]
#

sub is_zero {    # the argument is zero
    my $x = shift;

    return -4e-9 < $x && $x < 4e-9;
}


sub not_zero {    # the argument is not zero
    my $x = shift;

    return $x <= -4e-9 || $x >= 4e-9;
}


sub inner {       # inner product
    my ($ax, $ray) = @_;

    my $v = 0;
    for my $i (0 .. $VARNO - 1) { $v += $AX[$ax]->[$i] * $ray->[$i]; }
    return $v;
}


my @INBASE = ();    # define how inequalities are handled


# After processing by find_base()
#  1: it is in the base
#  0: to be processed later by the helper program

my @INITIAL_RAYS = ();    # initial rays are stored here


# Find a rank $VARNO subset of inequalities using the order given in $ORDER[]
sub find_base {
    @INBASE = (0) x $AXNO;
    my @cols       = (-1) x $VARNO;    # which columns are filled
    my $filled     = 0;
    my $nextlistno = 0;
    my @M;
    while ($filled < $VARNO) {
        if ($nextlistno == $AXNO) {
            die "Cannot generate base ($filled, $VARNO)\n";
        }
        my $nextax = $ORDER[$nextlistno];
        $nextlistno++;
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
        next if ($new < 0);    # all zero row
        $cols[$new] < 0 || die "next entry is not negative\n";
        $cols[$new] = $filled;
        my $pivot = $row[$new];
        for my $i (0 .. $VARNO - 1) { $row[$i] /= $pivot; }
        $row[$new] = 1;
        push @M, \@row;
        $filled++;
        $INBASE[$nextax] = 1;
    }
}


# --------------------------------------------------------------------
# integrify( $ray[] )
#  given a ray of size $VARNO, make all entries integer
#  this is quite tricky since the limited precision

sub lcm {    # 0< $v <1; find the smalles $d such that $d*$v is integer
    my ($v, $info) = @_;

    for my $dd (2 .. 250) {    ## exceptional values 792=8*9*11
        return $dd if (is_zero(($v * $dd - int($v * $dd + 2.4e-9)) / $info));
    }
    die "Strange value to make an integer from: v=$v, d=$info\n";
}


sub gcd {                      # x,y non-negative integers, fast Euclidean algorithm
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
sub store_initial_ray {
    my ($ray) = @_;

    # check that the ray is on the positive side of all axioms in the base
    for my $axn (0 .. $AXNO - 1) {
        next if ($INBASE[$axn] == 0);
        my $v = inner($axn, $ray);
        if ($v < -1e-10) { die "ray is not on axiom $axn (v=$v)\n"; }
    }
    push @INITIAL_RAYS, integrify($ray);
}


# invert the matrix consisting of the base vectors to generate the initial
# set of rays.
sub find_initial_rays {
    my @M   = ();
    my $idx = 0;
    for my $i (0 .. $AXNO - 1) {
        next if ($INBASE[$i] != 1);    # keep axioms with $INBASE[ax]==1
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
    print "inequalities to be processed in order:\n";
    for my $i (0 .. $AXNO - 1) {
        next if ($INBASE[$ORDER[$i]] != 1);
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

    scalar @INITIAL_RAYS == $VARNO || die "wrong number of initial rays\n";
    my $hcols = 0;
    open(FILE, ">", $file) || die "Cannot create file $file\n";
    print FILE "$VARNO";    # number of axioms handled so far
    for my $i (0 .. $AXNO - 1) {
        if ($INBASE[$i] == 1) { print FILE " $i"; $hcols++ }
    }
    print FILE "\n";
    $hcols == $VARNO || die "wrong number of handled axioms ($hcols,$VARNO)\n";

    # and the rays
    for my $r (@INITIAL_RAYS) {
        print FILE "", join(' ', @$r), "\n";
    }
    close(FILE);
}


# ==============================================================================

sub generate_solution {
    find_base();
    find_initial_rays();
    if ($dryrun) { print_axiom_order(); return ""; }
    my $itno = $VARNO;
    if (-e $stopfile) { unlink $stopfile; return "stopfile exists, aborting"; }
    my $total = $itno;
    for my $i (0 .. $AXNO - 1) { $total++ if ($INBASE[$i] == 0); }
    write_initial_rays("$filestub-$itno.txt");
    my $starttime = time();

    for my $i (0 .. $AXNO - 1) {
        next if ($INBASE[$ORDER[$i]] != 0);

        # check if we should quit
        if (-e $stopfile) {
            unlink $stopfile;
            $starttime = time() - $starttime;
            return "process stopped after $itno iterations, time=$starttime";
        }

        # call the helper program
        my $nextit = $itno + 1;
        print "Remaining: ", $total - $itno, ", axiom: ", $LABEL[$ORDER[$i]], ", elapsed=", time() - $starttime, "\n";
        system($helper, "-v", $itno, $ORDER[$i], "$filestub-$itno.txt", "$filestub-$nextit.txt", $stopfile);
        if (($? >> 8) != 0) {
            if (-e $stopfile) {
                unlink $stopfile;
                $starttime = time() - $starttime;
                return "iteration $itno is incomplete, total=$starttime";
            }
            else { die "Program $helper reported an error\n"; }
        }
        if ($purge) { unlink("$filestub-$itno.txt"); }
        $itno++;
    }
    $starttime = time() - $starttime;
    return "*** DONE, total time=$starttime ***";
}


# ===================================================================================

print "", generate_solution(), "\n";

__END__

