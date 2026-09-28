#!/usr/bin/env perl6

use Math::NIntegrate;

#my %res;
#my @sample-points = do gather {
#    %res =
#            nintegrate(
#                    sub ($x, $y) { take ($x, $y); ($x + $y) },
#                    ('x', 0, 1), ('y', 0, {$_}),
#                    precision-goal => 8,
#                    #method => ('global-adaptive', method => ('GaussKronrodRule', points => 3), singularity-depth => 2, singularity-handler => 'none'),
#                    #method => ('global-adaptive', method => 'multi-dimensional-rule', singularity-depth => 2, singularity-handler => 'none'),
#                    max-recursion => 7):pairs
#}



my $n = 10;
my $precision-goal = 6;
my $min-recursion = 0;
my $max-recursion = 20;
my &f = { 1 / ($^x + $^y).sqrt };

my $tStart1 = now;
my %res1;
for ^$n {
    %res1 = nintegrate(
            &f,
            <x 0 1>, <y 0 1>,
            method => ('global-adaptive', method => ('gauss-kronrod-rule', points => 5), singularity-depth => Inf),
            :$max-recursion,
            :$precision-goal):pairs;
}
my $tEnd1 = now;
say "Integration time {$tEnd1 - $tStart1} seconds, {($tEnd1 - $tStart1) / $n} per call";

say %res1;

my $tStart2 = now;
my %res2;
for ^$n {
    %res2 = nintegrate(
            &f,
            <x 0 1>, <y 0 1>,
            method => ('global-adaptive', singularity-depth => Inf),
            :$max-recursion,
            :$precision-goal):pairs;
}
my $tEnd2 = now;
say "Integration time {$tEnd2 - $tStart2} seconds, {($tEnd2 - $tStart2) / $n} per call";

say %res2;
