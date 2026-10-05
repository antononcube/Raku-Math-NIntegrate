#!/usr/bin/env perl6

use Math::NIntegrate;

my $n = 10;
my $precision-goal = 3;
my $min-recursion = 0;
my $max-recursion = 20;
my $max-points = 10_000;
my $random-seed = 12231;
my &f = { ($^x + $^y).sqrt };


# With Rakupp we get the same list here!!
srand($random-seed);
#say 10.rand xx 10;

my $tStart1 = now;
my %res1;
for ^$n {
    %res1 = nintegrate(
            &f,
            <x 0 1>, <y 0 1>,
            method => ('monte-carlo', method => ('monte-carlo-rule', points => 1000), :$max-points),
            :$precision-goal):pairs;
}
my $tEnd1 = now;
say "Integration time {$tEnd1 - $tStart1} seconds, {($tEnd1 - $tStart1) / $n} per call";

say %res1;