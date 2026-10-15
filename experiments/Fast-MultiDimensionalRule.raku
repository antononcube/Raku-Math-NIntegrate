#!/usr/bin/env raku
use v6.d;

#use Math::NIntegrate:ver<0.0.8>;
use Math::NIntegrate;

my $repeat = 10;
my %res;
my $start = now;
for ^$repeat {
    #%res = nintegrate({ ($^x ** 2 + $^y ** 2 < 1).Int }, ['x', -1, 1], ['y', -1, 1], :4max-recursion, :4precision-goal):pairs;
    %res = nintegrate({ 1 / ($^x ** 2 + $^y ** 2) }, ['x', 0, 2], ['y', 0.1, 2], :12max-recursion, :6precision-goal, singularity-depth => Inf):pairs;
#    %res = nintegrate({ $^x ** 2 + $^y ** 5 }, ['x', 0, 2], ['y', 0, 2],
#            method => ('multidimensional-rule', points => 9),
#            :6max-recursion, :6precision-goal, singularity-depth => Inf):pairs;
}
my $end = now;
say %res;
say "total time: {$end - $start} s, per call: {($end - $start) / $repeat} s";
