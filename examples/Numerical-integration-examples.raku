#!/usr/bin/env perl6

use Math::NIntegrate;
use Text::Plot;

#----------------------------------------------------------
# Simple one dimensional integral
say '-' x 100;
say nintegrate( {$_.sqrt }, <x 0 1>);

#----------------------------------------------------------
# One dimensional integral with singularity
say '-' x 100;
say nintegrate( { 1 / $_.sqrt }, <x 0 1>, method => 'gauss-kronrod-rule');

#----------------------------------------------------------
# Two dimensional integral -- using :pairs in order to get the error estimate and number subregions used.
say '-' x 100;
say nintegrate( { $^x + $^y ** 2 }, <x 0 3>, <y 0 5>):pairs;

#----------------------------------------------------------
# Two dimensional integral with singularity.
# Obtain integral and error estimates and the sampling points used.
# Using specifications for precision goal, max-recursion, and method
say '-' x 100;
my %res;
my @sample-points = do gather {
    %res =
            nintegrate(
                    sub ($x, $y) { take ($x, $y); 1 / sqrt($x + $y) },
                    ('x', 0, 1), ('y', 0, 1),
                    precision-goal => 2,
                    #method => ('global-adaptive', method => 'multi-dimensional-rule', singularity-depth => 2, singularity-handler => 'none'),
                    max-recursion => 3):pairs
}

say (:%res);
say text-list-plot(@sample-points, :30height, :100width);