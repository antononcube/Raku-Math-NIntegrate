#!/usr/bin/env raku
use v6.d;

my %*SUB-MAIN-OPTS = :named-anywhere;

use Math::NIntegrate::Utilities;
use Math::NIntegrate::VariableTransformer;
use Math::NIntegrate::VariableTransformer::Infinite;
use Math::NIntegrate::Region;

use Data::Generators;
use Data::TypeSystem;


sub MAIN(
        Str:D $type = 'Infinite', #= One of: "Infinite" or "InfiniteEnBloc"
        Int:D :$repeat = 100,     #= Number of benchmarking computation runs,
        Int:D :$n = 1_000,        #= Number of points
        Int:D :$dim = 2,          #= Dimension of the points
        Int:D :$seed = 3234,      #= Random seed
        Bool:D :$echo = False,    #= Should intermediate results be printed out or not?
         ) {

    my $working-precision = Num;

    srand($seed);
    my @points = |(random-real([0, 1], $dim) xx $n)>>.Array;
    say deduce-type(@points) if $echo;

    my ($a1, $b1) = 10, Inf;
    my ($a2, $b2) = -Inf, Inf;

    my @result;
    my $tstart = now;
    given $type {

        when 'Infinite' {
            my $region = Math::NIntegrate::Region.new(min => [$a1, $a2], max => [$b1, $b2], dimension => 2);

            my @min-transform-bounds = 0 xx $region.dimension;
            my @max-transform-bounds = 1 xx $region.dimension;
            my @min-original-bounds = $a1, $a2;
            my @max-original-bounds = $b1, $b2;

            my %args =
                    region => $region,
                    working-precision => Rat,
                    ;

            my $vt = Math::NIntegrate::VariableTransformer::Infinite.new(
                    :@min-original-bounds, :@max-original-bounds, :@min-transform-bounds, :@max-transform-bounds,
                    |%args
                    );

            for ^$repeat {
                @result = [];
                for @points -> @point {
                    my $p = $vt.transform(:@point, jacobian => 1, :!functional-bounds);
                    @result.push($p<point>)
                }
            }
        }

        when 'InfiniteEnBloc' {
            # The initial version used the class Math::NIntegrate::VariableTransformer::AffineEnBloc
            # but that class was refactored and now Math::NIntegrate::VariableTransformer::Affine has
            # the method transform-en-bloc.
            my $region = Math::NIntegrate::Region.new(min => [$a1, $a2], max => [$b1, $b2], dimension => 2);

            my @min-transform-bounds = 0 xx $region.dimension;
            my @max-transform-bounds = 1 xx $region.dimension;
            my @min-original-bounds = $a1, $a2;
            my @max-original-bounds = $b1, $b2;

            my %args =
                    region => $region,
                    working-precision => Rat,
                    ;

            my $vt = Math::NIntegrate::VariableTransformer::Infinite.new(
                    :@min-original-bounds, :@max-original-bounds, :@min-transform-bounds, :@max-transform-bounds,
                    |%args
                    );

            for ^$repeat {
                @result = |$vt.transform-en-bloc(:@points, jacobian => 1, :!functional-bounds)<points>;
            }
        }

        default {
            die "Unknown type '$_'."
        }
    }
    my $tend = now;

    say "Total time : { $tend - $tstart }";
    say "Mean time  : { ($tend - $tstart) / $repeat }";

    if $echo { .say for @result }
}