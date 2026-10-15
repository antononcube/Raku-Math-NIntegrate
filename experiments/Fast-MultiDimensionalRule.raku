#!/usr/bin/env raku
use Math::NIntegrate;
use JSON::Fast;

my %*SUB-MAIN-OPTS = :named-anywhere;

sub MAIN(
        UInt:D $iterations = 10,
        Str:D :f(:$format)='json',
        Bool:D :$echo=False,

         ) {
    die 'Iteration count must be greater than zero' unless $iterations > 0;

    my @integrals =
            {
                integral => '$^ ** 2 + $^y ** 2 < 1, x ∈ (-1, 1), y ∈ (-1, 1)',
                function => -> $x, $y {
                    ($x ** 2 + $y ** 2 < 1).Int
                },
                domains => [<x -1 1>, <y -1 1>],
                max-recursion => 4,
                precision-goal => 4,
            },
            {
                integral => '1 / ($x ** 2 + $y ** 2), x ∈ (0, 2), y ∈ (0.1, 2)',
                function => -> $x, $y {
                    1 / ($x ** 2 + $y ** 2)
                },
                domains => [<x 0 2>, <y 0.1 2>],
                max-recursion => 12,
                precision-goal => 6,
            },
            {
                integral => '$x ** 2 + $y ** 2, x ∈ (0, 2), y ∈ (0, 2)',
                function => -> $x, $y {
                    $x ** 2 + $y ** 2
                },
                domains => [<x 0 2>, <y 0 2>],
                max-recursion => 6,
                precision-goal => 6,
            },
            ;

    my @methods =
            {
                method => 'multidimensional-rule',
                value => 'multidimensional-rule',
            },
            {
                method => 'gauss-kronrod-rule',
                value => 'gauss-kronrod-rule',
            },
            {
                method => "local-adaptive, method => 'multidimensional-rule'",
                value => ('local-adaptive', method => 'multidimensional-rule'),
            },
            {
                method => "local-adaptive, method => gauss-kronrod-rule",
                value => ('local-adaptive', method => 'gauss-kronrod-rule'),
            },
            {
                method => "unit-cube-rescaling",
                value => 'unit-cube-rescaling',
            },
            ;

    my @benchmarks;

    for @integrals -> %integral {

        say (%integral<integral>) if $echo;

        for @methods -> %method {
            next if %integral<domains>.elems == 1 && (%method<value> ~~ / multi /);

            say "\t{%method<value>}" if $echo;

            my $result;
            my $start = now;
            my @args = %integral<function>, |%integral<domains>;

            for ^$iterations {
                $result = nintegrate(|@args,
                        method => %method<value>,
                        max-recursion => %integral<max-recursion>,
                        precision-goal => %integral<precision-goal>,
                        singularity-depth => Inf):pairs;
            }

            my $total-time = now - $start;

            @benchmarks.push({
                integral => %integral<integral>,
                method => %method<method>,
                'total-time' => $total-time,
                'mean-time' => $total-time / $iterations,
                result => $result,
            });
        }
    }

    if $format.lc eq 'json' {
        say to-json(@benchmarks);
    } else {
        say @benchmarks.raku;
    }
}