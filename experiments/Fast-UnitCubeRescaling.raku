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
                integral => '1 / (1 + x ** 2), x ∈ (-Inf, Inf)',
                function => -> $x {
                    1 / (1 + $x ** 2)
                },
                domains => [<x -Inf Inf>,],
            },
            {
                integral => '1 / (x ** 2 + y ** 2), x ∈ (1, Inf), y ∈ (0, 1)',
                function => -> $x, $y {
                    1 / ($x ** 2 + $y ** 2)
                },
                domains => [<x 1 Inf>, <y 0 1>],
            },
            {
                integral => 'exp(-(x ** 2 + y ** 2)), x ∈ (-Inf, Inf), y ∈ (2, Inf)',
                function => -> $x, $y {
                    exp(-($x ** 2 + $y ** 2))
                },
                domains => [<x -Inf Inf>, <y 2 Inf>],
            },
            ;

    my @methods =
            {
                method => 'multidimensional-rule',
                value => 'multidimensional-rule',
            },
            {
                method => 'unit-cube-rescaling',
                value => 'unit-cube-rescaling',
            },
            {
                method => 'gauss-kronrod-rule',
                value => 'gauss-kronrod-rule',
            },
            {
                method => 'unit-cube-rescaling, method => gauss-kronrod-rule',
                value => ('unit-cube-rescaling', method => 'gauss-kronrod-rule'),
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
                $result = nintegrate(|@args, :method(%method<value>)):pairs;
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