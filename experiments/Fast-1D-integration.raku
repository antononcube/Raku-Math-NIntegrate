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
                integral => '1 / √x, x ∈ (0, 1)',
                function => -> $x {
                    1 / (1 + $x.sqrt)
                },
                domains => [<x 0 1>,],
            },
            {
                integral => '1 / (1 + x ** 2), x ∈ (-Inf, Inf)',
                function => -> $x {
                    1 / (1 + $x ** 2)
                },
                domains => [<x -Inf Inf>,],
            },
            {
                integral => '1 / (x ** 4 + x ** 2 + 0.9), x ∈ (-1, 1)',
                function => -> $x {
                    1 / ($x ** 4 + $x ** 2 + 0.9)
                },
                domains => [<x -1 1>, ],
            },
            {
                integral => '(sin(100) / pi) * $x, x ∈ (0.1, 1)',
                function => -> $x {
                    (sin(100) / pi) * $x
                },
                domains => [<x -.1 1>, ],
            },
            ;

    my @methods =
            {
                method => 'gauss-kronrod-rule',
                value => 'gauss-kronrod-rule',
            },
            {
                method => 'clenshaw-curtis-rule',
                value => 'clenshaw-curtis-rule',
            },
            {
                method => 'global-adaptive, singularit-depth => Inf',
                value => ('global-adaptive', singularity-depth => Inf),
            },
            {
                method => 'local-adaptive, singularit-depth => Inf',
                value => ('local-adaptive', singularity-depth => Inf),
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