use v6.d;

use Math::NIntegrate::Utilities;
use Math::NIntegrate::Rule::General;

class Math::NIntegrate::Rule::NewtonCotes
        is Math::NIntegrate::Rule::General {
    has $.points;
    has Bool:D $.open = False;

    submethod BUILD(UInt:D :$!points, :$wprec = Num, Bool:D :$!open = False) {
        my %res = self.make-weights($!points, $wprec, $!open);
        self.abscissas = |%res<abscissas>;
        self.weights = |%res<weights>;
        self.error-weights = |%res<error-weights>;
    }

    multi method new($points, $wprec = Num, Bool :$open = False) {
        self.bless(:$points, :$wprec, :$open)
    }

    multi method new(:$points, :$wprec = Num, Bool :$open = False) {
        self.bless(:$points, :$wprec, :$open)
    }

    method make-weights(UInt:D $points, $wprec = Num, Bool:D $open = False) {
        die 'The first argument is expected to be a positive integer.' unless $points > 0;

        # The fine rule is the rule used for integration. Its odd number of
        # nodes lets the lower-order rule use every other node, so that its
        # result can be subtracted without a second set of evaluations.
        my $count = 2 * $points - 1;
        my @abscissas = $open
                ?? (^$count).map({ numerical(($_ + 1) / ($count + 1), $wprec) }).Array
                !! ($count == 1
                        ?? [numerical(0, $wprec)]
                        !! (^$count).map({ numerical($_ / ($count - 1), $wprec) }).Array);

        my @weights = rule-weights(@abscissas, $wprec);
        my @seed-abscissas = (^$points).map({ @abscissas[2 * $_] }).Array;
        my @seed-weights = rule-weights(@seed-abscissas, $wprec);
        my @error-weights = @weights.clone;

        for ^$points -> $i {
            @error-weights[2 * $i] -= @seed-weights[$i];
        }

        return %(:@abscissas, :@weights, :@error-weights);
    }

    # Integrate the Lagrange basis polynomials exactly.
    # Keeping this calculation rational until the final conversion
    # avoids accumulating round-off when Rat or FatRat weights are specified.
    sub rule-weights(@nodes, $wprec) {
        @nodes.kv.map: -> $i, $xi {
            my @polynomial = 1;       # coefficients in ascending degree
            my $denominator = 1;

            for @nodes.kv -> $j, $xj {
                next if $i == $j;
                my @next = 0 xx (@polynomial.elems + 1);
                for @polynomial.kv -> $degree, $coefficient {
                    @next[$degree] += -$xj * $coefficient;
                    @next[$degree + 1] += $coefficient;
                }
                @polynomial = @next;
                $denominator *= $xi - $xj;
            }

            numerical(
                    ([+] @polynomial.kv.map: -> $degree, $coefficient {
                        $coefficient / $denominator / ($degree + 1)
                    }),
                    $wprec)
        }
    }
}