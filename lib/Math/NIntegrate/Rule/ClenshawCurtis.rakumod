use v6.d;

use Math::NIntegrate::Rule::General;
use Math::Polynomial::Chebyshev;
use Data::Transformers;

class Math::NIntegrate::Rule::ClenshawCurtis
        is Math::NIntegrate::Rule::General {
    has $.points;

    submethod BUILD(UInt:D :$!points, :$working-precision = Num) {
        my %res = self.make-weights($!points, $working-precision);
        self.abscissas = |%res<abscissas>;
        self.weights = |%res<weights>;
        self.error-weights = |%res<error-weights>;
    }

    multi method new($points, $working-precision = Num) {
        self.bless(:$points, :$working-precision)
    }

    multi method new(:$points, :$working-precision = Num) {
        self.bless(:$points, :$working-precision)
    }

    method make-weights(UInt:D $points, $working-precision = Num) {

        die 'The first argument is expected to be a positive integer.' unless $points > 0;

        return { abscissas => [0], weights => [2], error-weights => [0] } if $points == 1;

        my $N = 2 * $points - 2;                    # interpolant degree
        my @x = (^($N + 1)).map({ cos($_ * π / $N) });

        sub moments($deg) {
            (0 .. $deg).map: -> $j {
                if $j % 2 { 0 }
                else { $j == 0 ?? 2 !! 2 / (1 - $j * $j) }
            }
        }

        sub weights($deg, @nodes, @a) {
            @nodes.keys.map: -> $k {
                my $s = [+] @a.keys.map: -> $j {
                    my $h = ($j == 0 || $j == $deg) ?? 1/2 !! 1;
                    $h * @a[$j] * chebyshev-t($j, @nodes[$k])
                }
                my $α = ($k == 0 || $k == $deg) ?? 1/2 !! 1;
                $α * (2 / $deg) * $s
            }
        }

        my @w = weights($N, @x, moments($N));

        my $M = $N div 2;                      # nested coarser degree
        my @x-c = (0 .. $M).map({ @x[2 * $_] });
        my @w-c = weights($M, @x-c, moments($M));
        my @e = @w.keys.map: -> $k {
            $k %% 2 ?? @w[$k] - @w-c[$k div 2] !! @w[$k]
        }

        @x .= reverse; @w .= reverse; @e .= reverse;

        return %(
            abscissas => rescale(@x, [-1, 1], [0, 1]),
            weights => @w <</>> 2,
            error-weights => @e <</>> 2)
    }
}
