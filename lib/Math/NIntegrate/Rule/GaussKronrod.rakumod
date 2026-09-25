use v6.d;

use Math::NIntegrate::Utilities;
use Math::NIntegrate::Rule::General;

class Math::NIntegrate::Rule::GaussKronrod
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

    multi method new(:$points, :prec(:$working-precision) = Num) {
        self.bless(:$points, :$working-precision)
    }

    method make-weights(UInt:D $points, $working-precision = Num) {
        die 'The first argument is expected to be a positive integer.' unless $points > 0;

        # All construction is done on [-1, 1]. The returned rule is then mapped to [0, 1].

        # Coefficients of P_n in increasing-power order.
        my @legendre = $points == 1 ?? [0e0, 1e0] !! do {
            my @previous = [1e0];
            my @current = [0e0, 1e0];
            for 2 .. $points -> $degree {
                my @next = 0e0 xx ($degree + 1);
                for @current.kv -> $index, $coefficient {
                    @next[$index + 1] += (2 * $degree - 1) * $coefficient / $degree;
                }
                for @previous.kv -> $index, $coefficient {
                    @next[$index] -= ($degree - 1) * $coefficient / $degree;
                }
                @previous = @current;
                @current = @next;
            }
            @current
        };

        # Newton iteration for the n roots of P_n.  The cosine estimates are
        # already in their respective root basins.
        my @gauss = (^$points).map: -> $index {
            my $x = cos(π * ($index + 0.75) / ($points + 0.5));
            for ^40 {
                my ($pn, $pn-minus-one) = legendre-values($points, $x);
                my $derivative = $points * ($x * $pn - $pn-minus-one) / ($x * $x - 1e0);
                my $next = $x - $pn / $derivative;
                last if ($next - $x).abs < 2e-15;
                $x = $next;
            }
            $x
        };
        @gauss .= sort;

        # The monic Stieltjes polynomial E_(n+1) is characterized by
        # integral(P_n(x) E_(n+1)(x) x^k) = 0, k = 0 .. n.
        # Its roots are exactly the n + 1 Kronrod extension nodes and interlace the roots of P_n.
        my @stieltjes-matrix;
        my @stieltjes-rhs;
        for 0 .. $points -> $row {
            @stieltjes-matrix.push: (0 .. $points).map: -> $column {
                [+] @legendre.kv.map: -> $degree, $coefficient {
                    $coefficient * moment($degree + $column + $row)
                }
            };
            @stieltjes-rhs.push: -([+] @legendre.kv.map: -> $degree, $coefficient {
                $coefficient * moment($degree + $points + 1 + $row)
            });
        }
        my @extension-polynomial = |solve(@stieltjes-matrix, @stieltjes-rhs), 1e0;

        # There is one extension root in each interval delimited by the
        # Gauss roots.  Bisection exploits that interlacing and avoids a
        # second, less reliable, general polynomial-root solver.
        my @boundaries = -1e0, |@gauss, 1e0;
        my @extension;
        for ^($points + 1) -> $interval {
            my $left = @boundaries[$interval];
            my $right = @boundaries[$interval + 1];
            my $left-value = polynomial-value(@extension-polynomial, $left);
            for ^80 {
                my $middle = ($left + $right) / 2e0;
                my $middle-value = polynomial-value(@extension-polynomial, $middle);
                if $left-value * $middle-value <= 0 {
                    $right = $middle;
                } else {
                    $left = $middle;
                    $left-value = $middle-value;
                }
            }
            @extension.push(($left + $right) / 2e0);
        }

        my @nodes = (|@gauss, |@extension).sort;
        my $total-points = 2 * $points + 1;

        # Integrate each Lagrange cardinal polynomial.
        # This is considerably better conditioned than solving a monomial Vandermonde system for
        # the weights, while still giving the unique interpolatory rule.
        my @weights = @nodes.kv.map: -> $index, $node {
            my @cardinal = [1e0];
            my $denominator = 1e0;
            for @nodes.kv -> $other-index, $other-node {
                next if $other-index == $index;
                my @next = 0e0 xx (@cardinal.elems + 1);
                for @cardinal.kv -> $degree, $coefficient {
                    @next[$degree] -= $other-node * $coefficient;
                    @next[$degree + 1] += $coefficient;
                }
                @cardinal = @next;
                $denominator *= $node - $other-node;
            }
            ([+] @cardinal.kv.map: -> $degree, $coefficient {
                $coefficient * moment($degree) / $denominator
            })
        };

        my @error-weights = @weights.Array;
        for @gauss -> $node {
            my $index = @nodes.first({ ($_ - $node).abs < 1e-12 }, :k);
            my ($pn, $pn-minus-one) = legendre-values($points, $node);
            my $derivative = $points * $pn-minus-one / (1e0 - $node * $node);
            my $gauss-weight = 2e0 / ((1e0 - $node * $node) * $derivative * $derivative);
            @error-weights[$index] -= $gauss-weight;
        }

        return %(
            abscissas => @nodes.map({ numerical(($_ + 1e0) / 2e0, $working-precision) }).Array,
            weights => @weights.map({ numerical($_ / 2e0, $working-precision) }).Array,
            error-weights => @error-weights.map({ numerical($_ / 2e0, $working-precision) }).Array);
    }

    sub moment(Int:D $degree) {
        $degree %% 2 ?? 2e0 / ($degree + 1) !! 0e0
    }

    sub polynomial-value(@coefficients, $x) {
        my $value = 0e0;
        for @coefficients.reverse -> $coefficient {
            $value = $value * $x + $coefficient;
        }
        $value
    }

    sub legendre-values(Int:D $degree, $x) {
        return (1e0, 0e0) if $degree == 0;
        my $previous = 1e0;
        my $current = $x;
        for 2 .. $degree -> $index {
            my $next = ((2 * $index - 1) * $x * $current - ($index - 1) * $previous) / $index;
            $previous = $current;
            $current = $next;
        }
        ($current, $previous)
    }

    # Solve a small dense system using partial pivoting.
    # This is used both for the Stieltjes polynomial and for its interpolatory weights.
    sub solve(@matrix, @right-hand-side) {
        my @a = @matrix.map(*.Array).Array;
        my @b = @right-hand-side.Array;
        my $size = @b.elems;

        for ^$size -> $column {
            my $pivot = $column;
            for $column + 1 ..^ $size -> $row {
                $pivot = $row if @a[$row][$column].abs > @a[$pivot][$column].abs;
            }
            die 'Unable to construct a Gauss-Kronrod rule.'
            if @a[$pivot][$column].abs < 1e-28;

            (@a[$column], @a[$pivot]) = (@a[$pivot], @a[$column]) if $pivot != $column;
            (@b[$column], @b[$pivot]) = (@b[$pivot], @b[$column]) if $pivot != $column;

            for $column + 1 ..^ $size -> $row {
                my $factor = @a[$row][$column] / @a[$column][$column];
                next if $factor == 0;
                for $column ..^ $size -> $index {
                    @a[$row][$index] -= $factor * @a[$column][$index];
                }
                @b[$row] -= $factor * @b[$column];
            }
        }

        my @answer = 0e0 xx $size;
        for reverse ^$size -> $row {
            my $sum = @b[$row];
            for $row + 1 ..^ $size -> $column {
                $sum -= @a[$row][$column] * @answer[$column];
            }
            @answer[$row] = $sum / @a[$row][$row];
        }
        @answer
    }
}