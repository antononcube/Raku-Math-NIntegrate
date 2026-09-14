use v6.d;

unit module Math::NIntegrate::Utilities;

use Math::NIntegrate::Codes;

#==========================================================
# Machine epsilon
#==========================================================

our $MACHINE_EPSILON is export = 2.220446049250313e-16;

#| Machine epsilon computation
our sub compute-machine-epsilon() {
    my $eps = 1.0e0;
    while 1.0e0 + $eps / 2.0e0 > 1.0e0 {
        $eps /= 2.0e0;
    }
    $MACHINE_EPSILON = $eps;
    return $MACHINE_EPSILON;
}


#==========================================================
# Effective zero
#==========================================================

#| Effective zero check
sub is-zero(Numeric:D $x, :tol(:$tolerance) = 2 * $MACHINE_EPSILON) is export {
    $x.abs ≤ $tolerance;
}

#==========================================================
# Numerical
#==========================================================

sub numerical(Numeric:D $n, $wprec = Num, $default = Num) is export {
    return do given $wprec {
        when FatRat { $n.FatRat }
        when Rat { $n.Rat }
        when Num { $n.Num }
        when Numeric { $n.Numeric }
        default {
            die 'Unknown working precision specification.'
        }
    }
}

#==========================================================
# Parse range boundaries
#==========================================================

#| Parse range boundaries.
sub parse-range-boundaries($orig-min,
                           $orig-max,
                           :$real-flag is copy = True,
                           :$working-precision = Num
        --> Map:D) is export {

        my $min;
        my $max;
        my $min-inf-dir = 0;
        my $max-inf-dir = 0;
        my $bounds-case;

        if $orig-min.isa(Inf) {
                $min = $orig-min;
                $min-inf-dir = $orig-min.sign
        } else {
                $min = numerical($orig-min, $working-precision);
        }

        if $orig-max.isa(Inf) {
                $max = $orig-max;
                $max-inf-dir = $orig-max.sign
        } else {
                $max = numerical($orig-max,  $working-precision);
        }

        # Zero intervals should be handled. E.g. ('x', Inf, Inf)

        # Complex number cases
        if $real-flag && ($orig-min ~~ Complex:D || $orig-max ~~ Complex:D) {
                $real-flag = False
        }

        $bounds-case = do given ($orig-min.isa(Inf), $orig-max.isa(Inf)) {
                when $_.head && $_.tail { VT_INF_INF }
                when !$_.head && $_.tail { VT_FIN_INF }
                when $_.head && !$_.tail { VT_INF_FIN }
                when !$_.head && !$_.tail { VT_FIN_FIN }
        }

        return  %(:$min, :$max, :$min-inf-dir, :$max-inf-dir, :$real-flag, :$bounds-case)
}