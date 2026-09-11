use v6.d;

unit module Math::NIntegrate::Utilities;

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
# Range boundaries cases
#==========================================================

#| Range boundaries cases
our enum RangeBoundsCases is export <VT_A_B VT_A_INF VT_INF_B VT_INF_INF>;

#==========================================================
# Parse range boundaries
#==========================================================

#| Parse range boundaries.
sub parse-range-boundaries($orig-min is copy,
                           $orig-max is copy,
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
                $orig-min = numerical($orig-min);
                $min = $orig-min
        }

        if $orig-max.isa(Inf) {
                $max = $orig-max;
                $max-inf-dir = $orig-max.sign
        } else {
                $orig-max = numerical($orig-max);
                $max = $orig-max
        }

        # Zero intervals should be handled. E.g. ('x', Inf, Inf)

        # Complex number cases
        if $real-flag && ($orig-min ~~ Complex:D || $orig-max ~~ Complex:D) {
                $real-flag = False
        }

        $bounds-case = do given ($orig-min.isa(Inf), $orig-max.isa(Inf)) {
                when (True, True) { VT_INF_INF }
                when (False, True) { VT_A_INF }
                when (True, False) { VT_INF_B }
                when (False, False) { VT_A_B }
        }

        return  %(:$min, :$max, :$min-inf-dir, :$max-inf-dir, :$real-flag, :$bounds-case)
}