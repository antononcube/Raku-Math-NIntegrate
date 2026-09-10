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