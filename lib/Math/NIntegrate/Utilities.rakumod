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

proto sub numerical(Numeric:D $n, |) is export {*}

multi sub numerical(Numeric:D $n, Num) { $n.Num }

multi sub numerical(Numeric:D $n, $working-precision = Num) {
    numerical($n, :$working-precision)
}

multi sub numerical(Numeric:D $n, :wprec(:$working-precision) = Num) {
    return do given $working-precision {
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
# Frequent values
#==========================================================
# This might be a good idea -- numerical(1, self.working-precision) is used a lot.
# But I do not see much difference when profiling with Raku++.
#`[
proto sub ONE($prec = Num) is export {*}

multi sub ONE(Num) { 1e0 }
multi sub ONE(Rat) { 1.Rat }
multi sub ONE(FatRat) { 1.FatRat }
multi sub ONE(Numeric) { 1 }
]

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

        if $orig-min ~~ Inf | -Inf {
            $min = numerical($orig-min, $working-precision);
            $min-inf-dir = numerical($orig-min.sign, $working-precision)
        } else {
            $min = $orig-min ~~ Callable:D ?? $orig-min !! numerical($orig-min, $working-precision);
        }

        if $orig-max ~~ Inf | -Inf {
            $max = numerical($orig-max, $working-precision);
            $max-inf-dir = numerical($orig-max.sign, $working-precision)
        } else {
            $max = $orig-max ~~ Callable:D ?? $orig-max !! numerical($orig-max, $working-precision);
        }

        # Zero intervals should be handled. E.g. ('x', Inf, Inf)

        # Complex number cases
        if $real-flag && ($orig-min ~~ Complex:D || $orig-max ~~ Complex:D) {
            $real-flag = False
        }

        $bounds-case = do given (($orig-min ~~ Inf | -Inf), ($orig-max ~~ Inf | -Inf) ) {
                when $_.head && $_.tail { VT_INF_INF }
                when !$_.head && $_.tail { VT_FIN_INF }
                when $_.head && !$_.tail { VT_INF_FIN }
                when !$_.head && !$_.tail { VT_FIN_FIN }
        }

        return  %(:$min, :$max, :$min-inf-dir, :$max-inf-dir, :$real-flag, :$bounds-case)
}

#==========================================================
# Ranges or dimension spec
#==========================================================

#| Determine dimension
our sub determine-dimension($ranges, $dimension) {
    return do given ($ranges, $dimension) {
        when $_.head ~~ (Array:D | List:D | Seq:D) && $_.tail.isa(Whatever) {
            $ranges.elems
        }
        when $_.head.defined && $_.tail ~~ Int:D && $_.tail > 0 {
            warn 'The specified $ranges values is ignored, since $dimension is a positive integer';
            $_.tail
        }
        when $_.tail ~~ Int:D && $_.tail > 0 {
            # No change
            $_.tail
        }
        when $_.tail.defined {
            die 'The argument $dimensions is expected to be a positive integer,';
        }
        default {
            die 'At least one of the arguments $ranges and $dimension have to be specified.'
        }
    }
}

#======================================================
# Jacobian calculation related
#======================================================

# This method probably should be only called by length-calc-md --
# the framework should work only with points that are lists.
# (After initial processing.)

# There is no reason these subs to be a class methods in VariableTransformer,
# although that is tempting from a certain encapsulation perspective.
# Here memoization can be (more easily) applied.

#| Calculation of the interval length with working precision
our sub length-calc(Numeric:D $a is copy, Numeric:D $b is copy,
                    Bool:D :$mid-point = False,
                    :$working-precision = Num
        -->Map:D) {

    # Numeric evaluation of $a
    $a = numerical($a, $working-precision);

    # Numeric evaluation of $b
    $b = numerical($b, $working-precision);

    # The length of the interval
    my $length = $b - $a;

    # Mid-point calculation if needed
    my $middle = do if $mid-point {
        numerical(($a + $b) / 2, :$working-precision)
    } else {
        Whatever
    }

    return %(:$length, min => $a, :$middle, jacobian => $length)
}

#| Calculation of the multidimensional interval lengths with working precision
our sub length-calc-md(@a, @b,
                       Bool:D :$mid-point = False,
                       :$working-precision = Num
        -->Map:D) {
    die 'The sizes of the first two arguments are expected to match' unless @a.elems == @b.elems;

    my $jacobian = numerical(1, $working-precision);
    my @length;
    my @middle;
    my @min;
    for ^@a.elems -> $i {
        my %res = length-calc(@a[$i], @b[$i], :$mid-point, :$working-precision);
        @length.push(%res<length>);
        @middle.push(%res<middle>);
        @min.push(%res<min>);
        $jacobian *= %res<length>
    }

    return %(:@length, :@min, :@middle, :$jacobian)
}