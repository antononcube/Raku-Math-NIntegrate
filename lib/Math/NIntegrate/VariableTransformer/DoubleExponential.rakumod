use v6.d;

use Math::NIntegrate::VariableTransformer;
use Math::NIntegrate::Utilities;
use Math::NIntegrate::Codes;

class Math::NIntegrate::VariableTransformer::DoubleExponential
        is Math::NIntegrate::VariableTransformer {

    method clone(-->Math::NIntegrate::VariableTransformer::DoubleExponential) {
        Math::NIntegrate::VariableTransformer::DoubleExponential.new(region => self.region).copy(self, :clone)
    }

    # Very similar to TWEAK of Infinite
    submethod TWEAK(*%args) {

        # For each dimension
        for ^self.min-original-bounds.elems -> $i {
            # Parse boundaries for cases --
            # Math::NIntegrate::Codes::RangeBoundsCases

            my %parsed = parse-range-boundaries(self.min-original-bounds[$i], self.max-original-bounds[$i], working-precision => self.working-precision);

            self.min-original-bounds[$i] = %parsed<min>;
            self.max-original-bounds[$i] = %parsed<max>;

            given %parsed<bounds-case> {
                when VT_FIN_INF {
                    # min + exp( pi / 2 * sinh(x))
                    self.transforms[$i] = -> $point, $min, $max { self!fin-inf-transform($point, $min, $max) };

                    # exp( pi/2 * sinh(x) ) * pi/2 * cosh(x)
                    self.jacobians[$i] = { exp( pi/2 * sinh($_) ) * pi/2 * cosh($_) }

                    self.min-transform-bounds[$i] = -Inf;
                    self.max-transform-bounds[$i] = Inf;
                    self.max-original-bounds[$i] = %parsed<max-inf-dir>;

                    self.jacobian-factors[$i] = %parsed<max-inf-dir>
                }

                when VT_INF_FIN {
                    # max + exp( pi / 2 * sinh(-x))
                    self.transforms[$i] = -> $point, $min, $max { self!fin-inf-transform($point, $min, $max) };

                    # - exp( - pi/2 * sinh(x) ) * pi/2 * cosh(-x)
                    self.jacobians[$i] = { - exp( -pi/2 * sinh($_) ) * pi/2 * cosh($_) }

                    self.min-transform-bounds[$i] = -Inf;
                    self.max-transform-bounds[$i] = Inf;
                    self.max-original-bounds[$i] = %parsed<min-inf-dir>;
                    self.min-original-bounds[$i] = %parsed<max>;

                    self.jacobian-factors[$i] = -1 * %parsed<min-inf-dir>
                }

                when VT_INF_INF {
                    # sinh( (pi * sinh(x) ) / 2 )
                    self.transforms[$i] = -> $point, $min, $max { self!neg-inf-inf-transform($point, $min, $max) };

                    # (pi * cosh(x) * cosh( (pi * sinh(x) ) / 2) ) / 2
                    self.jacobians[$i] = { (pi * cosh($_) * cosh( (pi * sinh($_) ) / 2) ) / 2 }
                    self.min-transform-bounds[$i] = 0;
                    self.max-transform-bounds[$i] = 1;

                    self.jacobian-factors[$i] = %parsed<max-inf-dir>
                }

                when VT_FIN_FIN {
                    # <1/2> * ( tanh(<1/2> * π * sinh(t)) + 1 )
                    self.transforms[$i] = -> $point, $min, $max { self!fin-fin-transform($point, $min, $max) };

                    self.jacobians[$i] = { (pi * cosh($_) * sech((pi * sinh($_) ) / 2) ** 2 ) / 4 };
                    self.min-transform-bounds[$i] = %parsed<min>;
                    self.max-transform-bounds[$i] = %parsed<max>;

                    self.jacobian-factors[$i] = 1
                }

                default {
                    fail 'WRONG_TYPE: unknown variable transformer case'
                }
            }
        }
    }

    method !fin-inf-transform(Numeric:D $point, Numeric:D $min, Numeric:D $max -->Map) {

        # This might not work with the regular pi(π) for Rat and FatRat
        my $pi2 = numerical(pi / 2, self.working-precision);

        my $s = numerical(sinh($point), self.working-precision);
        my $c = numerical(cosh($point), self.working-precision);
        my $z = $pi2 * $s;
        my $v = exp($z);

        my $exp = sign($max) * $v;
        my $tPoint = $min + $exp;
        my $tJacobian = $pi2 * $v * $c;

        return %(point => $tPoint, jacobian => $tJacobian)
    }

    method !neg-inf-fin-transform(Numeric:D $point, Numeric:D $min, Numeric:D $max -->Map) {
        my %res = self!fin-inf-transform($point, -1 * $max, $min);
        %res<point> = -1 * %res<point>;
        return %res
    }

    method !neg-inf-inf-transform(Numeric:D $point, Numeric:D $min, Numeric:D $max -->Map) {
        # This might not work with the regular pi(π) for Rat and FatRat
        my $pi2 = numerical(pi / 2, self.working-precision);

        my $s = numerical(sinh($point), self.working-precision);
        my $c = numerical(cosh($point), self.working-precision);
        my $z = sinh( $pi2 * $s );

        my $tPoint = $z;

        my $tJacobian = $pi2 * cosh( $pi2 * $s ) * $c;

        return %(point => $tPoint, jacobian => $tJacobian)
    }

    method !fin-fin-transform(Numeric:D $point, Numeric:D $min, Numeric:D $max -->Map) {

        # This might not work with the regular pi(π) for Rat and FatRat
        my $pi2 = numerical(pi / 2, self.working-precision);

        my ($u, $v, $w);
        my $pInv;
        my $tPoint;
        my $tJacobian;
        my $exp;

        $u = exp($point);
        $v = 1 / $u;
        $w = $u + $v;
        $u = exp( $u * $pi2 );
        $v = exp( $v * $pi2 );

        $pInv = 1 / ($u + $v);
        $tJacobian = $w * $u * $v * $pInv * $pInv;

        $exp = $point < 0 ?? ($max - $min) * $u * $pInv !! ($min - $max) * $v * $pInv;

        $tPoint = $point < 0 ?? $min + $exp !! $max + $exp;

        return %(point => $tPoint, jacobian => $tJacobian)
    }

    method transform(
            :@point is copy,
            :$jacobian is copy,
            Bool:D :fb(:$functional-bounds) = False,
            :$context = Nil,
            :$axis = Nil
            --> Map:D) {

        # From context or its own?
        #my %orig-bounds = self.get-original-bounds(:$context);
        my %orig-bounds = self.get-original-bounds();

        if $functional-bounds || (%orig-bounds.values.flat(:hammer).any ~~ Callable:D) {
            die 'For functional boundaries computations $axis is expected to be a non-negative integer within the integral dimensions.'
            unless $axis ~~ Int:D && 0 ≤ $axis < %orig-bounds<min>.elems;

            # Infinite transformation for a single axis
            my %res = :@point, :$jacobian;
            with self.transforms[$axis] {

                my $min = self.min-original-bounds[$axis];
                $min = $min(|@point.head($axis)) if $min ~~ Callable:D;

                my $max = self.max-original-bounds[$axis];
                $max = $max(|@point.head($axis)) if $max ~~ Callable:D;

                my %h = self.transforms[$axis](@point[$axis], $min, $max);
                %res<point>[$axis] = %h<point>;
                %res<jacobian> = %res<jacobian> * %h<jacobian>;
            }
            return %res

        } else {
            my %res = :@point, :$jacobian;
            for ^self.min-original-bounds.elems -> $i {
                with self.transforms[$i] {
                    my %h = self.transforms[$i](@point[$i], %orig-bounds<min>[$i], %orig-bounds<max>[$i]);
                    %res<point>[$i] = %h<point>;
                    %res<jacobian> = %res<jacobian> * %h<jacobian>;
                }
            }
            return %res
        }
    }

}