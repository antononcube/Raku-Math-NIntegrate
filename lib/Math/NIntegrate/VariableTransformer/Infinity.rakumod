use v6.d;

use Math::NIntegrate::VariableTransformer;
use Math::NIntegrate::Utilities;
use Math::NIntegrate::Codes;

class Math::NIntegrate::VariableTransformer::Infinity
        is Math::NIntegrate::VariableTransformer {

    method clone(-->Math::NIntegrate::VariableTransformer::Infinity) {
        Math::NIntegrate::VariableTransformer::Infinity.new(region => self.region).copy(self, :clone)
    }

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
                    # min + (1 - x) / x
                    #self.transforms[$i] = { self.min-original-bounds[$i] + (1 - $_) / $_ }
                    self.transforms[$i] = -> $point, $min, $max { self!fin-inf-transform($point, $min, $max) };

                    # 1 / x^2
                    self.jacobians[$i] = { 1 / $_ ** 2 }

                    self.min-transform-bounds[$i] = 0;
                    self.max-transform-bounds[$i] = 1;
                    self.max-original-bounds[$i] = %parsed<max-inf-dir>;

                    self.jacobian-factors[$i] = %parsed<max-inf-dir>
                }

                when VT_INF_FIN {
                    # min - (1 - x) / x
                    #self.transforms[$i] = { self.max-original-bounds[$i] - (1 - $_) / $_ }
                    self.transforms[$i] = -> $point, $min, $max { self!fin-inf-transform($point, $min, $max) };

                    # 1 / x^2
                    self.jacobians[$i] = { 1 / $_ ** 2 }

                    self.min-transform-bounds[$i] = 0;
                    self.max-transform-bounds[$i] = 1;
                    self.max-original-bounds[$i] = %parsed<min-inf-dir>;
                    self.min-original-bounds[$i] = %parsed<max>;

                    self.jacobian-factors[$i] = -1 * %parsed<min-inf-dir>
                }

                when VT_INF_INF {
                    # (1 - x) / x - 1 / x
                    #self.transforms[$i] = { self.max-original-bounds[$i] - (1 - $_) / $_ }
                    self.transforms[$i] = -> $point, $min, $max { self!neg-inf-inf-transform($point, $min, $max) };

                    # x^-2 + (1-x)^-2
                    self.jacobians[$i] = { 1 / $_ ** 2 + 1 / (1 - $_) ** 2 }
                    self.min-transform-bounds[$i] = 0;
                    self.max-transform-bounds[$i] = 1;

                    self.jacobian-factors[$i] = 1
                }

                when VT_FIN_FIN {
                    self.transforms[$i] = WhateverCode;
                    self.jacobians[$i] = Whatever;
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

        my $tPoint;
        my $tJacobian;
        my $pInv;
        if is-zero($point) {
            $tPoint = $max.sign * Inf;
            $tJacobian = 0;
        } else {
            $pInv = numerical(1 / $point, self.working-precision);
            $tPoint = $min + $max * ($pInv - 1);
            $tJacobian = $pInv * $pInv;
        }
        return %(point => $tPoint, jacobian => $tJacobian)
    }

    method !neg-inf-fin-transform(Numeric:D $point, Numeric:D $min, Numeric:D $max -->Map) {
        my %res = self!fin-inf-transform($point, -1 * $max, $min);
        %res<point> = -1 * %res<point>;
        return %res
    }

    method !neg-inf-inf-transform(Numeric:D $point, Numeric:D $min, Numeric:D $max -->Map) {
        return %(:$point, jacobian => 1)
    }

    method transform(
            :@point is copy,
            :$jacobian is copy,
            Bool:D :fb(:$functional-bounds) = False,
            :$context = Nil
            --> Map:D) {
        my %bounds = self.get-transform-bounds(:$context);

        if $functional-bounds {
            die 'Functional boundaries variable transformation is not implemented yet.'
        } else {
            my %res = :@point, :$jacobian;
            for ^self.min-original-bounds.elems -> $i {
                with self.transforms[$i] {
                    my %h = self.transforms[$i](@point[$i], self.min-original-bounds[$i], self.max-original-bounds[$i]);
                    %res<point>[$i] = %h<point>;
                    %res<jacobian> = %res<jacobian> * %h<jacobian>;
                }
            }
            return %res
        }
    }

}
