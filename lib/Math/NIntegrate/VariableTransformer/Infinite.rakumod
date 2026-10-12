use v6.d;

use Math::NIntegrate::VariableTransformer;
use Math::NIntegrate::Utilities;
use Math::NIntegrate::Codes;

class Math::NIntegrate::VariableTransformer::Infinite
        is Math::NIntegrate::VariableTransformer {

    method clone(-->Math::NIntegrate::VariableTransformer::Infinite) {
        Math::NIntegrate::VariableTransformer::Infinite.new(region => self.region).copy(self, :clone)
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
                    self.transforms[$i] = -> $point, $min, $max { self.fin-inf-transform($point, $min, $max) };

                    # 1 / x^2
                    self.jacobians[$i] = { 1 / $_ ** 2 }

                    self.min-transform-bounds[$i] = 0;
                    self.max-transform-bounds[$i] = 1;
                    self.max-original-bounds[$i] = %parsed<max-inf-dir>;

                    self.jacobian-factors[$i] = %parsed<max-inf-dir>
                }

                when VT_INF_FIN {
                    # max - (1 - x) / x
                    #self.transforms[$i] = { self.max-original-bounds[$i] - (1 - $_) / $_ }
                    self.transforms[$i] = -> $point, $min, $max { self.fin-inf-transform($point, $min, $max) };

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
                    self.transforms[$i] = -> $point, $min, $max { self.neg-inf-inf-transform($point, $min, $max) };

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

    #------------------------------------------------------
    # [finite-min, Inf)
    #------------------------------------------------------
    multi method fin-inf-transform(Numeric:D $point, Numeric:D $min, Numeric:D $max -->Map) {

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

    # En bloc transformation
    multi method fin-inf-transform(@points where @points.all ~~ Numeric:D, Numeric:D $min, Numeric:D $max -->Map) {

        my @zeroes = @points.map({ is-zero($_) ?? 0 !! 1 });
        my @pInvs = @points.kv.map( -> $i, $v { @zeroes[$i] ?? $max.sign * Inf !! numerical(1 / $v, self.working-precision) });
        my @tPoints = $min <<+>> ($max <<*>> (@pInvs <<->> 1));
        @pInvs = @pInvs.kv.map( -> $i, $v { @zeroes[$i] ?? 1 !! $v });
        my @tJacobians = @pInvs <<*>> @pInvs;

        return %(points => @tPoints, jacobian => [*] |@tJacobians)
    }

    #------------------------------------------------------
    # (-Inf, finite-max]
    #------------------------------------------------------
    multi method neg-inf-fin-transform(Numeric:D $point, Numeric:D $min, Numeric:D $max -->Map) {
        my %res = self.fin-inf-transform($point, -1 * $max, $min);
        %res<point> = -1 * %res<point>;
        return %res
    }

    # En bloc transformation
    multi method neg-inf-fin-transform(@points where @points.all ~~ Numeric:D, Numeric:D $min, Numeric:D $max -->Map) {
        my %res = self.fin-inf-transform(@points, -1 * $max, $min);
        %res<points> = -1 <<*>> %res<points>;
        return %res
    }

    #------------------------------------------------------
    # (-Inf, Inf)
    #------------------------------------------------------
    multi method neg-inf-inf-transform(Numeric:D $point, Numeric:D $min, Numeric:D $max -->Map) {
        my $tPoint;
        my $tJacobian;
        my $pInv;
        my $pInv1;
        if is-zero($point) {
            $tPoint = $min.sign * Inf;
            $tJacobian = 0;
        } elsif is-zero(1 - $point) {
            $tPoint = $max.sign * Inf;
            $tJacobian = 0;
        } else {
            $pInv = numerical(1 / $point, self.working-precision);
            $pInv1 = numerical(1 / (1 - $point), self.working-precision);
            $tPoint = $pInv1 - $pInv;
            $tJacobian = $pInv * $pInv + $pInv1 * $pInv1;
        }
        return %(point => $tPoint, jacobian => $tJacobian)
    }

    # En bloc transformation
    multi method neg-inf-inf-transform(@points where @points.all ~~ Numeric:D , Numeric:D $min, Numeric:D $max -->Map) {

        my @zeroes = @points.map({ is-zero($_) ?? 0 !! 1 });
        my @pInvs = @points.map(-> $i, $v { @zeroes[$i] ?? $min.sign * Inf !! numerical(1 / $v, self.working-precision) });
        my @pInvs1 = @points.map(-> $i, $v { @zeroes[$i] ?? $max.sign * Inf !! numerical(1 / (1 - $v), self.working-precision) });
        my @tPoints = @pInvs1 <<->> @pInvs;

        @pInvs = @pInvs.kv.map( -> $i, $v { @zeroes[$i] ?? 1 !! $v });
        @pInvs1 = @pInvs.kv.map( -> $i, $v { @zeroes[$i] ?? 1 !! $v });
        my @tJacobians = (@pInvs <<*>> @pInvs) <<+>> (@pInvs1 <<*>> @pInvs1);

        return %(points => @tPoints, jacobian => [*] |@tJacobians)
    }

    #------------------------------------------------------
    # [finite-min, finite-max]
    #------------------------------------------------------
    # No need

    #------------------------------------------------------
    # Transform point-wise
    #------------------------------------------------------
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

    method transform-en-bloc(
            :@points where @points.all ~~ (Array:D | List:D | Seq:D),
            :$jacobian is copy,
            Bool:D :fb(:$functional-bounds) = False,
            :$context = Nil,
            :$axis= Nil
            --> Map:D) {

        # Not using context
        # my %orig-bounds = self.get-original-bounds(:context);
        my %orig-bounds = self.get-original-bounds();

        # Should all points be checked?
        fail 'DIMENSIONS_DO_NOT_MATCH: for point argument and boundary' if @points.head.elems != %orig-bounds<min>.elems;

        if $functional-bounds || (%orig-bounds.values.flat(:hammer).any ~~ Callable:D) {
            die 'INCORRECT_ARGUMENTS: Affine en bloc works only with constant ranges.'
        } else {
            # Infinite range transformation for each axes separately
            my %res = points => ([] xx %orig-bounds<min>.elems).Array, :$jacobian;
            for ^self.min-original-bounds.elems -> $i {
                with self.transforms[$i] {
                    my @pointsLocal = @points.map(*[$i]);
                    my %h = self.transforms[$i](@pointsLocal, %orig-bounds<min>[$i], %orig-bounds<max>[$i]);
                    %res<points>[$i] = %h<points>;
                    %res<jacobian> = %res<jacobian> * %h<jacobian>;
                }
            }

            # Transpose the points
            %res<points> = |zip(|%res<points>);

            return %res
        }
    }
}
