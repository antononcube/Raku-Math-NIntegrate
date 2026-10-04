use v6.d;

use Math::NIntegrate::VariableTransformer;
use Math::NIntegrate::Utilities;
use Math::NIntegrate::Codes;

#| Iri, Moriguti and Takesawa transformation
class Math::NIntegrate::VariableTransformer::IMT
        is Math::NIntegrate::VariableTransformer {

    has @.derivatives;

    method clone(-->Math::NIntegrate::VariableTransformer::IMT) {
        Math::NIntegrate::VariableTransformer::IMT.new(region => self.region).copy(self, :clone)
    }

    method set-transform-axis(UInt:D $i) {
        self.transforms[$i] = -> $point, $min, $max { self!singfin-fin-transform($point, $min, $max) }
    }

    method drop-transform-axis(UInt:D $i) {
        self.transforms[$i] = WhateverCode
    }

    method !singfin-fin-transform(
            Numeric:D $point,
            Numeric:D $min, Numeric:D $max,
            Numeric:D :$a = 10, Numeric:D :$p = 1
            -->Map) {

        my $tPoint;
        my $tJacobian;
        my $pInv;
        if is-zero($point) {
            $tPoint = 0;
            $tJacobian = 0;
        } else {
            # This for p==1
            $pInv = numerical(1 / $point, self.working-precision);
            # This is not high precision. At some point continued fractions can be used.
            my $exp = exp($a * (1 - $pInv));
            $tPoint = $min + ($max - $min) * $exp;
            $tJacobian = ($max - $min) * $exp * $pInv * $pInv * $a;
        }
        return %(point => $tPoint, jacobian => $tJacobian)
    }

    method !fin-singinf-transform(Numeric:D $point, Numeric:D $min, Numeric:D $max -->Map) {!!!}

    method !singfin-inf-transform(Numeric:D $point, Numeric:D $min, Numeric:D $max -->Map) {!!!}

    # Similar as Infinity.
    # When IMT is added to the Composite stack, the axis for which IMT is applied is specified.
    # Additional IMT applications (for other axes) use that stack object.
    method transform(
            :@point is copy,
            :$jacobian is copy,
            Bool:D :fb(:$functional-bounds) = False,
            :$context = Nil,
            :$axis = Nil
            --> Map:D) {
        my %bounds = self.get-transform-bounds(:$context);

        # Derivative signs are used to known in which
        # direction to make a step back from singular points.
        if @!derivatives.elems == 0 {
            @!derivatives = 0 xx %bounds<min>.elems
        }

        if $functional-bounds || (%bounds.values.flat(:hammer).any ~~ Callable:D) {
            die 'For functional boundaries computations $axis is expected to be a non-negative integer within the integral dimensions.'
            unless 0 ≤ $axis < %bounds<min>.elems;

            # IMT transformation for a single axis
            my %res = :@point, :$jacobian;
            with self.transforms[$axis] {
                say "Infinity for $axis";
                my %h = self.transforms[$axis](@point[$axis], self.min-original-bounds[$axis], self.max-original-bounds[$axis]);
                %res<point>[$axis] = %h<point>;
                %res<jacobian> = %res<jacobian> * %h<jacobian>;
                @!derivatives[$axis] = %h<jacobian>;
            }
            return %res

        } else {
            my %res = :@point, :$jacobian;
            for ^self.min-original-bounds.elems -> $i {
                with self.transforms[$i] {
                    my %h = self.transforms[$i](@point[$i], self.min-original-bounds[$i], self.max-original-bounds[$i]);
                    %res<point>[$i] = %h<point>;
                    %res<jacobian> = %res<jacobian> * %h<jacobian>;
                    @!derivatives[$i] = %h<jacobian>;
                }
            }
            return %res
        }
    }
}
