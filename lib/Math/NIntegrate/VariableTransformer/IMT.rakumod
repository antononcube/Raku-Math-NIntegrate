use v6.d;

use Math::NIntegrate::VariableTransformer;
use Math::NIntegrate::Utilities;
use Math::NIntegrate::Codes;

#| Iri, Moriguti and Takesawa transformation
class Math::NIntegrate::VariableTransformer::IMT
        is Math::NIntegrate::VariableTransformer {

    method clone(-->Math::NIntegrate::VariableTransformer::IMT) {
        Math::NIntegrate::VariableTransformer::IMT.new(region => self.region).copy(self, :clone)
    }

    method set-transform-axis(UInt:D $i) {
        self.transforms[$i] = -> $point, $min, $max { self!singfin-fin-transform($point, $min, $max) }
    }

    method drop-transform-axis(UInt:D $i) {
        self.transforms[$i] = WhateverCode
    }

    method !singfin-fin-transform(Numeric:D $point, Numeric:D $min, Numeric:D $max -->Map) {

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
            my $exp = exp(1 - $pInv);
            $tPoint = $min + ($max - $min) * $exp;
            $tJacobian = ($max - $min) * $exp * $pInv * $pInv;
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
