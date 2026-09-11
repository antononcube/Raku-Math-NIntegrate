use v6.d;

use Math::NIntegrate::VariableTransformer;

class Math::NIntegrate::VariableTransformer::Affine
        is Math::NIntegrate::VariableTransformer {

    method transform(:@points, :$jacobian, Bool:D :fb(:$functional-bounds) = False --> Map:D) {
        my %bounds = self.get-transform-bounds();

        if $functional-bounds {
            die 'Functional boundaries variable transformation is not implemented yet.'
        } else {
            my %res;
            my %calc;
            if self.get-region.dimension == 1 {
                %calc = self!length-calc(%bounds<min>, %bounds<max>);
                %res<points> = %calc<min> <<+>> (@points >>*>> %calc<length>);
            } else {
                %calc = self!length-calc-md(%bounds<min>, %bounds<max>);
                %res<points> .= map({ %calc<min> <<+>> ($_ <<*>> %calc<length>) });
            }
            %res<jacobian> = $jacobian * %calc<jacobian>;
            return %res
        }
    }
}
