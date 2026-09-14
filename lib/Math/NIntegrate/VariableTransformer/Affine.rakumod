use v6.d;

use Math::NIntegrate::VariableTransformer;

class Math::NIntegrate::VariableTransformer::Affine
        is Math::NIntegrate::VariableTransformer {

    method transform(:@point, :$jacobian, Bool:D :fb(:$functional-bounds) = False --> Map:D) {
        my %bounds = self.get-transform-bounds();

        if $functional-bounds {
            die 'Functional boundaries variable transformation is not implemented yet.'
        } else {
            my %res;
            my %calc = self.length-calc-md(%bounds<min>, %bounds<max>);
            %res<point> = @point.kv.map( -> $i, $p { %calc<min>[$i] + ($p * %calc<length>[$i]) });
            %res<jacobian> = $jacobian * %calc<jacobian>;
            return %res
        }
    }
}
