use v6.d;

use Math::NIntegrate::VariableTransformer;

class Math::NIntegrate::VariableTransformer::Affine
        is Math::NIntegrate::VariableTransformer {

    method transform(:@point, :$jacobian, Bool:D :fb(:$functional-bounds) = False --> Map:D) {
        my %bounds = self.get-transform-bounds();

        fail 'DIMENSIONS_DO_NOT_MATCH: for point argument and boundary' if @point.elems != %bounds<min>.elems;

        if $functional-bounds {
            die 'Functional boundaries variable transformation is not implemented yet.'
        } else {
            my %calc = self.length-calc-md(%bounds<min>, %bounds<max>);
            my @res-point = @point.kv.map( -> $i, $p { %calc<min>[$i] + ($p * %calc<length>[$i]) });
            my $res-jacobian = $jacobian * %calc<jacobian>;
            return %(point => @res-point, jacobian => $res-jacobian)
        }
    }
}
