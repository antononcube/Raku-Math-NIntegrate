use v6.d;

use Math::NIntegrate::VariableTransformer;

class Math::NIntegrate::VariableTransformer::Affine
        is Math::NIntegrate::VariableTransformer {

    method transform(
            :@point is copy,
            :$jacobian is copy,
            Bool:D :fb(:$functional-bounds) = False,
            :$context = Nil
            --> Map:D) {
        # If the Affine object has a context (most likely, a Composite object)
        # then the transform boundaries of the context object are used.
        my %bounds = self.get-transform-bounds(:$context);

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
