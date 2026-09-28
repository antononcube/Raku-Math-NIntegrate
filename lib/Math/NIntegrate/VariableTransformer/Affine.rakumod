use v6.d;

use Math::NIntegrate::VariableTransformer;

class Math::NIntegrate::VariableTransformer::Affine
        is Math::NIntegrate::VariableTransformer {

    method clone(-->Math::NIntegrate::VariableTransformer::Affine) {
        Math::NIntegrate::VariableTransformer::Affine.new(region => self.region).copy(self, :clone)
    }

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

        if $functional-bounds || (%bounds.values.flat(:hammer).one ~~ Callable:D) {
            die 'Functional boundaries variable transformation is not implemented yet.';
            # We can verify that:
            # - dimension is higher than 1,
            # - the first range is with numbers
            # - at least one of the boundaries has a Callable
            # We have to match the variables and/or the arity of the callables.
            # This probably means using NumericalFunction objects for callable boundaries.
            for (1 ..^ %bounds<min>.elems) -> $i {
                if %bounds<min>[$i] ~~ Callable:D {
                    %bounds<min>[$i] = %bounds<min>[$i].(|@point.head($i))
                }

                if %bounds<max>[$i] ~~ Callable:D {
                    %bounds<max>[$i] = %bounds<max>[$i].(|@point.head($i))
                }
            }
        } else {
            my %calc = self.length-calc-md(%bounds<min>, %bounds<max>);
            my @res-point = @point.kv.map(-> $i, $p {%calc<min>[$i] + ($p * %calc<length>[$i])});
            my $res-jacobian = $jacobian * %calc<jacobian>;
            return %(point => @res-point, jacobian => $res-jacobian)
        }
    }
}
