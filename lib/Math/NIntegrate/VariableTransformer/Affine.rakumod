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
            :$context = Nil,
            :$axis= Nil
            --> Map:D) {
        # If the Affine object has a context (most likely, a Composite object)
        # then the transform boundaries of the context object are used.
        my %bounds = self.get-transform-bounds(:$context);

        fail 'DIMENSIONS_DO_NOT_MATCH: for point argument and boundary' if @point.elems != %bounds<min>.elems;

        if $functional-bounds || (%bounds.values.flat(:hammer).one ~~ Callable:D) {
            die 'For functional boundaries computations $axis is expected to be a non-negative integer within the integral dimensions.'
            unless 0 ≤ $axis < %bounds<min>.elems;

            # Affine transformation for a single axis, with Callable boundaries evaluation first
            my $min = %bounds<min>[$axis] ~~ Callable ?? %bounds<min>[$axis](|@point.head($axis)) !! %bounds<min>[$axis];
            my $max = %bounds<max>[$axis] ~~ Callable ?? %bounds<max>[$axis](|@point.head($axis)) !! %bounds<max>[$axis];

            my %calc = self.length-calc($min, $max);
            my @res-point = @point;
            @res-point[$axis] = %calc<min> + @point[$axis] * %calc<length>;
            my $res-jacobian = $jacobian * %calc<jacobian>;
            return %(point => @res-point, jacobian => $res-jacobian)

        } else {
            # Affine transformation for all axes
            my %calc = self.length-calc-md(%bounds<min>, %bounds<max>);
            my @res-point = @point.kv.map(-> $i, $p {%calc<min>[$i] + ($p * %calc<length>[$i])});
            my $res-jacobian = $jacobian * %calc<jacobian>;
            return %(point => @res-point, jacobian => $res-jacobian)
        }
    }
}
