use v6.d;

use Math::NIntegrate::VariableTransformer::Affine;
use Math::NIntegrate::Utilities;

class Math::NIntegrate::VariableTransformer::AffineEnBloc
        is Math::NIntegrate::VariableTransformer::Affine {

    method clone(-->Math::NIntegrate::VariableTransformer::AffineEnBloc) {
        Math::NIntegrate::VariableTransformer::AffineEnBloc.new(region => self.region).copy(self, :clone)
    }

    method transform(
            :@points is copy where @points.all ~~ Positional:D,
            :$jacobian is copy,
            Bool:D :fb(:$functional-bounds) = False,
            :$context = Nil,
            :$axis= Nil
            --> Map:D) {
        # If the Affine object has a context (most likely, a Composite object)
        # then the transform boundaries of the context object are used.
        my %bounds = self.get-transform-bounds(:$context);

        # Should all points be checked?
        fail 'DIMENSIONS_DO_NOT_MATCH: for point argument and boundary' if @points.head.elems != %bounds<min>.elems;

        if $functional-bounds || (%bounds.values.flat(:hammer).any ~~ Callable:D) {
            die 'INCORRECT_ARGUMENTS: AffineEnBloc works only with constant ranges.'
        } else {
            # Affine transformation for all axes
            my %calc = Math::NIntegrate::Utilities::length-calc-md(%bounds<min>, %bounds<max>, :!mid-point, working-precision => self.working-precision);
            # This can be optimized using Math::SparseMatrix::Native
            my @res-points = @points.map({ %calc<min> <<+>> ($_ <<*>> %calc<length>) });
            my $res-jacobian = $jacobian * %calc<jacobian>;
            return %(points => @res-points, jacobian => $res-jacobian)
        }
    }
}
