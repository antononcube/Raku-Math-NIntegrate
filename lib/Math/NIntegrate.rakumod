use v6.d;

unit module Math::NIntegrate;

use Math::NIntegrate::Builder;
use Math::NIntegrate::Spec;

#| Numerical integration of a function for given ranges using method- and precision specifications.
#| C<&f> -- integrand.
#| C<*@ranges> -- ranges like ('x', 0, 10) or ('$x', 1, 3), ('$y', 0, -> $x { $x } )
#| C<*%args> -- options like max-recursion => 12, etc.
proto sub nintegrate(&f, *@ranges,
                     :$method is copy = Whatever,
                     :prec(:$working-precision) = Whatever,
                     :acc(:$accuracy) = Whatever,
                     :p(:$pairs) = False,
                     *%args) is export {*}

multi sub nintegrate(&f, *@ranges,
                     :$method = Whatever,
                     :prec(:$working-precision) = Whatever,
                     :acc(:$accuracy) = Whatever,
                     :p(:$pairs) = False,
                     *%args) {

    # Check and normalize
    my $spec = Math::NIntegrate::Spec.new(&f, @ranges, %(:$method, :$working-precision, :$accuracy, |%args));

    # Builder object
    my $builder = Math::NIntegrate::Builder.new;

    # Integrator object
    my $integrator = $builder.make-integrator($spec.numerical-function, ranges => $spec.ranges, options => $spec.options);

    # Integrate
    my %res = $integrator.algorithm(
            relative-tolerance => $spec.options<relative-tolerance>,
            absolute-tolerance => $spec.options<absolute-tolerance>,
            );

    return $pairs ?? %res !! %res<integral>
}