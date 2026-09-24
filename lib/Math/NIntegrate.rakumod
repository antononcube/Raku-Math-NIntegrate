use v6.d;

unit module Math::NIntegrate;

use Math::NIntegrate::Builder;
use Math::NIntegrate::Spec;

# Are these needed or good ideas?
#| Synonym of &nintegrate.
# sub numerical-integration(*@args, *%args) is export { nintegrate(|@args, |%args) }

#| Synonym of &nintegrate.
# sub NIntegrate(*@args, *%args) is export { nintegrate(|@args, |%args) }

#| Numerical integration of a function for given ranges using method- and precision specifications.
#| C<&f> -- integrand.
#| C<*@ranges> -- ranges like ('x', 0, 10) or ('$x', 1, 3), ('$y', 0, -> $x { $x } )
#| C<*%args> -- options like max-recursion => 12, etc.
proto sub nintegrate(+@args,
                     :$method is copy = Whatever,
                     :prec(:$working-precision) = Whatever,
                     :acc(:$accuracy) = Whatever,
                     :p(:$pairs) = False,
                     *%args) is export {*}

multi sub nintegrate(+@args,
                     :$method = Whatever,
                     :prec(:$working-precision) = Whatever,
                     :acc(:$accuracy) = Whatever,
                     :p(:$pairs) = False,
                     *%args) {
    die 'At least two positional arguments are expected.'
    unless @args.elems ≥ 2;

    die 'The first argument is expected to be a callable.'
    unless @args.head ~~ Callable:D;

    my &f = @args.head;

    my @ranges = @args.tail(*-1);

    # Check and normalize
    my $spec = Math::NIntegrate::Spec.new(&f, @ranges, %(:$method, :$working-precision, :$accuracy, |%args));

    # Builder object
    my $builder = Math::NIntegrate::Builder.new;

    # Integrator object
    my $integrator = $builder.make-integrator(
            integrand => $spec.integrand,
            ranges => $spec.ranges,
            |$spec.options);

    # Integrate
    my %res = $integrator.algorithm(
            relative-tolerance => $spec.options<relative-tolerance>,
            absolute-tolerance => $spec.options<absolute-tolerance>,
            );

    return $pairs ?? %res !! %res<integral>
}