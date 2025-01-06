use v6.d;

unit module Math::NIntegrate;

#| Numerical integration of a function for given ranges using method- and precision specifications.
#| C<&f> -- integrand.
#| C<*%args> -- should have ranges in the form Str:D => (Numeric:D, Numeric:D).
proto sub nintegrate(&f, *%args, :$method = Whatever, :prec(:$working-precision) = Whatever) is export {*}