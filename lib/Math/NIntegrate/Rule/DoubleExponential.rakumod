use v6.d;

use Math::NIntegrate::Utilities;
use Math::NIntegrate::Rule::General;

# Should it inherit ::Rule::General or ::Rule ?

class Math::NIntegrate::Rule::DoubleExponential
        is Math::NIntegrate::Rule::General {

    # Should the class have a step attribute?
    # Or "just" a step argument for method integrate?
    # What is a good default value of the integration step? (E.g. 1.)
    # Should take working precision or tolerance as an argument --
    # in order to figure out when to stop accumulating the terms sum.

    sub fixed-point(&F, $start, $h, $increment) {
        my $i = $start;
        my $sum = 0;

        loop {
            my $next = $sum + &F($i * $h);
            $i += $increment;
            last if $next == $sum;
            $sum = $next;
        }

        $sum
    }


    method integrate($region, Numeric:D $step = 1) {

        # Computing terms of the infinite range trapezoidal formula in both directions:
        # LaTeX: \int_{-\infty }^{+\infty } f(\phi (t)) \phi '(t) \, dx
        # ASCII: ∫ [-∞,+∞] f(φ(t)) φ'(t) dx
        # where phi is one of the transformations in ::VariableTransform::DoubleExponential
    }

    method integrate-dim1($region, Numeric:D $step = 1) {!!!}

    # Essentially, dynamic Cartesian rule application
    method integrate-dimN($region, Numeric:D $step = 1) {!!!}
}
