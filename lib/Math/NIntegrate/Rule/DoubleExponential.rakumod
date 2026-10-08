use v6.d;

use Math::NIntegrate::Utilities;
use Math::NIntegrate::Rule::General;

# Should it inherit ::Rule::General or ::Rule ?

class Math::NIntegrate::Rule::DoubleExponential
        is Math::NIntegrate::Rule::General {

    # Should the class have a step attribute?
    # Or "just" a step argument for method integrate?
    # What is a good default value of the integration step? (E.g. 1.)

    method integrate($region, Numeric:D $step = 1) {

        # Computing terms of the infinite range trapezoidal formula in both directions:
        # LaTeX: \int_{-\infty }^{+\infty } f(\phi (t)) \phi '(t) \, dx
        # ASCII: int_{-infty}^{+infty} f(phi(t)) phi'(t) dx
        # where phi is one of the transformations in ::VariableTransform::DoubleExponential
    }
}
