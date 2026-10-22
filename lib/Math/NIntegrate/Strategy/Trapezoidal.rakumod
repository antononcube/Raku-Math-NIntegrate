use v6.d;

use Math::NIntegrate::Strategy::DoubleExponential;

# Trapezoidal turns into DoubleExponential for infinite ranges.
# In some sense, it is just choosing/manipulating the variable transformers.

#| Trapezoidal strategy class
class Math::NIntegrate::Strategy::Trapezoidal
        is Math::NIntegrate::Strategy::DoubleExponential {

}
