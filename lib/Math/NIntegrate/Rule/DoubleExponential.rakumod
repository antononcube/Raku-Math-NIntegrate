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

    sub fixed-point($region, $start, $h, $increment) {
        my $i = $start;
        my $sum = 0;

        loop {
            my $next = $sum + $region.eval-integrand([$i * $h, ]);
            $i += $increment;
            last if $next == $sum;
            $sum = $next;
        }

        return $sum
    }

    sub de-rule-estimate($region, $h, $oldSum = Nil) {
        my $step = $oldSum === Nil ?? 1 !! 2;

        my $s1 = fixed-point($region, 1, $h, $step);
        my $s2 = fixed-point($region, -1, $h, -$step);

        return $oldSum === Nil
                ?? $h * ($s1 + $region.eval-integrand([0,]) + $s2)
                !! $h * ($s1 + $s2) + $oldSum / 2
    }

    method integrate($region) {

        my ($sum, $h, $k);
        $sum = $region.reuse-values<sum> // Nil;
        $h = $region.reuse-values<h> // 1;
        $k = $region.reuse-values<k> // 0;

        my $previous = de-rule-estimate($region, $h, Nil);
        my $new = $previous;

        $h /= 2;
        $new = de-rule-estimate($region, $h, $previous);
        $k++;

        $region.reuse-values = %(:$sum, :$h, :$k);

        self.integral = $new;
        self.error = ($new - $previous).abs;
        self.largest-error-axis = 0;

        return self
    }

    #| Full, on the spot integration. (I.e. a "mini" strategy.)
    method integrate-full($region, Numeric:D $step = 1, Numeric:D :$tolerance = 10e-10, UInt:D :$max-recursion = 7) {

        # Computing terms of the infinite range trapezoidal formula in both directions:
        # LaTeX: \int_{-\infty }^{+\infty } f(\phi (t)) \phi '(t) \, dx
        # ASCII: ∫ [-∞,+∞] f(φ(t)) φ'(t) dx
        # where phi is one of the transformations in ::VariableTransform::DoubleExponential

        my $h = 1;
        my $previous = de-rule-estimate($region, $h, Nil);
        my $temp = $previous;
        my $k = 0;

        loop {
            $h /= 2;
            $temp = de-rule-estimate($region, $h, $previous);
            $k++;

            last if $k >= $max-recursion;
            last unless abs($previous) * $tolerance <= abs($previous - $temp);

            $previous = $temp;
        }

        $temp
    }

    method integrate-dim1($region, Numeric:D $step = 1) {!!!}

    # Essentially, dynamic Cartesian rule application
    method integrate-dimN($region, Numeric:D $step = 1) {!!!}
}
