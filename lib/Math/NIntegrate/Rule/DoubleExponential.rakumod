use v6.d;

use Math::NIntegrate::Utilities;
use Math::NIntegrate::Rule::General;

# Should it inherit ::Rule::General or ::Rule ?

class Math::NIntegrate::Rule::DoubleExponential
        is Math::NIntegrate::Rule::General {

    has @!hyper-plane;
    has $.evaluations-count = 0;

    method !eval-integrand($region, Numeric:D $h, Numeric:D $point) {
        return do if @!hyper-plane.elems + 1 == self.dimension {
            $!evaluations-count++;
            $region.eval-integrand(@!hyper-plane.clone.push($point))
        } else {
            @!hyper-plane .= push($point);
            my $res = self!de-rule-estimate($region, $h / self.dimension, Nil);
            @!hyper-plane .= head(*-1);
            $res
        }
    }

    method !fixed-point($region, $start, $h, $increment) {
        my $i = $start;
        my $sum = 0;

        loop {
            my $next = $sum + self!eval-integrand($region, $h, $i * $h);
            $i += $increment;

            # This checks might relax the stopping criteria:
            #   is-zero($point) || is-zero(1 - $point)
            # where $point is transformed by ::VariableTransformer::DoubleExponential from (-Inf, Inf) to (0, 1).
            # Speed-wise -- for 1D -- there might be not much time saving since there is another variable transformation.
            # After some experiments it did not perform that well, e.g., 1/sqrt(x) over [0, 1] takes 3 times more evaluations.
            # This checks might have a positive effect if higher precision calculations are implemented.
            #my $point = $i * $h;
            #my $axis = @!hyper-plane.elems;
            #$point = $region.variable-transformer.transforms[$axis]($point, 0, 1)<point>;
            #last if $next == $sum || is-zero($point) || is-zero(1 - $point);

            last if $next == $sum;
            $sum = $next;
        }

        return $sum
    }

    method !de-rule-estimate($region, $h, $oldSum = Nil) {
        my $step = $oldSum === Nil ?? 1 !! 2;

        my $s1 = self!fixed-point($region, 1, $h, $step);
        my $s2 = self!fixed-point($region, -1, $h, -$step);

        return $oldSum === Nil
                ?? $h * ($s1 + self!eval-integrand($region, $h, 0) + $s2)
                !! $h * ($s1 + $s2) + $oldSum / 2
    }

    method integrate($region) {

        my ($sum, $h, $k);
        $sum = $region.reuse-values<sum> // Nil;
        $h = $region.reuse-values<h> // 1;
        $k = $region.reuse-values<k> // 0;

        $!evaluations-count = 0;

        my $previous = self!de-rule-estimate($region, $h, Nil);
        my $new = $previous;

        $h /= 2;
        $new = self!de-rule-estimate($region, $h, $previous);
        $k++;

        $region.reuse-values = %(:$sum, :$h, :$k);

        self.integral = $new;
        self.error = ($new - $previous).abs;
        self.largest-error-axis = 0; # does not play a role of DoubleExponentialRule

        return self
    }

    # The "mini strategy" method integrate-full should take working precision
    # or tolerance as an argument, in order to figure out when to stop accumulating the terms sum.

    #| Full, on the spot integration. (I.e. a "mini" strategy.)
    method integrate-full($region, Numeric:D $step = 1, Numeric:D :$tolerance = 10e-10, UInt:D :$max-recursion = 7) {

        # Computing terms of the infinite range trapezoidal formula in both directions:
        # LaTeX: \int_{-\infty }^{+\infty } f(\phi (t)) \phi '(t) \, dx
        # ASCII: ∫ [-∞,+∞] f(φ(t)) φ'(t) dx
        # where phi is one of the transformations in ::VariableTransform::DoubleExponential

        my $h = 1;
        my $previous = self!de-rule-estimate($region, $h, Nil);
        my $temp = $previous;
        my $k = 0;

        loop {
            $h /= 2;
            $temp = self!de-rule-estimate($region, $h, $previous);
            $k++;

            last if $k >= $max-recursion;
            last unless abs($previous) * $tolerance <= abs($previous - $temp);

            $previous = $temp;
        }

        note "The DoubleExponential rule method integrate-full failed to converge to prescribed accuracy after $max-recursion interated refinements."
        if $k >= $max-recursion;

        $temp
    }

    method integrate-dim1($region, Numeric:D $step = 1) {!!!}

    # Essentially, dynamic Cartesian rule application
    method integrate-dimN($region, Numeric:D $step = 1) {!!!}
}
