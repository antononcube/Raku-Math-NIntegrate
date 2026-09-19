use v6.d;

use Math::NIntegrate::Strategy;
use LeftistHeap;

class Math::NIntegrate::Strategy::GlobalAdaptive
        is Math::NIntegrate::Strategy {

    #| GlobalAdaptive strategy's algorithm
    method algorithm(
            Numeric:D :$relative-tolerance = 1e-6,
            Numeric:D :$absolute-tolerance = 0,
            :$working-precision = Num
            -->Map:D) {

        # At this point the working precision should be known in the integrand and region.
        # Do we want to reset the working precision here?

        # Integral dimension shortcut
        my $dim = self.regions.head.dimension;

        # Divide regions according to min-recursion
        self.regions = self.min-recursion.regions();

        # First integration step
        self.regions>>.apply-rule;
        my $error = self.regions.map(*.error).sum;
        my $integral = self.regions.map(*.integral).sum;

        # Make the heap of regions
        my $heap = LeftistHeap.new(self.regions, comparator => { $^a > $^b });

        my Bool:D $done-tol = $error ≤ $relative-tolerance * $integral.abs;
        my Bool:D $done-accuracy = $error ≤ $absolute-tolerance;
        my Bool:D $done-max-recursion = $heap.top.level > self.max-recursion;
        my $step;
        my $topRegion;
        while !($done-tol || $done-accuracy || ) {
            $step++;

            # Delete from heap the region with largest error
            $topRegion = $heap.delete-top-element;

            $error -= $topRegion.error;
            $integral -= $topRegion.integral;
            my $axis = $topRegion.axis;

            # Application of singularity handler
            if $topRegion.level == self.singularity-depth {
                # Apply singularity handler
            }

            # Split the region
            my $newRegion = $topRegion.split(:$axis);

            # Reverse the variable if needed for the new region
            # TBD...

            # Integrate
            $topRegion.apply-rule;
            $newRegion.apply-rule;

            # Convergence monitoring
            # TBD ...

            # Estimate sums
            $error += $topRegion.error + $newRegion.error;
            $integral += $topRegion.integral + $newRegion.integral;

            # Add the split regions to the heap
            $heap.insert($topRegion);
            $heap.insert($newRegion);

            # Integration monitor
            # TBD...

            # Stopping criteria
            $done-tol = $error ≤ $relative-tolerance * $integral.abs;
            $done-accuracy = $error ≤ $absolute-tolerance;
            $done-max-recursion = $topRegion.level > self.max-recursion;
        }

        # Warning the max recursion was reached
        if $done-max-recursion {
            note "Failed to converge to prescribed accuracy after {self.max-recursion} recursive bisections in near {$topRegion.numerical-function.last-argument-values}."
        }

        # Put the regions in the heap in the object regions holder
        self.regions = $heap.values;

        # Result
        return %(:$integral, :$error, number-of-regions => self.regions.elems)
    }
}
