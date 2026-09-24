use v6.d;

use Math::NIntegrate::Strategy;
use LeftistHeap;

class Math::NIntegrate::Strategy::GlobalAdaptive
        is Math::NIntegrate::Strategy {

    my %no-result = integral => Whatever, error => Whatever;

    #| GlobalAdaptive strategy's algorithm
    method algorithm(
            Numeric:D :tol(:$relative-tolerance) = 1e-6,
            Numeric:D :acc(:$absolute-tolerance) = 0,
            :$working-precision = Num,
            :&integration-monitor = WhateverCode
            -->Map:D) {

        # At this point the working precision should be known in the integrand and region.
        # Do we want to reset the working precision here?

        # Integral dimension shortcut
        my $dim = self.regions.head.dimension;

        # Divide regions according to min-recursion
        # self.min-recursion-regions;

        # First integration step
        try self.regions>>.apply-rule;
        return %no-result if $!;

        my $error = self.regions.map(*.error).sum;
        my $integral = self.regions.map(*.integral).sum;

        # Make the heap of regions
        my $heap = LeftistHeap.new(self.regions, comparator => { $^a.error > $^b.error });

        # Integration monitor
        with &integration-monitor {
            try &integration-monitor($heap.values);
            warn 'Cannot apply integration monitor.' if $!
        }

        # Criteria variables and their first values
        my Bool:D $done-tol = $error ≤ $relative-tolerance * $integral.abs;
        my Bool:D $done-accuracy = $error ≤ $absolute-tolerance;
        my Bool:D $done-max-recursion = $heap.top.levels.max > self.max-recursion;
        my $step;
        my $topRegion;

        # Main loop
        while !($done-tol || $done-accuracy || $done-max-recursion) {
            $step++;

            #say ('start of loop:', :$integral, :$error, :$step, region-count => $heap.elems);
            # Delete from heap the region with largest error
            $topRegion = $heap.delete-top-element;


            $error -= $topRegion.error;
            $integral -= $topRegion.integral;
            my $axis = $topRegion.axis;

            #say ('after removing top region:', :$integral, :$error, :$step, region-count => $heap.elems);

            # Application of singularity handler
            if $topRegion.levels[$axis] == self.singularity-depth {
                # Apply singularity handler
            }

            # Split the region
            my $newRegion = $topRegion.split(:$axis);

            # Reverse the variable if needed for the new region
            # TBD...

            # Integrate
            $topRegion.apply-rule;
            return %no-result if $!;

            $newRegion.apply-rule;
            return %no-result if $!;

            # Convergence monitoring
            # TBD ...

            # Estimate sums
            $error += $topRegion.error + $newRegion.error;
            $integral += $topRegion.integral + $newRegion.integral;

            # Instead of compensated summation
            # LeftistHeap has a traverse method which can be used to get error and integral estimates
            # without making a new array of regions, only a new array of scalar values.
            #$error = [|$heap.values.map(*.error), $topRegion.error, $newRegion.error].sort(*.abs).sum;
            #$integral = [|$heap.values.map(*.integral), $topRegion.integral, $newRegion.integral].sort(*.abs).sum;

            # Add the split regions to the heap
            $heap.insert($topRegion);
            $heap.insert($newRegion);

            # Integration monitor
            with &integration-monitor {
                try &integration-monitor($heap.values);
                warn 'Cannot apply integration monitor.' if $!
            }

            # Stopping criteria
            $done-tol = $error ≤ $relative-tolerance * $integral.abs;
            $done-accuracy = $error ≤ $absolute-tolerance;
            $done-max-recursion = $topRegion.levels[$axis] > self.max-recursion;

            #say ('end of loop:', :$integral, :$error, relative-error => $error/$integral, region-count => $heap.elems, :$step)
        }

        # Warning the max recursion was reached
        if $done-max-recursion {
            note "Failed to converge to prescribed accuracy after {self.max-recursion} recursive bisections in near {$topRegion.numerical-function.last-argument-values}."
        }

        # Put the regions in the heap in the object regions holder
        self.regions = $heap.values;

        # Result
        return %(:$integral, :$error, region-count => self.regions.elems)
    }
}
