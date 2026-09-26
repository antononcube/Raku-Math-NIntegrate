use v6.d;

use LeftistHeap;
use Math::NIntegrate::Strategy;
use Math::NIntegrate::VariableTransformer::IMT;
use Math::NIntegrate::Codes;

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
        self.min-recursion-regions;

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

            # Application of singularity handler or region splitting
            if $topRegion.levels[$axis] == self.singularity-depth && $topRegion.range-end-cases[$axis] ne RE_NONE  {
                # Apply singularity handler
                if self.singularity-handler ~~ Str:D && self.singularity-handler.lc eq 'imt' {
                    $topRegion.add-variable-transformer('imt', :$axis, :$working-precision)
                }

                # Reset split level
                $topRegion.levels[$axis] = 0;

                # Integrate
                try $topRegion.apply-rule;
                return %no-result if $!;

                # Estimate sums
                $error += $topRegion.error;
                $integral += $topRegion.integral;

                # Add the region with a new variable transformer to the heap
                $heap.insert($topRegion);

            } else {
                # Split the region
                my $newRegion = $topRegion.split(:$axis);

                # Reverse the variable if needed for the new region
                $newRegion.add-variable-transformer('reverse', :$axis, :$working-precision)
                if $newRegion.range-end-cases[$axis] eq RE_RIGHT;

                # Not that during the region splitting the new region can marked as a "middle" region
                # with the enum value RE_NONE and that prevents IMT not be applied to it.

                # Integrate
                try $topRegion.apply-rule;
                return %no-result if $!;

                try $newRegion.apply-rule;
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
            }

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
