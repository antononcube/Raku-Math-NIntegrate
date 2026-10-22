use v6.d;

use Math::NIntegrate::Strategy;
use Math::NIntegrate::Strategy::ErrorHandlish;
use Math::NIntegrate::Rule::DoubleExponential;
use LeftistHeap;

class Math::NIntegrate::Strategy::DoubleExponential
        is Math::NIntegrate::Strategy
        does Math::NIntegrate::Strategy::ErrorHandlish {

    #| DoubleExponential strategy's algorithm
    method algorithm(
            Numeric:D :tol(:$relative-tolerance) = 1e-6,
            Numeric:D :acc(:$absolute-tolerance) = 0,
            :$working-precision = Num,
            :&integration-monitor = WhateverCode
            -->Map:D) {

        # This algorithm is very similar to that of crude MonteCarlo.
        # The differences are:
        # (i) in the stopping criteria (integral and error estimates) and
        # (ii) initialization (partitioning and min-recursion, variable transformers).
        # (Maybe, also convergence monitoring when it is implemented.)
        # Hence, there is the idea to have ::Strategy::Saturation class.
        # But such class might make algorithms harder to perceive or understand
        # because of the delegated steps to the children (like this class).
        # Hence, each of these, conceptual at this point children has its own full algorithm.

        # At this point the working precision should be known in the integrand and region.
        # Do we want to reset the working precision here?

        # Integral dimension shortcut
        my $dim = self.regions.head.dimension;

        # Divide regions according to min-recursion
        self.min-recursion-regions;

        # Make sure all regions have DoubleExponential variable transformer
        # First integration step
        try self.regions>>.apply-rule;

        # For now the error is printed out.
        # Something more graceful has to be done later.
        #return %no-result if $!;
        if $! { note "{self.msgProbOrig}\n{$!}"; return self.no-result }

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
        my UInt:D $number-of-points = $heap.values.map(*.rule.abscissas.elems).sum;
        my $done-max-points = False;
        $done-max-points = $number-of-points > self.max-points if self.max-points ~~ Numeric:D;
        my $step;
        my $topRegion;

        # Main loop
        while !($done-tol || $done-accuracy || $done-max-recursion || $done-max-points) {
            $step++;

            # say ('start of loop:', :$integral, :$error, :$step, region-count => $heap.elems);

            # Delete from heap the region with largest error
            $topRegion = $heap.delete-top-element;

            $error -= $topRegion.error;
            $integral -= $topRegion.integral;
            my $axis = $topRegion.axis;

            # say ('after removing top region:', :$integral, :$error, :$axis, :$step, region-count => $heap.elems);

            # Integrate
            try $topRegion.apply-rule;
            if $! { note "{self.msgProbSplit}\n{$!}"; return self.no-result }

            # Convergence monitoring
            # TBD ...

            # Estimate sums
            $error += $topRegion.error;
            $integral += $topRegion.integral;
            for ^$dim -> $i { $topRegion.levels[$i] += 1}

            # Add the split regions to the heap
            $heap.insert($topRegion);

            # Needs to be properly implemented
            $number-of-points = $topRegion.rule.last-values.elems;

            # Integration monitor
            with &integration-monitor {
                try &integration-monitor($heap.values);
                warn 'Cannot apply integration monitor.' if $!
            }

            # Stopping criteria
            $done-tol = $error ≤ $relative-tolerance * $integral.abs;
            $done-accuracy = $error ≤ $absolute-tolerance;
            $done-max-recursion = $topRegion.levels[$axis] > self.max-recursion;
            $done-max-points = $number-of-points > self.max-points if self.max-points ~~ Numeric:D;

            # say ('end of loop:', :$integral, :$error, relative-error => $error/$integral, region-count => $heap.elems, :$step)
        }

        # Warning the max recursion was reached
        if $done-max-recursion {
            note self.max-iterative-refinements-reached($topRegion.min, $topRegion.max, self.max-recursion)
        }

        if $done-max-points {
            note "Failed to converge to prescribed accuracy after {self.max-points} integrand evaluations."
        }

        # Put the regions in the heap in the object regions holder
        self.regions = $heap.values;

        # Result
        return %(:$integral, :$error, region-count => self.regions.elems)
    }
}
