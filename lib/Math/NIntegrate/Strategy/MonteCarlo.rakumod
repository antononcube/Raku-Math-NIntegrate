use v6.d;

# Inheriting of Saturating is to be done
#use Math::NIntegrate::Strategy::Saturating;
use Math::NIntegrate::Strategy;
use LeftistHeap;

class Math::NIntegrate::Strategy::MonteCarlo
        is Math::NIntegrate::Strategy {

    my %no-result = integral => Whatever, error => Whatever;

    has $.partitioning = Whatever;
    has $.random-seed = Whatever;

    #| MonteCarlo strategy's algorithm
    method algorithm(
            Numeric:D :tol(:$relative-tolerance) = 1e-6,
            Numeric:D :acc(:$absolute-tolerance) = 0,
            :$working-precision = Num,
            :&integration-monitor = WhateverCode
            -->Map:D) {
        if $!random-seed ~~ Numeric:D {
            srand($!random-seed.round)
        }

        # At this point the working precision should be known in the integrand and region.
        # But does high precision matter for Monte-Carlo methods.

        # Integral dimension shortcut
        my $dim = self.regions.head.dimension;

        # Does the partitioning option override max-recursion option?
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
        my UInt:D $number-of-points = $heap.values.map(*.rule.abscissas.elems).sum;
        my $done-max-points = $number-of-points > self.max-points;
        my $step;
        my $topRegion;

        # Main loop
        while !($done-tol || $done-accuracy || $done-max-recursion || $done-max-points) {
            $step++;

            #say ('start of loop:', :$integral, :$error, :$step, region-count => $heap.elems);
            # Delete from heap the region with largest error
            $topRegion = $heap.delete-top-element;


            $error -= $topRegion.error;
            $integral -= $topRegion.integral;
            my $axis = $topRegion.axis;

            #say ('after removing top region:', :$integral, :$error, :$axis, :$step, region-count => $heap.elems);

            # More points to the top region

            # Integrate
            try $topRegion.apply-rule;
            return %no-result if $!;

            # Estimate sums
            $error += $topRegion.error;
            $integral += $topRegion.integral;

            # Add the region with a new variable transformer to the heap
            $heap.insert($topRegion);

            # Integration monitor
            with &integration-monitor {
                try &integration-monitor($heap.values);
                warn 'Cannot apply integration monitor.' if $!
            }

            # Stopping criteria
            $done-tol = $error ≤ $relative-tolerance * $integral.abs;
            $done-accuracy = $error ≤ $absolute-tolerance;
            $done-max-recursion = $topRegion.levels[$axis] > self.max-recursion;

            $number-of-points += $topRegion.rule.abscissas.elems;
            $done-max-points = $number-of-points > self.max-points;

            #say ('end of loop:', :$integral, :$error, relative-error => $error/$integral, region-count => $heap.elems, :$step)
        }

        # Warning the max recursion was reached
        if $done-max-recursion {
            note "Failed to converge to prescribed accuracy after {self.max-recursion} recursive bisections in near {$topRegion.numerical-function.last-argument-values}."
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
