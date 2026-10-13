use v6.d;

use LeftistHeap;
use Math::NIntegrate::Strategy;
use Math::NIntegrate::Strategy::ErrorHandlish;
use Math::NIntegrate::Spec;
use Math::NIntegrate::Codes;

class Math::NIntegrate::Strategy::AdaptiveMonteCarlo
        is Math::NIntegrate::Strategy
        does Math::NIntegrate::Strategy::ErrorHandlish {

    has $.partitioning = Whatever;
    has $.random-seed = Whatever;
    has $.bisection-dithering = 0;

    submethod TWEAK(:$!partitioning, :$!random-seed, *%args) {

        die 'The value of $random-seed is expected to be an integer or Whatever.'
        unless $!random-seed ~~ Int:D || $!random-seed.isa(Whatever);

        note 'The (crude) Monte Carlo strategy does not use min-recursion.'
        unless %args<min-recursion>.isa(Whatever) || %args<min-recursion> ~~ Numeric:D && %args<min-recursion> == 0;
    }

    #| GlobalAdaptive strategy's algorithm
    method algorithm(
            Numeric:D :tol(:$relative-tolerance) = 1e-2,
            Numeric:D :acc(:$absolute-tolerance) = 0,
            :$working-precision = Num,
            :&integration-monitor = WhateverCode
            -->Map:D) {
        if $!random-seed ~~ Numeric:D {
            srand($!random-seed.round)
        }

        # At this point the working precision should be known in the integrand and region.
        # For Monte Carlo high precision does not make sense.
        # Should the precision be reset to Num?

        # Integral dimension shortcut
        my $dim = self.regions.head.dimension;

        say (:$dim, :$!partitioning, :$!bisection-dithering);

        # The min-recursion option is has to be taken into account.
        # Should min-recursion be applied first and then partitioning or the other way around?
        # Is only one of them applied?
        # Divide regions according to partitioning option (normalized at this point)
        $!partitioning = Math::NIntegrate::Spec.normalize-partitioning($!partitioning, :$dim);
        self.regions = |self.regions.map({ $_.divide($!partitioning) }).flat(:hammer) with $!partitioning;

        # Divide regions according to min-recursion
        self.min-recursion-regions;

        # First integration step
        try self.regions>>.apply-rule;

        # For now the error is printed out.
        # Something more graceful has to be done later.
        #return %no-result if $!;
        if $! { note "{self.msgProbOrig}\n{$!}"; return self.no-result }

        my $error = self.regions.map(*.error ** 2).sum;
        my $integral = self.regions.map(*.integral).sum;

        # Make the heap of regions
        my $heap = LeftistHeap.new(self.regions, comparator => { $^a.error > $^b.error });

        # Integration monitor
        with &integration-monitor {
            try &integration-monitor($heap.values);
            warn 'Cannot apply integration monitor.' if $!
        }

        # Criteria variables and their first values
        my Bool:D $done-tol = $error ≤ ($relative-tolerance * $integral.abs) ** 2;
        my Bool:D $done-accuracy = $error ≤ $absolute-tolerance ** 2;
        my Bool:D $done-max-recursion = $heap.top.levels.max > self.max-recursion;
        my UInt:D $number-of-points = $heap.values.map(*.rule.abscissas.elems).sum;
        my $done-max-points = False;
        $done-max-points = $number-of-points > self.max-points if self.max-points ~~ Numeric:D;
        my $step;
        my $topRegion;

        # Main loop
        while !($done-tol || $done-accuracy || $done-max-recursion || $done-max-points) {
            $step++;

            #say ('start of loop:', :$integral, :$error, :$step, region-count => $heap.elems);
            # Delete from heap the region with largest error
            $topRegion = $heap.delete-top-element;

            $error -= $topRegion.error ** 2;
            $integral -= $topRegion.integral;
            my $axis = $topRegion.axis;

            #say ('after removing top region:', :$integral, :$error, :$axis, :$step, region-count => $heap.elems);

            # Split the region
            my $newRegion = $topRegion.split(:$axis);

            # No reuse
            $topRegion.reuse-values = %();
            $newRegion.reuse-values = %();

            # Integrate
            try $topRegion.apply-rule;
            if $! { note "{self.msgProbSplit}\n{$!}"; return self.no-result }

            try $newRegion.apply-rule;
            if $! { note "{self.msgProbSplit}\n{$!}"; return self.no-result }

            # Convergence monitoring
            # TBD ...

            # Estimate sums
            $error += $topRegion.error ** 2 + $newRegion.error ** 2;
            $integral += $topRegion.integral + $newRegion.integral;

            # Instead of compensated summation
            # LeftistHeap has a traverse method which can be used to get error and integral estimates
            # without making a new array of regions, only a new array of scalar values.
            #$error = [|$heap.values.map(*.error), $topRegion.error, $newRegion.error].sort(*.abs).sum;
            #$integral = [|$heap.values.map(*.integral), $topRegion.integral, $newRegion.integral].sort(*.abs).sum;

            # Add the split regions to the heap
            $heap.insert($topRegion);
            $heap.insert($newRegion);

            $number-of-points += $topRegion.rule.abscissas.elems;
            $number-of-points += $newRegion.rule.abscissas.elems;

            # Integration monitor
            with &integration-monitor {
                try &integration-monitor($heap.values);
                warn 'Cannot apply integration monitor.' if $!
            }

            # Stopping criteria
            $done-tol = $error ≤ ($relative-tolerance * $integral.abs) ** 2;
            $done-accuracy = $error ≤ $absolute-tolerance ** 2;
            $done-max-recursion = $topRegion.levels[$axis] > self.max-recursion;
            $done-max-points = $number-of-points > self.max-points if self.max-points ~~ Numeric:D;

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
