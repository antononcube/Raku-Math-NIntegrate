use v6.d;

# Inheriting of Saturating is to be done
#use Math::NIntegrate::Strategy::Saturating;
use Math::NIntegrate::Strategy;
use Math::NIntegrate::Strategy::ErrorHandlish;
use Math::NIntegrate::Spec;
use LeftistHeap;

class Math::NIntegrate::Strategy::MonteCarlo
        is Math::NIntegrate::Strategy
        does Math::NIntegrate::Strategy::ErrorHandlish {

    has $.partitioning = Whatever;
    has $.random-seed = Whatever;

    submethod TWEAK(:$!partitioning, :$!random-seed, *%args) {

        die 'The value of $random-seed is expected to be an integer or Whatever.'
        unless $!random-seed ~~ Int:D || $!random-seed.isa(Whatever);

        note 'The (crude) Monte Carlo strategy does not use min-recursion.'
        unless %args<min-recursion>.isa(Whatever) || %args<min-recursion> ~~ Numeric:D && %args<min-recursion> == 0;
    }

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
        # But does high precision matter for Monte-Carlo methods?

        # Integral dimension shortcut
        my $dim = self.regions.head.dimension;

        # The min-recursion option is not respected, but should some warning be given if it is larger than 0?
        # Divide regions according to partitioning option (normalized at this point)
        $!partitioning = Math::NIntegrate::Spec.normalize-partitioning($!partitioning, :$dim);
        self.regions = |self.regions.map({ $_.divide($!partitioning) }).flat(:hammer) with $!partitioning;

        # First integration step
        try self.regions>>.apply-rule;
        if $! { say "{self.msgProbOrig}\n{$!}"; return self.no-result }

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

            # Delete from the heap the region with largest error
            $topRegion = $heap.delete-top-element;

            $error -= $topRegion.error;
            $integral -= $topRegion.integral;
            my $axis = $topRegion.axis;

            # More points to the top region -- see how MonteCarloRule reuses values.
            # Integrate
            try $topRegion.apply-rule;
            if $! { note "{self.msgProbRegion} {$topRegion.levels.max}\n{$!}"; return self.no-result }

            # Estimate sums
            $error += $topRegion.error;
            $integral += $topRegion.integral;
            for ^$dim -> $i { $topRegion.levels[$i] += 1}

            # Add the (former) top region to the heap
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
        }

        # Warning the max recursion was reached
        if $done-max-recursion {
            my @range-info = |($topRegion.min>>.Numeric Z $topRegion.max>>.Numeric).map({ "[{$_.head}, {$_.tail}]" });
            note "Failed to converge to prescribed accuracy after {self.max-recursion} " ~
                    "recursive refinements in region with range{ @range-info == 1 ?? '' !! 's' } {@range-info.join(' x ')}."
        }

        # Warning the max points was reached
        if $done-max-points {
            note "Failed to converge to prescribed accuracy after {self.max-points} integrand evaluations."
        }

        # Put the regions in the heap in the object regions holder
        self.regions = $heap.values;

        # Result
        return %(:$integral, :$error, region-count => self.regions.elems)
    }
}
