use v6.d;

use Math::NIntegrate::Strategy;
use Math::NIntegrate::VariableTransformer::IMT;
use Math::NIntegrate::Utilities;
use Math::NIntegrate::Codes;

class Math::NIntegrate::Strategy::LocalAdaptive
        is Math::NIntegrate::Strategy {
    my %no-result = integral => Whatever, error => Whatever;

    has $.partitioning = Whatever;
    has Bool:D $.initial-estimate-relaxation = True;


    #| LocalAdaptive strategy's algorithm
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

        # Estimate the precision tolerance based on the user-specified precision tolerance
        # and the first integral estimate.
        my $eps = $MACHINE_EPSILON;
        # I am not sure what a good values here for Rat and FatRat
        my $est-rel-tol = $relative-tolerance / $eps;
        my $rel-est = $integral * $est-rel-tol;

        my $acc-est = $relative-tolerance / $eps;

        my $reference-estimate = max($rel-est, $acc-est);

        # Maybe the integral estimate is zero -- this does not mean the integral is zero.
        # Use the sum of the volumes of the regions.
        if is-zero($integral) {
            $reference-estimate = [+] self.regions>>.volume;
        }

        # Integration monitor
        with &integration-monitor {
            try &integration-monitor(self.regions);
            warn 'Cannot apply integration monitor.' if $!
        }

        # Recursive computation
        my @rec-res = self.regions.map({
            self.recursive-step(
                    region => $_,
                    :$reference-estimate,
                    :$relative-tolerance,
                    :$absolute-tolerance,
                    :$working-precision,
                    :&integration-monitor)
        });

        # Result
        return {
            integral => @rec-res.map(*<integral>).sum,
            error => @rec-res.map(*<error>).sum,                # sum or max?
            region-count => @rec-res.map(*<region-count>).sum,
        }
    }

    #| LocalAdaptive recursive step
    method recursive-step(
            :$region,
            :$reference-estimate,
            :$relative-tolerance = 1e-6,
            :$absolute-tolerance = 0,
            :$working-precision = Num,
            :&integration-monitor = WhateverCode,
            -->Map:D) {

        # It might happen that the precision is exhausted.
        # Hence, check the middle points.
        # TBD...

        # The region is either post-integration from method algorithm,
        # or it is clone of an integrated region obtained by Region::divide.
        # Hance it has a largest error axis.
        # In both cases, the levels for all axes should be the same.
        my $axis = $region.axis;
        if $region.levels[$axis] == self.singularity-depth && $region.range-end-cases[$axis] ne RE_NONE {
            # Apply singularity handler
            if self.singularity-handler ~~ Str:D && self.singularity-handler.lc eq 'imt' {
                $region.add-variable-transformer('imt', :$axis, :$working-precision);
                # Prevent another application of a singularity handler on that axis
                $region.range-end-cases[$axis] = RE_NONE
            }

            # Reset split level
            $region.levels[$axis] = 0;
        }

        # Integrate
        try $region.apply-rule;
        return %no-result if $!;

        my $error = $region.error;
        my $integral = $region.integral;
        $axis = $region.axis;

        if is-zero(numerical($reference-estimate + $error, Num) - $reference-estimate) {
            # Cannot see the error
            return {:$integral, :$error, region-count => 1}
        } else {

            # Warning the max recursion was reached
            if $region.levels.max ≥ self.max-recursion {
                note "Failed to converge to prescribed accuracy after {self.max-recursion} recursive bisections in near {$region.numerical-function.last-argument-values}.";
                return {:$integral, :$error, region-count => 1}
            }

            # If $!partitioning = Whatever the region should partitioned
            # in order to reuse the integrand values, if the integration rule is closed.

            # If one of the transformers is a singularity handler the region is just split.
            # TBD...

            # Special handling of closed rules for integrand values reuse.
            # TBD...

            my @divisions = do given $!partitioning {
                when Whatever { 2 xx $region.dimension }
                when $_ ~~ Int:D && $_ > 1 { $_ xx $region.dimension }
                default {
                    die 'The partitioning option is expected to be an integer greater than 1 or Whatever.'
                }
            }

            # Divide the region
            my @regions = $region.divide(@divisions);

            # Reverse the variable if needed for the new regions
            for @regions -> $r {
                $r.add-variable-transformer('reverse', :$axis, :$working-precision)
                if $r.range-end-cases[$axis] eq RE_RIGHT;
            }

            # Recursive computation
            my %result = integral => 0, error => 0, region-count => 0;

            # Using a for loop in order to facilitate early bailout
            for @regions -> $region {
                my %recRes = self.recursive-step(
                        :$region,
                        :$reference-estimate,
                        :$relative-tolerance,
                        :$absolute-tolerance,
                        :&integration-monitor,
                        :$working-precision);

                return %no-result unless %recRes<integral> ~~ Numeric:D;

                %result<integral> += %recRes<integral>;
                %result<error> += %recRes<error>;
                %result<region-count> += %recRes<region-count>;
            }

            return %result
        }

        # Should not be reached
        return %no-result
    }
}
