use v6.d;

#| Error message and no-result role
role Math::NIntegrate::Strategy::ErrorHandlish {

    has %.no-result = integral => Whatever, error => Whatever;

    has Bool:D $.has-reached-max-recursion is rw = False;

    has $.slow-convergence-detection-split-count = 7;
    has UInt:D $.slow-convergence-msg-count is rw = 0;

    has $.msgProbOrig = 'Problems integrating the original set of regions.';
    has $.msgProbSplit = 'Problems integrating the split regions.';
    has $.msgProbIMT = 'Problems integrating the region after IMT application.';
    has $.msgProbRegion = 'Problems integrating the region at level';
    has $.msgSuppress = 'Suppressing further messages of this type.';

    method max-iterative-refinements-reached(@min, @max, $max-recursion) {
        my @range-info = |(@min>>.Numeric Z @max>>.Numeric).map({ "[{$_.head}, {$_.tail}]" });
        "Failed to converge to prescribed accuracy after {$max-recursion} " ~
                "recursive refinement{ $max-recursion == 1 ?? '' !! 's'} " ~
                "in the region {@range-info.join(' x ')}."
    }

    method max-recursion-reached(@args, $max-recursion) {
        # Maybe this should be set explicitly in the strategy object.
        # This is needed in order to "immediate" stop in LocalAdaptive.
        self.has-reached-max-recursion = True;
        "Failed to converge to prescribed accuracy after $max-recursion recursive bisections in near {@args.raku}."
    }

    method max-points-reached($max-points) {
        "Failed to converge to prescribed accuracy after $max-points integrand evaluations."
    }

    method slow-convergence-detected() {
        self.slow-convergence-msg-count += 1;
        # When high precision is fully implemented the message should also say "or working-precision too small."
        "Numerical integration converging too slowly; suspect one of the following: singularity, value of the integration is 0, or highly oscillatory integrand."
    }
}
