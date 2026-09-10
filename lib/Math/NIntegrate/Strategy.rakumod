use v6.d;

class Math::NIntegrate::Strategy {

    has Numeric $.absolute-tolerance;
    has Numeric $.relative-tolerance;

    # Precision epsilon
    has Numeric $.epsilon;

    # Current integral estimate
    has Numeric $.integral-estimate;

    # Options
    has $.singularity-depth;
    has $.min-recursion;
    has $.max-recursion;
    has $.max-points;
    has $.max-number-of-error-increases;

    # State
    has $.error-increase-count;

    # Regions, all Math::NIntegrate::Region
    has @.regions;

    #======================================================
    # Methods
    #======================================================

    method algorithm(-->Math::NIntegrate::Strategy) {!!!}
}