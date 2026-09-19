use v6.d;

use Math::NIntegrate::Region;

class Math::NIntegrate::Strategy {

    # Do we need these here or they should be arguments to the method algorithm?
    #has Numeric $.absolute-tolerance;
    #has Numeric $.relative-tolerance;

    # Precision epsilon
    has Numeric $.epsilon;

    # Current integral estimate
    has Numeric $.integral-estimate;

    # Options
    has Numeric $.singularity-depth;
    has UInt $.min-recursion;
    has UInt $.max-recursion;
    has $.max-points;
    has $.max-number-of-error-increases;

    # State
    has UInt $.error-increase-count;

    # Regions
    has Math::NIntegrate::Region @.regions is rw;

    #======================================================
    # Creators
    #======================================================
    submethod BUILD(
            UInt :$!min-recursion = 0,
            UInt :$!max-recursion = 12,
            :$!max-points = Whatever,
            :$!singularity-depth = 4,
            :$!max-number-of-error-increases = Whatever,
            :@!regions
                    ) {
        without @!regions {
            fail 'MISSING_OBJECT: at least one region is expected for integration startegy creation.'
        }

        # :12max-recursion is going to produce combinatorial explosion in high dimensional integrals.
        # Hence, it can be Whatever and determined by regions' dimension.

        # The option max-points is "soft", i.e., hard to respect precisely.
        # Being Whatever means that it is ignored.
    }

    #======================================================
    # Methods
    #======================================================

    #| Template method
    method min-recursion-regions(-->Array:D) {
        # Divide -- not split -- the regions $!min-recursion number of times
        warn 'Processing of regions by min-recursion are is not implemented yet.';
        return @!regions;
    }

    #| Strategy's initialization
    method intialize(-->Math::NIntegrate::Strategy) {!!!}

    #| Strategy's stopping criteria.
    #| Returns a hashmap with Boolean values of different stopping criteria:
    #| <precision accuracy max-recursion max-points max-number-of-error-increases>
    method stopping-criteria(-->Map:D) {!!!}

    #| Strategy's algorithm
    method algorithm(
            Numeric:D :$relative-tolerance = 1e-6,
            Numeric:D :$absolute-tolerance = 0,
            :$working-precision = Num
            -->Map:D) {!!!}
}
