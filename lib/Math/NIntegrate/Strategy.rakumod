use v6.d;

use Math::NIntegrate::Region;
use LeftistHeap;

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
    method min-recursion-regions(-->Math::NIntegrate::Strategy) {

        without @!regions {
            fail 'MISSING_OBJECT: no regions to partition according to min-recursion spec.'
        }

        return self if $!min-recursion == 0;

        # Get dimension
        my $dim = @!regions.head.dimension;

        # Splitting each dimension $!min-recursion times is the
        # how min-recursion spec is interpreted:
        # :2min-recursion on a 2D region is going to produce 16 (4 x 4) subregions;
        # :2min-recursion on a 3D region is going to produce 64 (4 x 4 x 4) subregions.

        if (2 ** $!min-recursion) ** $dim > 1024 {
            # Preventing combinatorial explosion
            $!min-recursion = max((1024 ** (1/$dim)).log(2).floor, 1);
            note "Combinatorial explosion with the specified min-recursion; using {$!min-recursion} instead."
        }

        # Using region splitting for now -- simple, elegant and slow.
        # Using partitioning should be considered.
        my $heap = LeftistHeap.new(@!regions, comparator => { $^a.levels.min < $^b.levels.min });
        while $heap.top.levels.min < $!min-recursion {
            my $reg = $heap.delete-top-element;
            my $newReg = $reg.split(axis => $reg.levels.min(:k).head);
            $heap.insert($reg);
            $heap.insert($newReg);
        }

        @!regions = $heap.values;

        return self
    }

    #| Strategy's initialization
    method intialize(-->Math::NIntegrate::Strategy) {!!!}

    #| Strategy's stopping criteria.
    #| Returns a hashmap with Boolean values of different stopping criteria:
    #| <precision accuracy max-recursion max-points max-number-of-error-increases>
    method stopping-criteria(-->Map:D) {!!!}

    #| Strategy's algorithm
    method algorithm(
            Numeric:D :tol(:$relative-tolerance) = 1e-6,
            Numeric:D :acc(:$absolute-tolerance) = 0,
            :$working-precision = Num
            -->Map:D) {!!!}
}
