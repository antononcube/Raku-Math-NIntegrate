use v6.d;

use Math::NIntegrate::Strategy::Decorator;

class Math::NIntegrate::Strategy::UnitCubeRescaling
        is Math::NIntegrate::Strategy::Decorator {

    # Whether the unit-cube rescaling is applied already or not.
    # I am not sure that it is needed since unit-cube rescaling an idempotent.
    has Bool:D $!rescaled = False;

    # This decorator class replaces the integrands and ranges of component's regions.
    # Since regions can have different integrands and/or ranges,
    # each region is treated separately.

    # In principle if at creation the component has regions, they can be modified at creation.
    # But it is more flexible to it as an added behavior.
    #submethod TWEAK(:$component, :@ranges) {
    #}

    # This can be an utility.
    multi method ranges-to-cube(
            @ranges,
            @cube-sides where @cube-sides.all ~~ Numeric:D) {
        self.ranges-to-cube(@ranges, @cube-sides xx @ranges.elems)
    }

    multi method ranges-to-cube(
            @ranges,
            @cube-sides where @cube-sides.all ~~ Positional:D) {
        die 'The number ranges and cube sides is expected to be the same.'
        unless @ranges.elems == @cube-sides.elems;

    }

    method rescale() {
        # For each region transform the integrand and ranges
        my @newRegions = do for self.component.regions -> $reg {

        }

        self.component.regions = |@newRegions;
    }

    method algorithm(*%args --> Map:D) {
        # Conditional rescaling should not be used since since the regions are rw in Math::NIntegrate::Strategy.
        # Unit cube rescaling, though, is idempotent, hence it can be bravely applied/invoked.
        #self.rescale unless $!rescaled;
        #self.rescale;
        note 'Call in Math::NIntegrate::Strategy::UnitCubeRescaling::algorithm';
        self.Math::NIntegrate::Strategy::Decorator::algorithm(|%args)
    }
}
