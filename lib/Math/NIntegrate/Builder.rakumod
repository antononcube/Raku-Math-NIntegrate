use v6.d;

use Math::NIntegrate::Strategy;
use Math::NIntegrate::Region;
use Math::NIntegrate::NumericalFunction;
use Math::NIntegrate::Rule;
use Math::NIntegrate::VariableTransformer::Composite;

class Math::NIntegrate::Builder {

    has Math::NIntegrate::Strategy $.strategy;

    method make-variable-transformer-composite(
            :$region
            --> Math::NIntegrate::VariableTransformer::Composite) {

        my $obj = Math::NIntegrate::VariableTransformer::Composite.new(:$region);

        # Should the stack have at least one variable transformer?
        # The Infinity transform can be seen as finite-range proxy of the region.
        # The Affine transform is for mapping or the finite boundaries to integration rules [0, 1] abscissas ranges.

        my @min-original-bounds = $region.min;
        my @max-original-bounds = $region.max;
        my $vtInf = Math::NIntegrate::VariableTransformer::Infinity.new(:@min-original-bounds, :@max-original-bounds);
        $obj.add($vtInf);
        return $obj
    }

    method make-region(
            Math::NIntegrate::NumericalFunction:D $numerical-function,
            @min,
            @max,
            :$rule = Nil,
            :$working-precision = Num
            --> Math::NIntegrate::Region) {

        fail 'DIMENSIONS_DO_NOT_MATCH: region boundaries dimensions do not match.'
        if @min.elems != @max.elems;

        # First a region object is made
        my $region = Math::NIntegrate::Region.new(:$numerical-function, :@min, :@max, :$rule, dimension => @min.elems);

        # Variable transformer (Composite) is made with $region
        my $var-trans = self.make-variable-transformer-composite(:$region);

        # Variable transformer object is set in $region
        $region.variable-transformer = $var-trans;

        return $region
    }

    method make-numerical-function(&function, :$working-precision = Num) {

        # Numerical function
        my $nf = Math::NIntegrate::NumericalFunction(:&function, $working-precision);

        return $nf
    }

    method make-rule() {

        # At this point $!strategy has regions

        # Integration rule is created according to specified strategy name and integration rule
        # $rule = ...

        # Integration rule object is set to all regions
    }
}
