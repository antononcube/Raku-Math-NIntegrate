use v6.d;

use Math::NIntegrate::Strategy;
use Math::NIntegrate::Region;
use Math::NIntegrate::NumericalFunction;
use Math::NIntegrate::Rule;
use Math::NIntegrate::VariableTransformer::Composite;

class Math::NIntegrate::Builder {

    has Math::NIntegrate::Strategy $.strategy;

    method make-region(Math::NIntegrate::NumericalFunction:D $numerical-function, @min, @max, :$working-precision = Num) {

        fail 'DIMENSIONS_DO_NOT_MATCH: region boundaries dimensions do not match.'
        if @min.elems != @max.elems;

        # First a region object is made
        my $region = Math::NIntegrate::Region.new(:$numerical-function, :@min, :@max);

        # Variable transformer (Composite) is made with $region
        my $var-trans = Math::NIntegrate::VariableTransformer::Composite.new(:$region);

        # Variable transformer object is set in $region
        $region.variable-transformer = $var-trans;

        return $region
    }

    method make-numerical-function(&function, :$working-precision = Num) {

        # Numerical function
        my $nf = Math::NIntegrate::NumericalFunction(:&function, $working-precision);

        return $nf
    }

    method make-rules() {

        # At this point $!strategy has regions

        # Integration rule is created according to specified strategy name and integration rule
        # $rule = ...

        # Integration rule object is set to all regions
    }
}
