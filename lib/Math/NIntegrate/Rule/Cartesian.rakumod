use v6.d;

use Math::NIntegrate::Rule::General;

# Composite design pattern application over 1D integration rules.
class Math::NIntegrate::Rule::Cartesian
        is Math::NIntegrate::Rule::General {

    has @.components;

    #======================================================
    # Components management methods
    #======================================================
    #| Add a transformer
    method add(Math::NIntegrate::Rule:D $obj) {
        # Should we check for the rule 1D dimensional?
        @!components.push($obj);
        return self
    }

    #| Remove a transformer object
    multi method remove(Math::NIntegrate::Rule:D $obj) {
        @!components .= grep({ $_ !=== $obj });
        return self
    }

    #| Remove a transformer at give index
    multi method remove(Int:D $index) {
        @!components = @!components.splice($index, 1);
        return self
    }

    #======================================================
    # Integration
    #======================================================

    method integrate($region) {
        my $integralLocal = 0;
        my $errorLocal = 0;
        my @values;

        for @!components -> $rule {

        }

        return self;
    }
}
