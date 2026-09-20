use v6.d;

use Math::NIntegrate::VariableTransformer;

# It would be nice to just reverse the variables at given indexes and be "done with it".
# For now we follow the typical variable transformer application.
# This is for top-level integration region usage/manipulation.
# The internal variable inversion is handled by Infinity.

class Math::NIntegrate::VariableTransformer::Reverse
        is Math::NIntegrate::VariableTransformer {

    method clone(-->Math::NIntegrate::VariableTransformer::Reverse) {
        # fail 'INVOKING_TEMP_OBJECT cloning reverse variable transformer is meaningless';
        Math::NIntegrate::VariableTransformer::Reverse.new(region => self.region).copy(self, :clone)
    }

    submethod TWEAK(*%args) {
        my @indexes = %args<indexes> // %args<indices>;
    }
}
