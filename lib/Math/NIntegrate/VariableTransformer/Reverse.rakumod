use v6.d;

use Math::NIntegrate::VariableTransformer;

# It would be nice to just reverse the variables at given indexes and be "done with it" --
# see the Region.add-variable-transformer
# Below the typical variable transformer application is followed.

# In principle @!transformers can be used but using @!indexes is cleaner.

# This class is for top-level integration region usage/manipulation.
# The internal variable inversion is handled by Infinity and Region.

class Math::NIntegrate::VariableTransformer::Reverse
        is Math::NIntegrate::VariableTransformer {

    # Indexes of variables to revert
    has @.indexes;

    method clone(-->Math::NIntegrate::VariableTransformer::Reverse) {
        # fail 'INVOKING_TEMP_OBJECT cloning reverse variable transformer is meaningless';
        Math::NIntegrate::VariableTransformer::Reverse.new(region => self.region).copy(self, :clone)
    }

    submethod TWEAK(*%args) {
        @!indexes = %args<indexes> // %args<indices> // [];

        # There should be a check that the indexes are with regions dimensions.

        # Should there be a check that the original and transform boundaries are the same?
        # In general, those are ignored if there is a provided context in the method transform.
        if self.min-original-bounds ne self.min-transform-bounds || self.max-original-bounds ne self.max-transform-bounds {
            fail 'INCORRECT_ARGUMENTS: the original and transform arguments are expected to be the same.'
        }
    }

    method set-transform-axis(UInt:D $i) {
        @!indexes.push($i)
    }

    method drop-transform-axis(UInt:D $i) {
        @!indexes .= grep(* != $i)
    }

    method transform(
            :@point is copy,
            :$jacobian is copy,
            Bool:D :fb(:$functional-bounds) = False,
            :$context = Nil
            --> Map:D) {
        my %bounds = self.get-transform-bounds(:$context);

        if $functional-bounds {
            die 'Functional boundaries variable transformation is not implemented yet.'
        } else {
            my %res = :@point, :$jacobian;
            for @!indexes -> $i {
                %res<point>[$i] = %bounds<min>[$i] + %bounds<max>[$i] - %res<point>[$i];
                %res<jacobian> = %res<jacobian> * -1;
            }
            return %res
        }
    }
}
