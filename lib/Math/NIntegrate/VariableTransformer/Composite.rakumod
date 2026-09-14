use v6.d;

use Math::NIntegrate::VariableTransformer;
use Math::NIntegrate::VariableTransformer::Affine;

class Math::NIntegrate::VariableTransformer::Composite {
    #| Stack of Math::NIntegrate::VariableTransformer objects
    has Math::NIntegrate::VariableTransformer @.stack;
    has Math::NIntegrate::VariableTransformer::Affine $.vtAffine;

    #======================================================
    # Creators
    #======================================================
    submethod TWEAK (*%args) {
        # Composite always uses the Affine variable transformer
        # in order to put the rule abscissas within the boundaries of region.
        $!vtAffine = Math::NIntegrate::VariableTransformer::Affine.new(context => self);

        # All variable transformers in the stack have this object as context.
        # The method get-region() is going return this object's region.
        with @!stack {
            for @!stack -> $vt { $vt.context = self }
        }

        # Should Composite be made to always have a region object?
        # without self.region {
        #     fail 'MISSING_OBJECT: composite initialized without region object'
        # }
    }

    #======================================================
    # Stack management methods
    #======================================================
    #| Add a transformer
    method add(Math::NIntegrate::VariableTransformer:D $obj) {

        if $obj.region.dimension != self.region.dimension {
            fail 'DIMENSIONS_DO_NOT_MATCH: Non-equal dimension when adding a variable transformer.'
        }

        @!stack.push($obj);
        return self
    }

    #| Remove a transformer object
    multi method remove(Math::NIntegrate::VariableTransformer:D $obj) {
        @!stack .= grep({ $_ ne $obj });
        return self
    }

    #| Remove a transformer at give index
    multi method remove(Int:D $index) {
        @!stack = @!stack.splice($index, 1);
        return self
    }

    #======================================================
    # Composing methods
    #======================================================
    method clone(-->Math::NIntegrate::VariableTransformer::Composite) {
        Math::NIntegrate::VariableTransformer::Composite.new(stack => @!stack>>.clone)
    }

    method get-original-bounds(-->Map:D) {
        return @!stack.tail.get-original-bounds
    }

    method get-transform-bounds(-->Map:D) {
        return @!stack.tail.get-tranform-bounds
    }

    method get-scale(-->Numeric:D) {
        return reduce({$^a * $^b.scale}, self.scale, |@!stack>>.scale )
    }

    method transform(:@point, :$jacobian, Bool:D :fb(:$functional-bounds) = False --> Map:D) {
        # Special treatment is needed for functional boundaries
        if $functional-bounds {
            die 'Functional boundaries variable transformation is not implemented yet.'
        } else {
            # Affine transformation is always done with Composite
        }

        return %(:@point, :$jacobian)
    }
}
