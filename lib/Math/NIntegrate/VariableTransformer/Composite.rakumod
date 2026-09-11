use v6.d;

use Math::NIntegrate::VariableTransformer;

class Math::NIntegrate::VariableTransformer::Composite {
    #| Stack of Math::NIntegrate::VariableTransformer objects
    has Math::NIntegrate::VariableTransformer @.stack;

    #======================================================
    # Stack management methods
    #======================================================
    #| Add a transformer
    method add(Math::NIntegrate::VariableTransformer:D $obj) {
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

    method transform(:@points, :$jacobian, Bool:D :fb(:$functional-bounds) = False --> Map:D) {
        # Special treatment is needed for functional boundaries
        if $functional-bounds {
            die 'Functional boundaries variable transformation is not implemented yet.'
        } else {
            # Affine transformation is always done with Composite
        }

        return %(:@points, :$jacobian)
    }
}
