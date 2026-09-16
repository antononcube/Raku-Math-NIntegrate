use v6.d;

use Math::NIntegrate::VariableTransformer;
use Math::NIntegrate::VariableTransformer::Affine;
use Math::NIntegrate::VariableTransformer::Infinity;
use Math::NIntegrate::Region;

class Math::NIntegrate::VariableTransformer::Composite
        is Math::NIntegrate::VariableTransformer {
    #| Stack of Math::NIntegrate::VariableTransformer objects
    has Math::NIntegrate::VariableTransformer @.stack;
    has Math::NIntegrate::VariableTransformer::Affine $.vtAffine;

    #======================================================
    # Creators
    #======================================================
    submethod TWEAK(*%args) {
        # Composite always uses the Affine variable transformer
        # in order to put the rule abscissas within the boundaries of region.
        $!vtAffine = Math::NIntegrate::VariableTransformer::Affine.new(context => self);

        # Should the stack have at least one variable transformer?
        # The Infinity transform can be seen as finite-range proxy of the region.
        # The Affine transform is for mapping or the finite boundaries to integration rules [0, 1] abscissas ranges.

        # All variable transformers in the stack have this object as context.
        # The method get-region() is going return this object's region.
        for @!stack -> $vt { $vt.context = self }

        # Should Composite be made to always have a region object?
        # without self.region {
        #     fail 'MISSING_OBJECT: composite initialized without region object'
        # }
    }

#    method new(Math::NIntegrate::Region:D $region) {
#
#        #        if @!stack.elems == 0 {
#        #            my $vtInf = Math::NIntegrate::VariableTransformer::Infinity.new(context => self, region => self.region);
#        #            say (:$vtInf);
#        #            @!stack.push($vtInf)
#        #        }
#    }

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

    #| Copy from object
    method copy(Math::NIntegrate::VariableTransformer::Composite:D $from,  Bool:D :deep(:deep-copy(:$clone)) = False) {
        # Delegate to parent class
        self.Math::NIntegrate::VariableTransformer::copy($from, :$clone);

        # No need to copy the Affine object
        @!stack = $clone ?? $from.stack>>.clone !! $from.stack;

        return self
    }

    method clone(-->Math::NIntegrate::VariableTransformer::Composite) {

        # No need to clone the Affine object
        Math::NIntegrate::VariableTransformer::Composite.new().copy(self)
    }

    method get-original-bounds(-->Map:D) {
        return @!stack.tail.get-original-bounds
    }

    method get-transform-bounds(-->Map:D) {
        return @!stack.elems
                ?? @!stack.tail.get-transform-bounds
                !! self.Math::NIntegrate::VariableTransformer::get-transform-bounds
    }

    method get-scale(-->Numeric:D) {
        return reduce({$^a * $^b.scale}, self.scale, |@!stack>>.scale )
    }

    method transform(:@point is copy, :$jacobian is copy, Bool:D :fb(:$functional-bounds) = False --> Map:D) {

        $jacobian = 1;
        # Special treatment is needed for functional boundaries
        if $functional-bounds {
            die 'Functional boundaries variable transformation is not implemented yet.'
        } else {
            # Affine transformation is always done with Composite
            if $.vtAffine {
                # Affine transform is multidimensional
                my %res = $.vtAffine.transform(:@point, :$jacobian);
                @point = |%res<point>;
                $jacobian = %res<jacobian>
            }

            # Apply the stack of transformations in reverse order
            for @!stack.reverse -> $vt {
                my %res = $vt.transform(:@point, :$jacobian);
                @point = |%res<point>;
                $jacobian = %res<jacobian>
            }
        }

        return %(:@point, :$jacobian)
    }
}
