use v6.d;

use Math::NIntegrate::Utilities;

class Math::NIntegrate::VariableTransformer {
    #| List of variable transformations for each dimension
    has @.transforms;

    #| List of Jacobian values for each dimension corresponding to the transformations
    has @.jacobians;

    #| Original lower boundary;
    #| for each variable before its transformation is applied.
    has @.min-original-bounds;

    #| Original upper boundary;
    #| for each variable before its transformation is applied.
    has @.max-original-bounds;

    #| Transformed lower boundary;
    #| for each variable after its transformation is applied.
    has @.min-transform-bounds;

    #| Transformed upper boundary;
    #| for each variable after its transformation is applied.
    has @.max-transform-bounds;

    #| Working precision
    has $.working-precision = Num;

    # Each jacobian factor should be either 1 or -1.
    # They simplify the reverse of variables, which is fundamental for singularity handling with IMT rule.
    #| Jacobian factors for each dimension
    has @.jacobian-factors;

    # Contextual variable transformer (e.g. Composite)
    # Cannot be used -- Raku hangs when being set in TWEAK or a Builder method.
    # has Math::NIntegrate::VariableTransformer $.context is rw;

    #| Integration region
    has $.region is rw;

    #======================================================
    # Creators
    #======================================================

#    submethod BUILD() {
#
#       # At some a dedicated Exception class have to be made
#       #fail 'REGION_NOT_SET' unless self.region;
#
#       #self.scale = self.region.dimension xx 1;
#    }

    #| Copy attributes from an object
    method copy(Math::NIntegrate::VariableTransformer:D $from, Bool:D :deep(:deep-copy(:$clone)) = False) {
        # @!transforms are functions/subs, so, $from.transforms>>.clone.Array seems inappropriate.
        @!transforms = $clone ?? $from.transforms.clone !! $from.transforms;
        @!jacobians = $clone ?? $from.jacobians.clone !! $from.jacobians;
        @!min-original-bounds = $clone ?? $from.min-original-bounds>>.clone.Array !! $from.min-original-bounds;
        @!max-original-bounds = $clone ?? $from.max-original-bounds>>.clone.Array !! $from.max-original-bounds;
        @!min-transform-bounds = $clone ?? $from.min-transform-bounds>>.clone.Array !! $from.min-transform-bounds;
        @!max-transform-bounds = $clone ?? $from.max-transform-bounds>>.clone.Array !! $from.max-transform-bounds;
        $!working-precision = $from.working-precision;
        @!jacobian-factors = $clone ?? $from.jacobian-factors>>.clone.Array !! $from.jacobian-factors;
        $!region = $from.region;
            
        return self
    }

    #| Clone the object
    method clone(-->Math::NIntegrate::VariableTransformer:D) {
        Math::NIntegrate::VariableTransformer.new.copy(self, :clone)
    }

    #======================================================
    # Template Method methods
    #======================================================

    #| Get region object
    method get-region(:$context = Nil) {
        return do if $context {
            $context.get-region()
        } else {
            self.region
        }
    }

    #| Get original bounds
    method get-original-bounds(:$context = Nil -->Map:D) {
        return $context ?? $context.get-original-bounds !! %(min => @!min-original-bounds, max => @!max-original-bounds)
    }

    #| Get transformation bounds
    method get-transform-bounds(:$context = Nil -->Map:D) {
        return $context ?? $context.get-transform-bounds() !! %(min => @!min-transform-bounds, max => @!max-transform-bounds)
    }

    #| Set min original bound(s)
    method set-min-original-bound(UInt:D $axis, Numeric:D $value) {
        @!min-original-bounds[$axis] = $value;
        return self
    }

    #| Set max original bound(s)
    method set-max-original-bound(UInt:D $axis, Numeric:D $value) {
        @!max-original-bounds[$axis] = $value;
        return self
    }

    #| Set min transformation bound(s)
    method set-min-transform-bound(UInt:D $axis, Numeric:D $value) {
        @!min-transform-bounds[$axis] = $value;
        return self
    }

    #| Set max transformation bound(s)
    method set-max-transform-bound(UInt:D $axis, Numeric:D $value) {
        @!max-transform-bounds[$axis] = $value;
        return self
    }

    #| Functional boundaries determination
    method has-functional-bounds(:$context = Nil -->Bool:D) {
        my %bounds = self.get-transform-bounds(:$context);
        my %orig-bounds = self.get-original-bounds(:$context);
        return %bounds.values.flat(:hammer).any ~~ Callable:D || %orig-bounds.values.flat(:hammer).any ~~ Callable:D;
    }

    #| Abstract transform method
    method transform(
            :@point is copy,
            :$jacobian is copy,
            Bool:D :fb(:$functional-bounds) = False,
            :$context = Nil
            --> Map:D) {!!!}
}