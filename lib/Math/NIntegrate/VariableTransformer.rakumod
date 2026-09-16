use v6.d;

use Math::NIntegrate::Utilities;

class Math::NIntegrate::VariableTransformer {
    #| List of variable transformations for each dimension
    has @.transforms;

    #| List of Jacobians for each dimension corresponding to the transformations
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

    #| Scaling factors for each dimension; each scale should be either 1 or -1
    has @.scales;

    #| Contextual variable transformer (e.g. Composite)
    has Math::NIntegrate::VariableTransformer $.context is rw;

    #| Integration region
    has $.region;

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
        @!scales = $clone ?? $from.scales>>.clone.Array !! $from.scales;
        $!context = $from.context;
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

    #| Get original bounds
    method get-original-bounds(-->Map:D) {
        return do if $!context {
            $!context.get-original-bounds()
        } else {
            %(min => @!min-original-bounds, max => @!max-original-bounds)
        }
    }

    #| Get transformation bounds
    method get-transform-bounds(-->Map:D) {
        return do if $!context {
            $!context.get-transform-bounds()
        } else {
            %(min => @!min-transform-bounds, max => @!max-transform-bounds)
        }
    }

    #| Get region object
    method get-region(-->Map:D) {
        return do if $!context {
            $!context.get-region()
        } else {
            self.region
        }
    }

    #| Abstract transform method
    method transform(:@point, :$jacobian, Bool:D :fb(:$functional-bounds) = False --> Map:D) {!!!}

    #======================================================
    # Jacobian
    #======================================================

    # This method probably should be only called by length-calc-md --
    # the framework should work only with points that are lists.
    # (After initial processing.)
    #| Calculation of the interval length with working precision
    method length-calc(Numeric:D $a is copy, Numeric:D $b is copy, Bool:D $mid-point = False -->Map:D) {

        # Numeric evaluation of $a
        $a = numerical($a, self.working-precision);

        # Numeric evaluation of $b
        $b = numerical($b, self.working-precision);

        # The length of the interval
        my $length = $b - $a;

        # Mid-point calculation if needed
        my $middle = do if $mid-point {
            numerical(($a + $b) / 2, self.working-precision)
        } else {
            numerical( 1 / 2, self.working-precision)
        }

        return %(:$length, min => $a, :$middle, jacobian => $length)
    }

    #| Calculation of the multidimensional interval lengths with working precision
    method length-calc-md(@a, @b, Bool:D $mid-point = False -->Map:D) {
        die 'The sizes of the first two arguments are expected to match' unless @a.elems == @b.elems;

        my $jacobian = numerical(1, self.working-precision);
        my @length;
        my @middle;
        my @min;
        for ^@a.elems -> $i {
            my %res = self.length-calc(@a[$i], @b[$i], $mid-point);
            @length.push(%res<length>);
            @middle.push(%res<middle>);
            @min.push(%res<min>);
            $jacobian *= %res<length>
        }

        return %(:@length, :@min, :@middle, :$jacobian)
    }
}