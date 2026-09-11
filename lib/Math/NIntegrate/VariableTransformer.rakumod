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

    #| Scaling factor should be either 1 or -1
    has Numeric:D $.scale = 1;

    #| Contextual variable transformer (e.g. Composite)
    has Math::NIntegrate::VariableTransformer $.context;

    #| Integration region
    has $.region;

    #======================================================
    # Creators
    #======================================================

    #| Clone the object
    method clone(-->Math::NIntegrate::VariableTransformer:D) {
        Math::NIntegrate::VariableTransformer.new(
                transforms => @!transforms>>.clone.Array,
                jacobians => @!jacobians>>.clone.Array,
                min-original-bounds => @!min-original-bounds>>.clone.Array,
                max-original-bounds => @!max-original-bounds>>.clone.Array,
                min-transform-bounds => @!min-transform-bounds>>.clone.Array,
                max-transform-bounds => @!max-transform-bounds>>.clone.Array,
                :$!working-precision
                :$!scale,
                :$!context,
                :$!region
                )
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

    #| Abstract transform method
    method transform(:@points, :$jacobian, Bool:D :fb(:$functional-bounds) = False --> Map:D) {!!!}

    #======================================================
    # Jacobian
    #======================================================

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
            @min.push(%res<start>);
            $jacobian *= %res<length>
        }

        return %(:@length, :@min, :@middle, :$jacobian)
    }
}