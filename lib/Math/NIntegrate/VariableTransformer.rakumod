use v6.d;

use Math::NIntegrate::Utilities;

class Math::NIntegrate::VariableTransformer {
    #| List of variable transformations for each dimension
    has @.transforms;

    #| List of jacobians for each dimension corresponding to the transformations
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

    #======================================================
    # Methods
    #======================================================
    method clone(-->Math::NIntegrate::VariableTransformer:D) {
        Math::NIntegrate::VariableTransformer.new(
                transforms => @!transforms>>.clone.Array,
                jacobians => @!jacobians>>.clone.Array,
                min-original-bounds => @!min-original-bounds>>.clone.Array,
                max-original-bounds => @!max-original-bounds>>.clone.Array,
                min-transform-bounds => @!min-transform-bounds>>.clone.Array,
                max-transform-bounds => @!max-transform-bounds>>.clone.Array,
                :$!working-precision
                :$!scale
                )
    }

    method !length-calc(Numeric:D $a is copy, Numeric:D $b is copy, Bool:D $mid-point = False) {

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

        return %(:$length, start => $a, :$middle, jacobian => $length)
    }

    method !length-calc-multidimensional(@a, @b, Bool:D $mid-point = False) {
        die 'The sizes of the first two arguments are expected to match' unless @a.elems == @b.elems;

        my $jacobian = numerical(1, self.working-precision);
        my @length;
        my @middle;
        my @start;
        for ^@a.elems -> $i {
            my %res = self!length-calc(@a[$i], @b[$i], $mid-point);
            @length.push(%res<length>);
            @middle.push(%res<middle>);
            @start.push(%res<start>);
            $jacobian *= %res<length>
        }

        return %(:@length, :@start, :@middle, :$jacobian)
    }
}