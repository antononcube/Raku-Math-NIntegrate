use v6.d;

use Math::NIntegrate::Rule::General;
use Math::NIntegrate::NumericalFunction;
use Math::NIntegrate::VariableTransformer;
use Math::NIntegrate::VariableTransformer::IMT;
use Math::NIntegrate::VariableTransformer::DoubleExponential;
use Math::NIntegrate::VariableTransformer::Infinite;
use Math::NIntegrate::Utilities;
use Math::NIntegrate::Codes;

class Math::NIntegrate::Region {

    #--------------------------------------
    # Data attributes
    #--------------------------------------

    # True if @!A or @!B contain functions
    # E.g. the range spec
    #   ((x,0,1), (y,0,1-x), (z,0,1-x-y))
    # is kept as
    #   [[0,0,0], [1, {1-$^a},{1-$^a-$^b}]]

    has Bool:D $.has-functional-bounds = False;
    has RangeEndCases @.range-end-cases;

    # Lists of interval/simplex bounds
    has @.min;
    has @.max;

    # Lists of values at the end points
    has @.fmin;
    has @.fmax;

    # Last abscissas
    has @.abscissas;

    # Last values
    has @.values;

    # Last axis to split on (for dimensions greater than 1)
    has UInt $.axis;

    # dimension == Length(@!A) == Length(@!B)
    has UInt $.dimension;

    # Maybe this should be just delegated to $!rule ?
    # See apply-rule()
    has $.integral = Whatever;
    has $.error = Whatever;

    # Levels of partitioning -- start from zero;
    #   @!levels[i] says how many times the original simplex
    #   has been partitioned along the axis i in order to reach the
    #   current values of $!vt.get-transforms-bounds
    has @.levels;

    # Boolean values whether one of the range points is an end point
    has @.has-end-point;

    # Accumulated number of levels of recursion for which the error
    # has failed to decrease by a factor of at least 7.
    # If this gets as high as 4 a message regarding convergence rate is issued.
    has UInt $.no-error-decrease-count;

    # Region type
    has Str $.type;

    # Reuse values
    has %.reuse-values;

    #--------------------------------------
    # Object attributes
    #--------------------------------------

    # Quadrature rule object
    has Math::NIntegrate::Rule::General $.rule is rw;

    # Numerical function object
    has Math::NIntegrate::NumericalFunction $.numerical-function is rw;

    # Variable transformator, usually a Math::NIntegrate::VariableTransformer::Composite object
    has Math::NIntegrate::VariableTransformer $.variable-transformer is rw;

    # Reference to the "main" integration object;
    # it should be Math::NIntegrate::Strategy or Whatever
    has $.strategy = Whatever;

    #--------------------------------------
    # Creators
    #--------------------------------------
    submethod TWEAK(*%args) {
        without $!dimension {
            $!dimension = @!min.elems
        }
        if @!levels.elems == 0 {
            @!levels = 0 xx $!dimension
        }
        without $!axis {
            $!axis = 0
        }
        @!range-end-cases = RE_BOTH xx $!dimension;
    }

    method copy(
            Math::NIntegrate::Region:D $from,
            Bool:D :deep(:deep-copy(:$clone)) = False
                ) {
        @!min = $clone ?? $from.min.clone !! $from.min;
        @!max = $clone ?? $from.max.clone !! $from.max;
        @!fmin = $clone ?? $from.fmin.clone !! $from.fmin;
        @!fmax = $clone ?? $from.fmax.clone !! $from.fmax;
        @!abscissas = $clone ?? $from.abscissas.clone !! $from.abscissas;
        @!values = $clone ?? $from.values.clone !! $from.values;
        $!axis = $from.axis;
        $!integral = $from.integral;
        $!error = $from.error;
        @!levels = $clone ?? $from.levels.clone !! $from.levels;
        @!has-end-point = $clone ?? $from.has-end-point.clone !! $from.has-end-point;
        $!no-error-decrease-count = $from.no-error-decrease-count;
        $!type = $from.type;
        %!reuse-values = $clone ?? $from.reuse-values.clone !! $from.reuse-values;
        @!range-end-cases = $clone ?? $from.range-end-cases.clone !! $from.range-end-cases;

        # This is somewhat redundant
        $!dimension = $from.dimension;

        # The rule must be cloned
        $!rule = $clone ?? $from.rule.clone !! $from.rule;

        # The numerical function is not cloned
        $!numerical-function = $from.numerical-function;

        # Note that the variable transformation object is cloned too,
        # but the region has to be replaced.
        $!variable-transformer = $clone ?? $from.variable-transformer.clone !! $from.variable-transformer;

        $!strategy = $from.strategy;

        return self
    }

    method clone() {
        my $obj = Math::NIntegrate::Region.new.copy(self, :clone);
        $obj.variable-transformer.region = $obj;
        return $obj
    }

    #--------------------------------------
    # Public
    #--------------------------------------

    #| Get working precision
    method get-working-precision() {
        without $!numerical-function {
            fail 'MISSING_OBJECT: no numerical function object for working precision request'
        }
        $!numerical-function ?? $!numerical-function.working-precision !! Whatever
    }

    #| Volume of the region
    method volume() {
        my %bounds = $!variable-transformer.get-transform-bounds();
        my @lengths = |Math::NIntegrate::Utilities::length-calc-md(%bounds<min>, %bounds<max>, :!mid-point, working-precision => $!variable-transformer.working-precision)<length>;
        my $volume = [*] |@lengths;
        return numerical($volume, $!variable-transformer.working-precision)
    }

    #| Is a region and its variable transformer en bloc computations ready.
    method is-en-bloc-ready(-->Bool:D) {

        # Constant ranges
        my $noFuncBounds = ! $!variable-transformer.has-functional-bounds;

        # Next check that the region has a Composite variable transformer with only
        # an Infinite transformer in its stack that has no concrete transformers Callable:D -- i.e. all ranges are finite.

        # This check cannot be used because of circular dependencies
        # $!variable-transformer ~~ Math::NIntegrate::VariableTransformer::Composite:D

        return
                ($!variable-transformer.^name ~~ / Composite / ) &&
                        ( $!variable-transformer.stack.elems == 0 ||
                                $!variable-transformer.stack.elems == 1 &&
                                        ($!variable-transformer.stack.head ~~ Math::NIntegrate::VariableTransformer::Infinite:D) &&
                                !($!variable-transformer.stack.head.transforms.any ~~ Callable:D)
                        )
                && $noFuncBounds
    }

    #--------------------------------------
    # Split
    #--------------------------------------

    #| Split region across given axis and dithering
    multi method split(Int:D $axis, Numeric:D $dithering = 0) {
        self.split(:$axis, :$dithering)
    }

    multi method split(Int:D :$axis!, Numeric:D :$dithering = 0) {
        die "The value of \$axis is expected to be an integer between 0 and {self.dimension}."
        unless 0 ≤ $axis ≤ self.dimension;

        # Make the right-side region of the split
        my $new-obj = self.clone;

        # Change the boundaries of the transformation object.
        # This is needed in order to map the integration rule abscissas to into the transformed half-ranges.
        without $!variable-transformer {
            fail 'MISSING_OBJECT: no variable transformer in region spliting.'
        }
        my %bounds = $!variable-transformer.get-transform-bounds();
        my $min = %bounds<min>[$axis] ~~ Callable:D ?? %bounds<min>[$axis](|%bounds<min>.head($axis)) !! %bounds<min>[$axis];
        my $max = %bounds<max>[$axis] ~~ Callable:D ?? %bounds<max>[$axis](|%bounds<max>.head($axis)) !! %bounds<max>[$axis];

        my $mid = ($min + $max) / 2 + $dithering * ($max - $min);

        $mid = numerical($mid, self.get-working-precision);
        $!variable-transformer.set-max-transform-bound($axis, $mid);
        $new-obj.variable-transformer.set-min-transform-bound($axis, $mid);

        # Level of splitting
        @!levels[$axis] += 1;
        $new-obj.levels[$axis] += 1;

        # Region mid-point to split over
        # This not an operational mid-point -- it is an "info".
        # See the variable transformer split below.
        # Special care is needed for functional boundaries.
        #my $mid-info = (@!min[$axis] + @!max[$axis]) / 2 + $dithering * (@!max[$axis] - @!min[$axis]);
        # Since Affine is always used 1/2 just have to transformed by $!variable-transformer
        #`[
        my $mid-info = (@!min[$axis] + @!max[$axis]) / 2 + $dithering * (@!max[$axis] - @!min[$axis]);
        # Should this precision setting be before or after the computation of the mid point?
        $mid-info = numerical($mid-info, self.get-working-precision);
        self.max[$axis] = $mid-info;
        $new-obj.min[$axis] = $mid-info;
        ]

        # Special treatment of reuse-values
        # TBD...
        # Full blown Rule class has to be implemented first.

        # Range end cases
        # These are used to decide should the singularity handlers IMT and DoubleExponent be applied or not
        given @!range-end-cases[$axis] {
            when RE_BOTH {
                @!range-end-cases[$axis] = RE_LEFT;
                $new-obj.range-end-cases[$axis] = RE_RIGHT
            }
            when RE_LEFT {
                @!range-end-cases[$axis] = RE_LEFT;
                $new-obj.range-end-cases[$axis] = RE_NONE
            }
            when RE_RIGHT {
                @!range-end-cases[$axis] = RE_NONE;
                $new-obj.range-end-cases[$axis] = RE_RIGHT
            }
            default {
                @!range-end-cases[$axis] = RE_NONE;
                $new-obj.range-end-cases[$axis] = RE_NONE
            }
        }

        # Which of these results is most useful?
        # return {left => self, right => $new-obj}
        # return (self, $new-obj)
        return $new-obj
    }

    #--------------------------------------
    # Divide
    #--------------------------------------

    #| Divide the region with a given array of number of divisions for each dimension.
    method divide(@divisions where @divisions.all ~~ UInt:D --> Array) {

        die 'The number of divisions is expected equal the region dimension.'
        unless @divisions.elems == self.dimension;

        die 'INCORRECT_ARGUMENTS: Each element of the divisions argument is expected to be an integer greater than 0.'
        unless @divisions.min ≥ 1;

        # Change the boundaries of the transformation object.
        without $!variable-transformer {
            fail 'MISSING_OBJECT: no variable transformer in region spliting.'
        }

        my $working-precision = self.variable-transformer.working-precision;

        my %bounds = $!variable-transformer.get-transform-bounds();

        # For each axis form partition pairs
        my @range-pairs = do for ^self.dimension -> $axis {
            my $min = %bounds<min>[$axis];
            my $max = %bounds<max>[$axis];
            my $h = numerical(numerical($max - $min, $working-precision) / @divisions[$axis], $working-precision);

            # Using ($min, $min + $h ... $max) can produce incomplete set boundaries
            # because of the precision manipulations
            my @bounds = $min;
            for ^@divisions[$axis] {
                @bounds.push( numerical(@bounds.tail + $h, $working-precision) )
            }

            @bounds.map({ numerical($_, $working-precision) }).rotor(2 => -1);
        }

        # Same code as in method partition.
        @range-pairs = @range-pairs.elems > 1 ?? |cross(|@range-pairs>>.Array) !! |@range-pairs.head.map({ [$_,] });

        # Cartesian product of the range pairs per axis
        my @regions = @range-pairs.map({
            my $obj = self.clone;
            $obj.levels = |self.levels.map(* + 1);
            # This transposes the min-max values, could be done beforehand as in Builder::make-region-ranges
            $_.map(*.head).kv.map(-> $axis, $b { $obj.variable-transformer.set-min-transform-bound($axis, $b) });
            $_.map(*.tail).kv.map(-> $axis, $b { $obj.variable-transformer.set-max-transform-bound($axis, $b) });

            $obj.range-end-cases = RE_NONE xx self.dimension;

            for ^self.dimension -> $axis {
                given @!range-end-cases[$axis] {
                    when $_ eq RE_BOTH {
                        if $obj.variable-transformer.get-transform-bounds<min>[$axis] == %bounds<min>[$axis] {
                            $obj.range-end-cases[$axis] = RE_LEFT
                        } elsif $obj.variable-transformer.get-transform-bounds<max>[$axis] == %bounds<max>[$axis] {
                            $obj.range-end-cases[$axis] = RE_RIGHT
                        }
                    }
                    when RE_LEFT {
                        if $obj.variable-transformer.get-transform-bounds<min>[$axis] == %bounds<min>[$axis] {
                            $obj.range-end-cases[$axis] = RE_LEFT
                        }
                    }
                    when RE_RIGHT {
                        if $obj.variable-transformer.get-transform-bounds<max>[$axis] == %bounds<max>[$axis] {
                            $obj.range-end-cases[$axis] = RE_RIGHT
                        }
                    }
                    default {
                        # It was already set above
                        $obj.range-end-cases[$axis] = RE_NONE
                    }
                }
            }

            $obj
        });

        return @regions
    }

    #--------------------------------------
    # Partition
    #--------------------------------------

    # The argument points can be a) user specified integral ranges points, or b) integration rule abscissas.
    # These points reside on different domains: a) on the original integral domain, and b) on the transformed domain.
    # For the latter, this method can facilitate the reuse of integral computations, by strategies, like, LocalAdaptive.
    # This means that partitioning over points (integration rule nodes) makes sense for 1D regions.
    # For nD regions (n > 1) the partitioning makes sense for Cartesian rules.
    # In order this method to "make sense" the reuse of integration values have to be implemented.

    #| Partition the region with a given array of partition points.
    method partition(@points --> Array) {

        # Change the boundaries of the transformation object.
        without $!variable-transformer {
            fail 'MISSING_OBJECT: no variable transformer in region spliting.'
        }

        # Should we check if the points are within region's ranges?
        my %bounds = $!variable-transformer.get-transform-bounds();

        # For each axis for partition pairs
        my @range-pairs = do for ^self.dimension -> $axis {
            my $min = %bounds<min>[$axis];
            my $max = %bounds<max>[$axis];
            [$min, |@points.map(*[$axis]), $max].squish(with => {($^x - $^y).abs ≤ 2 * $MACHINE_EPSILON}).rotor(2 => -1);
        }

        @range-pairs = @range-pairs.elems > 1 ?? cross(|@range-pairs) !! |@range-pairs.head.map({ [$_,] });

        # Cartesian product of the range pairs per axis
        my @regions = @range-pairs.map({
            my $obj = self.clone;
            # Should the level of the new regions be increased?
            # $obj.levels = |self.levels.map(* + 1);
            $_.map(*.head).kv.map(-> $axis, $b { $obj.variable-transformer.set-min-transform-bound($axis, $b) });
            $_.map(*.tail).kv.map(-> $axis, $b { $obj.variable-transformer.set-max-transform-bound($axis, $b) });
            $obj
        });

        return @regions
    }

    #--------------------------------------
    # Evaluate integrand
    #--------------------------------------

    #| Evaluate integrand over transformed arguments and multiply by the Jacobian
    method eval-integrand(@point is copy) {
        my $jacobian = 1;
        my $value;

        # Apply variable transformation
        with $!variable-transformer {
            my %res = $!variable-transformer.transform(:@point, :$jacobian);
            @point = |%res<point>;
            $jacobian = %res<jacobian>
        }

        # If the Jacobian is near zero return zero
        if is-zero($jacobian) {
            return 0
        }

        # Preparer derivative-signs
        my $derivative-signs = Whatever;
        if $!variable-transformer && $!variable-transformer.stack.tail ~~ Math::NIntegrate::VariableTransformer::IMT {
            $derivative-signs = $!variable-transformer.stack.tail.derivatives>>.sign;
        }

        # Evaluate the integrand functor over the transformed point
        try {
            $value = $!numerical-function.eval(@point, :$derivative-signs);
        }

        if $! || $value !~~ Numeric:D {
            fail 'NOT_A_NUMERICAL_FUNCTION: non-numerical integrand value is obtained.';
            #die 'Non-numerical integrand value is obtained.'
            return Whatever
        }

        # Should a check be made that $!variable-transformer.jacobian-factors is not empty?
        # If $!variable-transformer.jacobian-factors is empty we get 1, because [*] |() == 1
        return do if $!variable-transformer {
            my $factor = [*] |$!variable-transformer.jacobian-factors;
            $value * $jacobian * $factor
        } else {
            $value * $jacobian
        }
    }

    #--------------------------------------
    # Integrate
    #--------------------------------------

    #| Get integration estimate
    method apply-rule(-->Math::NIntegrate::Region) {
        $!rule.integrate(self);
        $!integral = $!rule.integral;
        $!error = $!rule.error;
        $!axis = $!rule.largest-error-axis;
        return self
    }

    #--------------------------------------
    # Add variable transformer
    #--------------------------------------
    method add-variable-transformer(Str:D $vtID, Int:D :$axis!, :$working-precision = Num) {

        # Cannot make this check because this would introduce a cyclic file dependency. Hence using class name matching.
        #if !($!variable-transformer ~~ Math::NIntegrate::VariableTransformer::Composite:D) {
        without $!variable-transformer.^name ~~ / Composite / {
            fail 'WRONG_TYPE: the region variable transformer is expected to be of type Math::NIntegrate::VariableTransformer::Composite.'
        }

        # These options are for any variable transformer would use.
        my @min-original-bounds = 0 xx $!dimension;
        my @max-original-bounds = 1 xx $!dimension;
        my @min-transform-bounds = 0 xx $!dimension;
        my @max-transform-bounds = 1 xx $!dimension;

        my %args =
                region => self,
                :$working-precision
                ;

        given $vtID.lc {
            when 'imt' {

                # Check if the last transformer in the Composite stack is IMT.
                # If not make a new object and add it.
                if $!variable-transformer.stack.tail ~~ Math::NIntegrate::VariableTransformer::IMT {

                    $!variable-transformer.stack.tail.set-transform-axis($axis);

                    my %bounds = self.variable-transformer.tail.get-transform-bounds;
                    $!variable-transformer.stack.tail.set-min-original-bound($axis, %bounds<min>[$axis]);
                    $!variable-transformer.stack.tail.set-min-original-bound($axis, %bounds<min>[$axis]);
                    $!variable-transformer.stack.tail.set-min-transform-bound($axis, 0);
                    $!variable-transformer.stack.tail.set-max-transform-bound($axis, 1);

                } else {
                    my $vt = Math::NIntegrate::VariableTransformer::IMT.new(
                            :@min-original-bounds,
                            :@max-original-bounds,
                            :@min-transform-bounds,
                            :@max-transform-bounds,
                            |%args);

                    $vt.set-transform-axis($axis);

                    $!variable-transformer.add($vt)
                }

            }

            when 'double-exponential' {

                # Applying this variable transformer requires DoubleExponentialRule.

                # Check if the region is one-dimensional
                fail 'DIMENSIONS_DO_NOT_MATCH: DoubleExponential variable transformer is only for 1D regions.'
                unless $!dimension == 1;

                # Check if the last transformer in the Composite stack is Infinite or the stack is empty.
                # DoubleExponential should replace the Infinite transformer.
                if $!variable-transformer.stack.elems == 0
                        || $!variable-transformer.stack.tail ~~ Math::NIntegrate::VariableTransformer::Infinite {

                    my %bounds = self.variable-transformer.tail.get-transform-bounds;

                    @min-original-bounds = |%bounds<min>;
                    @max-original-bounds = |%bounds<max>;

                    my $vt = Math::NIntegrate::VariableTransformer::DoubleExponential.new(
                            :@min-original-bounds,
                            :@max-original-bounds,
                            :@min-transform-bounds,
                            :@max-transform-bounds,
                            |%args);

                    $!variable-transformer.add($vt)
                }

            }

            when 'reverse' {
                # This have a corresponding implementation with Math::NIntegrate::VariableTransformer::Reverse.
                # That class is most of top-level use, not internally -- this internal handling is simple,
                # but since it is not very bureaucratic it is probably not that good.

                my %bounds = $!variable-transformer.get-transform-bounds();

                # Swap boundaries at $axis for the last variable transformer of the Composite stack
                my $tmp = $!variable-transformer.stack.tail.min-transform-bounds[$axis];
                $!variable-transformer.stack.tail.min-transform-bounds[$axis] = $!variable-transformer.stack.tail.max-transform-bounds[$axis];
                $!variable-transformer.stack.tail.max-transform-bounds[$axis] = $tmp;

                # Jacobian factor
                if $!variable-transformer.jacobian-factors.elems == 0 {
                    # This is somewhat too weak -- setting the factors here.
                    # But on the other hand it is the only place where it is set.
                    $!variable-transformer.jacobian-factors = 1 xx $!dimension
                }
                $!variable-transformer.jacobian-factors[$axis] *= -1;

                @!range-end-cases[$axis] = do given @!range-end-cases[$axis] {
                    when $_ eq RE_RIGHT { RE_LEFT }
                    when $_ eq RE_LEFT { RE_RIGHT }
                    default { $_ }
                }
            }
        }
    }

    #--------------------------------------
    # Future private methods
    #--------------------------------------

    method duffy-transform() {!!!}

}