use v6.d;

use Math::NIntegrate::Rule::General;
use Math::NIntegrate::NumericalFunction;
use Math::NIntegrate::VariableTransformer;
use Math::NIntegrate::VariableTransformer::IMT;
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
    has @.reuse-values;

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
        $!integral = $from.integral;
        $!error = $from.error;
        @!levels = $clone ?? $from.levels.clone !! $from.levels;
        @!has-end-point = $clone ?? $from.has-end-point.clone !! $from.has-end-point;
        $!no-error-decrease-count = $from.no-error-decrease-count;
        $!type = $from.type;
        @!reuse-values = $clone ?? $from.reuse-values.clone !! $from.reuse-values;

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

        # Region mid-point to split over
        # This not an operational mid-point -- it is an "info".
        # See the variable transformer split below.
        my $mid = (@!min[$axis] + @!max[$axis]) / 2 + $dithering * (@!max[$axis] - @!min[$axis]);
        # Should this precision setting be before or after the computation of the mid point?
        $mid = numerical($mid, self.get-working-precision);
        self.max[$axis] = $mid;
        $new-obj.min[$axis] = $mid;

        # Change the boundaries of the transformation object.
        # This is needed in order to map the integration rule abscissas to into the transformed half-ranges.
        without $!variable-transformer {
            fail 'MISSING_OBJECT: no variable transformer in region spliting.'
        }
        my %bounds = $!variable-transformer.get-transform-bounds();
        my $min = %bounds<min>[$axis];
        my $max = %bounds<max>[$axis];

        $mid = ($min + $max) / 2 + $dithering * ($max - $min);

        $mid = numerical($mid, self.get-working-precision);
        $!variable-transformer.set-max-transform-bound($axis, $mid);
        $new-obj.variable-transformer.set-min-transform-bound($axis, $mid);

        # Level of splitting
        @!levels[$axis] += 1;
        $new-obj.levels[$axis] += 1;

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

    #| Divide the region with a given array of number of divisions for each dimension
    method divide(@divisions --> Array) {
        die 'The region method divide is not implemented yet.'
    }

    #--------------------------------------
    # Partition
    #--------------------------------------

    # The argument points are, likely, integration rule abscissas.
    # This method facilitates reuse of integral computations, by strategies, like, LocalAdaptive.

    #| Partition the region with a given array of partition points.
    method parition(@points --> Array) {
        die 'The region method paritition is not implemented yet.'
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

        # Evaluate the integrand functor over the transformed point
        try {
            $value = $!numerical-function.eval(@point);
        }

        if $! || $value !~~ Numeric:D {
            fail 'NOT_A_NUMERICAL_FUNCTION: non-numerical integrand value is obtained'
        }

        # Should a check be made that $!variable-transformer.scales is not empty?
        # If $!variable-transformer.scales is empty we get 1, because [*] |() == 1
        return do if $!variable-transformer {
            my $scale = [*] |$!variable-transformer.scales;
            $value * $jacobian * $scale
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
        return self
    }

    #--------------------------------------
    # Add variable transformer
    #--------------------------------------
    method add-variable-transformer(Str:D $vtID, Int:D :$axis!, :$working-precision = Num) {

        # These are options are for any variable transformer would use.
        my @min-original-bounds = 0;
        my @max-original-bounds = 1;
        my @min-transform-bounds = 0 xx $!dimension;
        my @max-transform-bounds = 1 xx $!dimension;

        my %args =
                region => self,
                working-precision => Rat
                ;

        given $vtID.lc {
            when 'imt' {

                # Check if the last transformer in the Composite stack is IMT.
                # If not make a new object and add it.
                if $!variable-transformer.stack.tail ~~ Math::NIntegrate::VariableTransformer::IMT {

                    $!variable-transformer.stack.tail.set-transform-axis($axis)

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
        }
    }

    #--------------------------------------
    # Future private methods
    #--------------------------------------

    method partition(@nodex) {!!!}
    method duffy-transform() {!!!}
    method reverse-variable(Int:D $var-index) {!!!}

}