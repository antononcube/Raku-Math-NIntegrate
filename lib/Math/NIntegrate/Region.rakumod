use v6.d;

use Math::NIntegrate::Rule::General;
use Math::NIntegrate::NumericalFunction;
use Math::NIntegrate::VariableTransformer;
use Math::NIntegrate::Utilities;

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
        without self.dimension {
            $!dimension = @!min.elems
        }
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

        # The rule is not cloned
        $!rule = $from.rule;

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

    #| Get integration estimate
    method integrate(-->Math::NIntegrate::Region) {
        die 'Region.integrate is not implemented yet.'
    }

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

        # Which of these results is most useful?
        # return {left => self, right => $new-obj}
        # return (self, $new-obj)
        return $new-obj
    }

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
    # Future private methods
    #--------------------------------------

    method partition(@nodex) {!!!}
    method duffy-transform() {!!!}
    method reverse-variable(Int:D $var-index) {!!!}
    method add-variable-transform(Math::NIntegrate::NumericalFunction:D $nf, Int:D $var-index) {!!!}
}