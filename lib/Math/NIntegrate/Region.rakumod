use v6.d;

use Math::NIntegrate::Rule::General;
use Math::NIntegrate::NumericalFunction;
use Math::NIntegrate::VariableTransformer;

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
    has UInt $!no-error-decrease-count;

    # Region type
    has Str $.type;

    # Reuse values
    has @.reuse-values;

    #--------------------------------------
    # Object attributes
    #--------------------------------------

    # Quadrature rule object
    has Math::NIntegrate::Rule::General $.rule;

    # Numerical function object
    has Math::NIntegrate::NumericalFunction $.nf;

    # Variable transformator, usually a Math::NIntegrate::VariableTransformer::Composite object
    has Math::NIntegrate::VariableTransformer $.var-trans;

    # Reference to the "main" integration object;
    # it should be Math::NIntegrate::Strategy or Whatever
    has $.strategy = Whatever;

    #--------------------------------------
    # Creators
    #--------------------------------------

    #method new() {!!!}

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
        $!rule = $from.rule;
        $!nf = $from.nf;
        $!var-trans = $clone ?? $from.var-trans.clone !! $from.var-trans;
        $!strategy = $from.strategy;

        return self
    }

    method clone() {
        Math::NIntegrate::Region.new.copy(self, :clone)
    }

    #--------------------------------------
    # Public
    #--------------------------------------

    #| Get integration estimate
    method integrate(-->Math::NIntegrate::Region) {
        die 'Region.integrate is not implemented yet.'
    }

    method split(Int:D :$axis!, Numeric:D :$dithering = 0) {
        die "The value of \$axis is expected to be an integer between 0 and {self.dimension}."
        unless 0 ≤ $axis ≤ self.dimension;

        my $new-obj = self.clone;


    }

    #--------------------------------------
    # Future private methods
    #--------------------------------------
    #| Evaluate object's integrand.
    method eval-integrand(@args) {
        # Apply variable transformation
        # Make sure get the
        die 'Region.eval-integrand is not implemented yet.'
    }

    method partition(@nodex) {!!!}
    method duffy-transform() {!!!}
    method reverse-variable(Int:D $var-index) {!!!}
    method add-variable-transform(Math::NIntegrate::NumericalFunction:D $nf, Int:D $var-index) {!!!}
}