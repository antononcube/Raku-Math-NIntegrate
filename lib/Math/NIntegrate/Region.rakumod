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
    has Int $.axis;

    # dimension == Length(@!A) == Length(@!B)
    has Int $.dimension;

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
    has Int $!no-error-decrease-count;

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

    method clone() {
        Math::NIntegrate::Region.new(
                min => @!min.clone,
                max => @!max.clone,
                fmin => @!fmin.clone,
                fmax => @!fmax.clone,
                abscissas => @!abscissas.clone,
                values => @!values.clone,
                integral => $!integral,
                error => $!error,
                levels => @!levels.clone,
                has-no-end-point => @!has-end-point.clone,
                no-error-decrease-count => $!no-error-decrease-count,
                type => $!type,
                reuse-values => @!reuse-values.clone,
                rule => $!rule.clone,
                nf => $!nf.clone,
                var-trans => $!var-trans,
                strategy => $!strategy
                )
    }

    #--------------------------------------
    # Public
    #--------------------------------------

    #| Get integration estimate
    method integrate(-->Math::NIntegrate::Region) {
        die 'Region.integrate is not implemented yet.'
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

    method split(Int:D $axis, Numeric:D $dithering = 0) {!!!}
    method partition(@nodex) {!!!}
    method duffy-transform() {!!!}
    method reverse-variable(Int:D $var-index) {!!!}
    method add-variable-transform(Math::NIntegrate::NumericalFunction:D $nf, Int:D $var-index) {!!!}
}