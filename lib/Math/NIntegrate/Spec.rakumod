use v6.d;

use Math::NIntegrate::Utilities;
use Math::NIntegrate::NumericalFunction;
use Math::NIntegrate::Processing::Grammar;
use Math::NIntegrate::Processing::Actions::MethodSpec;
use Hash::Merge;

class Math::NIntegrate::Spec {
    has $.integrand;
    has %.ranges;
    has %.options;
    my %default-options =
            method => Whatever,
            min-recursion => 0,
            max-recursion => 12,
            max-points => Whatever,
            working-precision => Num,
            precision-goal => Whatever,
            accuracy-goal => Whatever,
            evaluation-monitor => Nil,
            singularity-depth => 4
            ;


    #======================================================
    # Creators
    #======================================================
    submethod BUILD(:&integrand, :$ranges, :$options = Whatever) {
        self.normalize(&integrand, $ranges, $options);
    }

    multi method new(&integrand, $ranges, $options = Whatever) {
        self.bless(:&integrand, :$ranges, :$options)
    }

    multi method new(:&integrand, :$ranges, :$options = Whatever) {
        self.bless(:&integrand, :$ranges, :$options)
    }

    #======================================================
    # Normalization methods
    #======================================================

    #------------------------------------------------------
    #| Validate & normalize ranges
    proto method normalize-ranges($ranges) {*}

    multi method normalize-ranges(@ranges) {
        my $msg = 'A ranges list is expected to be a list of lists, each list having a string (variable name) as a first element';

        die $msg unless @ranges.all ~~ (Array:D | List:D | Seq:D);

        die $msg unless @ranges>>.head.all ~~ Str:D && ([&&] @ranges>>.elems <<≥>> 3);

        return @ranges.kv.map( -> $i, @r { @r.head => %(index => $i, var => @r.head, min => @r[1], max => @r.tail, interval-points => @r[2... @r.elems - 2]) }).Hash
    }

    multi method normalize-ranges(%ranges is copy) {
        my $msg1 = 'A ranges hashmap is expected to have string varible names as keys and hashmaps as values.';
        my $msg2 = 'Each value of a ranges hashmap is expected to have the keys: "min", "max", and "index".';

        die $msg1 unless %ranges.keys.all ~~ Str:D;
        die $msg1 unless %ranges.values.all ~~ Map:D;

        # At some point no-index should be allowed and index => Whatever be handled.
        # The indexes can be concluded from integrand's arguments names and positions.
        die $msg2 unless [&&] %ranges.map({ ($_.value.keys (&) <max min index>).elems == 3});

        %ranges .= map({ $_.key => merge-hash($_.value, %(var => $_.head) ) });

        # Check range indexes are integers between 0 and %ranges.elems and unique
        my @indexes = %ranges.values.map({ $_<index> }).unique;
        die 'All ranges indexes are expected to be integers.' unless @indexes.all ~~ Int:D;
        die 'Range indexes are expected to be unique.' if @indexes.elems < %ranges.elems;
        die 'Range indexes are expected to be between 0 and the ranges spec length.' unless @indexes.min == 0 && @indexes.max == %ranges.elems - 1;

        return %ranges
    }

    #------------------------------------------------------
    #| Normalize method spec
    method normalize-method($method is copy --> Map:D) {
        # The normalizations of the non-Whatever specs are not needed because
        # they are handled by Math::NIntegrate::Processing::Grammar .
        # Some other normalizations might be implemented later.

        # This is the normalized variant of $method and has to be supported at some point.
        # %( name => 'GlobalAdaptive', method => %( name => 'ClenshawCurtisRule', points => 6 ) )

        $method = do given $method {
            when $_.isa(WhateverCode) || $_.isa(Whatever) { 'GlobalAdaptive' }
            when $_ ~~ Str:D && $_.lc ∈ <automatic auto> { 'GlobalAdaptive' }
            default { $method }
        }

        my $gr = Math::NIntegrate::Processing::Grammar.new;
        my $actions = Math::NIntegrate::Processing::Actions::MethodSpec.new;

        my %method-spec;
        try {
            %method-spec = $gr.parse($method.List.raku, rule => 'TOP', :$actions).made
        }

        if $! {
            die 'Cannot parse method option.'
        }

        return %method-spec
    }

    #------------------------------------------------------
    #| Verify options
    method normalize-options(%options is copy --> Map:D) {

        %options = merge-hash(%default-options, %options);

        # Working precision
        %options<working-precision> = Num without %options<working-precision>;

        die 'The value of working precision is expected to be Num, Rat, FatRat, or Whatever.'
        unless %options<working-precision> ~~ (Num | Rat | FatRat);

        # Precision goal
        my $msg-pg = 'The value of precision-goal is expected to be a positive number or Whatever.';
        if %options<precision-goal> ~~ Numeric:D {
            die $msg-pg unless %options<precision-goal> > 0;
        } elsif %options<precision-goal>.isa(Whatever) {
            %options<precision-goal> = do given %options<working-precision> {
                when Num { 6 }
                when FatRat { 20 }
                when Rat { 6 }
                default {
                    die 'Cannot deduce the value of precision-goal.'
                }
            }
        } else {
            die $msg-pg
        }

        %options<relative-tolerance> = 10 ** -%options<precision-goal>;

        # Accuracy
        %options<accuracy-goal> = Inf if %options<accuracy-goal>.isa(Whatever);

        die 'The value of accuracy-goal is expected to be a positive integer, Inf, or Whatever.'
        unless %options<accuracy-goal> ~~ Int:D && %options<accuracy-goal> > 0 || %options<accuracy-goal> ~~ Inf;

        %options<absolute-tolerance> = 10 ** -%options<accuracy-goal>;

        # Recursion options
        die 'The value of max-recursion is expected to be a non-negative integer.'
        unless %options<max-recursion> ~~ Int:D && %options<max-recursion> ≥ 0;

        die 'The value of min-recursion is expected to be a non-negative integer.'
        unless %options<min-recursion> ~~ Int:D && %options<min-recursion> ≥ 0;

        die 'The value of min-recursion is expected to be less or equal to the value of max-recursion.'
        unless %options<min-recursion> ≤ %options<max-recursion>;

        # Max points
        die 'The value of max-points is expected to be a positive integer or Whatever.'
        unless %options<max-points> ~~ Int:D && %options<max-points> > 0 || %options<max-points>.isa(Whatever);

        # Singularity depth
        %options<singularity-depth> = 4 if %options<singularity-depth>.isa(Whatever);

        die 'The value of singularity-depth is expected to be a non-negative integer, Inf, or Whatever.'
        unless %options<singularity-depth> ~~ Int:D && %options<singularity-depth> ≥ 0 || %options<singularity-depth> ~~ Inf;

        # Method option
        %options<method> = self.normalize-method(%options<method> // Whatever);

        # Assign
        %!options = %options;
    }


    #------------------------------------------------------
    #| Normalize integrand and ranges
    proto method normalize(&func, $ranges, $options) {*}

    multi method normalize(&func, @ranges, $options) {
        self.normalize(&func, self.normalize-ranges(@ranges), $options)
    }

    multi method normalize(&func, %ranges, $options is copy) {

        # The arity of the function should equal %ranges.elems.
        %!ranges = self.normalize-ranges(%ranges);

        # Options
        $options = %default-options if $options.isa(Whatever);

        die 'The value of options is expected to be a hashmap or Whatever.' unless $options ~~ Map:D;
        %!options = self.normalize-options($options);

        # Numerical function object for integrand
        $!integrand = Math::NIntegrate::NumericalFunction.new(function => &func, working-precision => %!options<working-precision>);

        # Check argument names correspond to variable names in ranges
        die 'None of the integrand arguments are found in the ranges spec.'
        unless ($!integrand.argument-names (&) %ranges.keys).elems == 0;

        die 'Not all integrand arguments are found in the ranges spec.'
        unless ($!integrand.argument-names (&) %ranges.keys).elems < $!integrand.argument-names;

        self
    }
}
