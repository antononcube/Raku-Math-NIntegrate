use v6.d;

use Math::NIntegrate::Utilities;
use Math::NIntegrate::NumericalFunction;
use Hash::Merge;

class Math::NIntegrate::Spec {
    has &.integrand;
    has %.ranges;
    has %.options =
            method => Whatever,
            min-recursion => 1,
            max-recursion => 12,
            max-points => Whatever,
            working-precision => Num,
            precision-goal => Whatever,
            accuracy-goal => Whatever,
            evaluation-monitor => Nil
            ;


    #======================================================
    # Creators
    #======================================================
    submethod BUILD(:&!integrand, :%!ranges, :%!options) {
        self.normalize(&!integrand, %!ranges);
        self.normalize-options(%!options)
    }

    multi method new(&integrand, %ranges, %options) {
        self.bless(:&integrand, :%ranges, :%options)
    }

    multi method new(:&integrand, :%ranges, :%options) {
        self.bless(:&integrand, :%ranges, :%options)
    }

    #======================================================
    # Normalization methods
    #======================================================

    #| Validate & normalize ranges
    proto method normalize-ranges($ranges) {*}

    multi method normalize-ranges(@ranges) {
        my $msg = 'A ranges list is expected to be a list of lists, each list having a string (variable name) as a first element';

        die $msg unless @ranges.all ~~ (Array:D | List:D | Seq:D);

        die $msg unless @ranges>>.head.all ~~ Str:D && ([&&] @ranges>>.elems <<≥>> 3);

        return @ranges.kv.map( -> $i, @r { @r.head => %(index => $i, var => @r.head, a => @r[1], b => @r.tail, interval-points => @r[2... @r.elems - 2]) }).Hash
    }

    multi method normalize-ranges(%ranges is copy) {
        my $msg1 = 'A ranges hashmap is expected to have string varible names as keys and hashmaps as values.';
        my $msg2 = 'Each value of a ranges hashmap is expected to have the keys: "a", "b", and "index".';

        die $msg1 unless %ranges.keys.all ~~ Str:D;
        die $msg1 unless %ranges.values.all ~~ Map:D;

        # At some point no-index should be allowed and index => Whatever be handled.
        # The indexes can be concluded from integrand's arguments names and positions.
        die $msg2 unless [&&] %ranges.map({ ($_.value.keys (&) <a b index>).elems == 3});

        %ranges .= map({ $_.key => merge-hash($_.value, %(var => $_.head) ) });

        # Check range indexes are integers between 0 and %ranges.elems and unique
        my @indexes = %ranges.values.map({ $_.value<index> }).unique;
        die 'All ranges indexes are expected to be integers.' unless @indexes.all ~~ Int:D;
        die 'Range indexes are expected to be unique.' unless @indexes.elems < %ranges.elems;
        die 'Range indexes are expected to be between 0 and the ranges spec length.' unless @indexes.min == 0 && @indexes.max == %ranges.elems - 1;

        return %ranges
    }

    #| Verify options
    method normalize-options(%options is copy) {

        %options = merge-hash(%!options, %options);

        # Working precision
        %options<working-precision> = Num if %options<working-precision>.isa(Whatever);

        die 'The value of working precision is expected to be Num, Rat, FatRat, or Whatever.'
        unless %options<working-precision> ~~ (Num | Rat | FatRat);

        # Precision goal
        if %options<precision-goal> ~~ Numeric:D {
            die 'The values of precision-goal is expected to be a positive number or Whatever.'
            unless %options<precision-goal> > 0;
        } else {
            %options<precision-goal> = do given %options<working-precision> {
                when Num { 6 }
                when FatRat { 20 }
                when Rat { 6 }
                default {
                    die 'Cannot process the value of precision-goal.'
                }
            }
        }

        # Accuracy
        %options<accuracy-goal> = Inf if %options<accuracy-goal>.isa(Whatever);

        die 'The value of accuracy-goal is expected to be a positive integer, Inf, or Whatever.'
        unless %options<accuracy-goal> ~~ Int:D && %options<accuracy-goal> > 0 || %options<accuracy-goal> ~~ Inf;

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

        # Assign
        %!options = %options;
    }

    #| Normalize integrand and ranges
    proto method normalize(&func, $ranges) {*}
    multi method normalize(&func, @ranges) {
        self.normalize(&func, self.normalize-ranges(@ranges))
    }
    multi method normalize(&func, %ranges, %options) {
        # The arity of the function should equal %ranges.elems.
        %!ranges = self.normalize-ranges(%ranges);

        # Check argument names correspond to variable names in ranges
        &!integrand = Math::NIntegrate::NumbericalFunction.new(function => &func);

        %!options = self.normalize-options(%options);
    }
}
