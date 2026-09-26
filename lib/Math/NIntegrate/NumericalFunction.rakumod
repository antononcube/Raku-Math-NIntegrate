use v6.d;

class Math::NIntegrate::NumericalFunction {
    has &.function is rw;
    has @.argument-dimensions;
    has @.argument-names;
    has @.last-argument-values;
    has $.working-precision is rw;
    has $.monitor is rw;
    has $.value;
    has UInt:D $!nan-count = 0;

    submethod BUILD(
                    :&!function!,
                    :@!argument-dimensions!,
                    :@!argument-names!,
                    :$!working-precision = Num,
                    :$!monitor = Nil,
                    :$!value = Nil) {
        @!last-argument-values = Whatever xx @!argument-names;
    }

    multi method new(:f(:func(:&function)), :prec(:$working-precision) = Num) {
        my $signature = &function.signature;
        my $num-args = $signature.params.elems;
        my @argument-names = $signature.params.map(*.name);
        my @argument-dimensions = $signature.params.map({ $_.name => $_.slurpy ?? 'List' !! 'Scalar' });
        self.bless(
                :&function,
                :@argument-dimensions,
                :@argument-names,
                :$working-precision,
                monitor => Nil,
                value => Nil
                );
    }

    multi method new(&function, $working-precision = Num) {
        self.new(:&function, :$working-precision)
    }

    multi method evaluate(*@args) { self.eval(@args) }
    multi method eval(*@args) {
        try {
            $!value = &!function(|@args);
        }

        # It seems it is a better to keep argument values with which the function failed to evaluate
        @!last-argument-values = @args;

        if $! {
            warn "Cannot evaluate numerical function at {@args.raku}.";
            return Whatever
        }

        if $!value ~~ Inf | -Inf || $!value.isNaN {
            # There should be counter how many times this warning was given.
            note "$!value was obtained for the arugments {@args.raku}." if $!nan-count < 3;
            $!nan-count++;
            # This message can be printed multiple times if this object is cloned or
            # similar Inf/NaN results are obtained by for other NumericalFunction objects
            # (associated with regions or not.)
            note "Suppressing further warnings for Inf, -Inf, and NaN results." if $!nan-count ≥ 3;

            # Is this a good idea?
            # It will avoid doing similar is-Inf / is-NaN checks, hence speed-up the computations.
            # return Whatever
        }

        return $!value
    }
}