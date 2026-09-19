use v6.d;

class Math::NIntegrate::NumericalFunction {
    has &.function is rw;
    has @.argument-dimensions;
    has @.argument-names;
    has @.last-argument-values;
    has $.working-precision is rw;
    has $.monitor is rw;
    has $.value;

    submethod BUILD(
                    :&!function!,
                    :@!argument-dimensions!,
                    :@!argument-names!,
                    :$!working-precision = Num,
                    :$!monitor = Nil,
                    :$!value = Nil) {
        @!last-argument-values = Whatever xx @!argument-names;
    }

    multi method new(:f(:func(:&function)), :wprec(:$working-precision) = Num) {
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

        return $!value
    }
}