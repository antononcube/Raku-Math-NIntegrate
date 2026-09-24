use v6.d;

use Math::NIntegrate::Rule;

# Leaf in the Composite design pattern
class Math::NIntegrate::Rule::General
        is Math::NIntegrate::Rule {
    # Integration rule data
    has @.abscissas is rw;
    has @.weights is rw;
    has @.error-weights is rw;
    has UInt:D $.dimension is rw = 1;

    # Evaluation data
    has Numeric $.integral;
    has Numeric $.error;
    has UInt $.largest-error-axis;

    #======================================================
    # Creators
    #======================================================

    submethod BUILD(:@!abscissas = Empty, :@!weights = Empty, :@!error-weights = Empty) {}

    multi method new(@abscissas, @weights, @error-weights) {
        self.bless(:@!abscissas, :@!weights, :@!error-weights);
    }

    multi method new(:@abscissas, :@weights, :@error-weights) {
        self.bless(:@!abscissas, :@!weights, :@!error-weights);
    }

    #| Copy attributes from an object
    method copy(Math::NIntegrate::Rule::General:D $from, Bool:D :deep(:deep-copy(:$clone)) = False) {
        # @!transforms are functions/subs, so, $from.transforms>>.clone.Array seems inappropriate.
        @!abscissas = $clone ?? $from.abscissas.clone !! $from.abscissas;
        @!weights = $clone ?? $from.weights.clone !! $from.weights;
        @!error-weights = $clone ?? $from.error-weights.clone !! $from.error-weights;

        return self
    }

    #| Clone the object
    method clone(-->Math::NIntegrate::Rule::General:D) {
        Math::NIntegrate::Rule::General.new.copy(self, :clone)
    }

    #======================================================
    # Integration
    #======================================================

    method integrate($region) {
        my $integralLocal = 0;
        my $errorLocal = 0;
        my @values;

        for ^@!abscissas.elems -> $i {
            my @point = @!abscissas[$i] ~~ Numeric:D ?? [@!abscissas[$i], ] !! @!abscissas[$i];
            my $value = $region.eval-integrand(@point);
            @values.push($value);
            $integralLocal += @!weights[$i] * $value;
            $errorLocal += @!error-weights[$i] * $value;
        }

        $!integral = $integralLocal;
        $!error = $errorLocal.abs;
        $!largest-error-axis = 0;

        return self;
    }
}