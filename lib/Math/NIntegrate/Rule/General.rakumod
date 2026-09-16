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