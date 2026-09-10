use v6.d;

use Math::NIntegrate::Rule;

class Math::NIntegrate::Rule::General
        is Math::NIntegrate::Rule {
    has @.abscissas is rw;
    has @.weights is rw;
    has @.error-weights is rw;

    submethod BUILD(:@!abscissas = Empty, :@!weights = Empty, :@!error-weights = Empty) { }

    multi method new(@abscissas, @weights, @error-weights) {
        self.bless(:@!abscissas, :@!weights, :@!error-weights);
    }

    multi method new(:@abscissas, :@weights, :@error-weights) {
        self.bless(:@!abscissas, :@!weights, :@!error-weights);
    }

    method integrate($simplex --> Map:D) {
        return %();
    }
}
