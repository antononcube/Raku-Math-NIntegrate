use v6.d;

class Math::NIntegrate::Rule::General {
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

    method approximate-integral($simplex --> Map:D) {
        return %();
    }
}
