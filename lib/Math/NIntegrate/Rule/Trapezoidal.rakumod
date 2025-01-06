use v6.d;

use Math::NIntegrate::Utilities;
use Math::NIntegrate::Rule::General;

class Math::NIntegrate::Rule::Trapezoidal
        is Math::NIntegrate::Rule::General {

    has $.points;

    submethod BUILD(UInt:D :$!points, :$wprec = Num) {
        my %res = self.make-weights($!points, $wprec);
        self.abscissas = |%res<abscissas>;
        self.weights = |%res<weights>;
        self.error-weights = |%res<error-weights>;
    }

    multi method new($points, $wprec = Num) {
        self.bless(:$points, :$wprec)
    }

    multi method new(:$points, :$wprec = Num) {
        self.bless(:$points, :$wprec)
    }

    method make-weights(UInt:D $points, $wprec) {
        my $n1 = $points;

        my $n = 2 * $n1 - 1;

        my @abscissas;
        my @weights;
        my @error-weights;

        my $h = 1 / (1 * ($n - 1));
        $h = numerical($h, $wprec);

        for (^$n) -> $i {
            my $vf = numerical($i * $h, $wprec);
            @abscissas[$i] = $vf;
        }

        @weights = 1 xx $n;
        for (1 ... ($n - 2)) -> $i {
            @weights[$i] = numerical(1 / ($n - 1), $wprec);
        }

        @weights[0] = numerical(1 / 2 / ($n - 1), $wprec);
        @weights[$n - 1] = numerical(1 / 2 / ($n - 1), $wprec);

        @error-weights = 0 xx $n;
        for (2, 4 ... ($n - 2)) -> $i {
            @error-weights[$i] = numerical(1 / ($n1 - 1), $wprec);
        }
        @error-weights[0] = numerical(1 / 2 / ($n1 - 1), $wprec);
        @error-weights[$n-1] = numerical(1 / 2 / ($n1 - 1), $wprec);

        @error-weights = @weights Z- @error-weights;

        return %(:@abscissas, :@weights, :@error-weights);
    }
}
