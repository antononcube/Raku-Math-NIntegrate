use v6.d;

use Math::NIntegrate::Utilities;
use Math::NIntegrate::Rule::General;

class Math::NIntegrate::Rule::Trapezoidal
        is Math::NIntegrate::Rule::General {

    has $.points is required;
    has Bool:D $.romberg-quadrature = True;

    submethod BUILD(UInt:D :$!points, :$working-precision = Num, :$!romberg-quadrature = True) {
        my %res = self.make-weights($!points, $working-precision, $!romberg-quadrature);
        self.abscissas = |%res<abscissas>;
        self.weights = |%res<weights>;
        self.error-weights = |%res<error-weights>;
    }

    multi method new($points, $working-precision = Num, :romberg(:$romberg-quadrature) = True) {
        self.bless(:$points, :$working-precision, :$romberg-quadrature)
    }

    multi method new(:$points, :$working-precision = Num, :romberg(:$romberg-quadrature) = True) {
        self.bless(:$points, :$working-precision, :$romberg-quadrature)
    }

    method make-weights(UInt:D $points, $working-precision = Num, Bool:D $romberg-quadrature = True) {

        my $n1 = numerical($points, $working-precision);
        my $n = numerical(2 * $n1 - 1, $working-precision);

        my @abscissas;
        my @weights;
        my @error-weights;

        my $h = 1 / (1 * ($n - 1));
        $h = numerical($h, $working-precision);

        for (^$n) -> $i {
            my $vf = numerical($i * $h, $working-precision);
            @abscissas[$i] = $vf;
        }

        @weights = 1 xx $n;
        for (1 ... ($n - 2)) -> $i {
            @weights[$i] = numerical(1 / ($n - 1), $working-precision);
        }

        @weights[0] = numerical(1 / 2 / ($n - 1), $working-precision);
        @weights[$n - 1] = numerical(1 / 2 / ($n - 1), $working-precision);

        @error-weights = 0 xx $n;
        for (2, 4 ... ($n - 2)) -> $i {
            @error-weights[$i] = numerical(1 / ($n1 - 1), $working-precision);
        }
        @error-weights[0] = numerical(1 / 2 / ($n1 - 1), $working-precision);
        @error-weights[$n-1] = numerical(1 / 2 / ($n1 - 1), $working-precision);

        if $romberg-quadrature {
            my @ws = @weights;
            @weights = (( 4 <<*>> @weights ) Z- @error-weights) <</>> 3;
            @error-weights = @weights Z- @ws
        } else {
            @error-weights = @weights Z- @error-weights;
        }

        return %(:@abscissas, :@weights, :@error-weights);
    }

    #======================================================
    # Representation
    #======================================================
    multi method Hash(::?CLASS:D:-->Map:D) {
        my %h = self.Math::NIntegrate::Rule::General::Hash;
        %h = |%h, :$!romberg-quadrature;
        return %h
    }
}
