use v6.d;

use Math::NIntegrate::Rule::General;

class Math::NIntegrate::Rule::NewtonCotes
        is Math::NIntegrate::Rule::General {
    has $.points;
    has Bool:D $is-open = False;

    submethod BUILD(UInt:D :$!points, Bool:D :$!is-open = False, :$wprec = Num) {
        my %res = self.make-weights($!points, $!is-open, $wprec);
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

    method make-weights(UInt:D $points, Bool:D $is-open = False, $wprec = Num) {!!!}
}
