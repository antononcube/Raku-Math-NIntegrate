use v6.d;

use Math::NIntegrate::Rule::General;

class Math::NIntegrate::Rule::GaussKronrod
        is Math::NIntegrate::Rule::General {
    has Int $.points;

    method new(Int :$!points) {
        my %data = self!generate-data($!points);
        Math::NIntegrate::Rule::General.new(:%data<abscissas>, %data<weights>, %data<error-weights>);
    }

    method !generate-data(Int $points) {
        given $points {
            when 5 {
                return %(
                    abscissas => [0.00795732, 0.0469101, 0.122917, 0.230765, 0.360185, 0.5, 0.639815, 0.769235, 0.877083, 0.95309, 0.992043],
                    weights => [0.021291, 0.0576167, 0.0934004, 0.12052, 0.136425, 0.141494, 0.136425, 0.12052, 0.0934004, 0.0576167, 0.021291],
                    error-weights => [0.021291, -0.0608468, 0.0934004, -0.118794, 0.136425, -0.142951, 0.136425, -0.118794, 0.0934004, -0.0608468, 0.021291]
                );
            }
            when 6 {
                return %(
                    abscissas => [0.0056484, 0.0337652, 0.0893133, 0.169395, 0.268441, 0.38069, 0.5, 0.61931, 0.731559, 0.830605, 0.910687, 0.966235, 0.994352],
                    weights => [0.0151981, 0.0418472, 0.0686603, 0.090536, 0.106605, 0.116885, 0.120536, 0.116885, 0.106605, 0.090536, 0.0686603, 0.0418472, 0.0151981],
                    error-weights => [0.0151981, -0.043815, 0.0686603, -0.0898448, 0.106605, -0.117072, 0.120536, -0.117072, 0.106605, -0.0898448, 0.0686603, -0.043815, 0.0151981]
                );
            }
            default {
                die "Unsupported number of points: $points";
            }
        }
    }
}
