use v6.d;

class Math::NIntegrate::Processing::Actions::MethodSpec {

    method !canonical-name(Str:D $name --> Str:D) {
        my $compact = $name.trim.subst(/<[-_\s]>/, '', :g).lc;

        my %names =
                globaladaptive => 'GlobalAdaptive', localadaptive => 'LocalAdaptive', doubleexponential => 'DoubleExponential',
                montecarlo => 'MonteCarlo', adaptivemontecarlo => 'AdaptiveMonteCarlo', quasimontecarlo => 'QuasiMonteCarlo',
                adaptivequasimontecarlo => 'AdaptiveQuasiMonteCarlo',
                symbolicpiecewisesubdivision => 'SymbolicPiecewiseSubdivision', evenoddsubdivision => 'EvenOddSubdivision',
                trapezoidalrule => 'TrapezoidalRule', gausskronrodrule => 'GaussKronrodRule', lobattokronrodrule => 'LobattoKronrodRule',
                clenshawcurtisrule => 'ClenshawCurtisRule', montecarlorule => 'MonteCarloRule', cartesianrule => 'CartesianRule',
                maxpoints => 'max-points', maxrecursion => 'max-recursion', minrecursion => 'min-recursion',
                singularitydepth => 'singularity-depth', singularityhandler => 'singularity-handler',
                symbolicprocessing => 'symbolic-processing', nosingularityhandler => 'None',
                duffycoordinates => 'DuffyCoordinates', imt => 'IMT', irimorigutitakesawa => 'IMT';

        return %names{$compact} // $name.trim
    }

    method !merge-options(%spec is copy, $/ --> Hash:D) {
        for flat $<method-option>, $<numeric-option-spec>, $<singularity-handler-option>, $<cartesian-rule-method>, $<method-rule-spec>, $<option> -> $option {
            next unless $option.defined;
            my $pair = $option.made;
            %spec{$pair.key} = $pair.value;
        }
        %spec.Hash
    }

    method strategy-symbol($/) {
        make { name => self!canonical-name($/.Str), type => 'strategy' };
    }

    method decorator-strategy-symbol($/) {
        make { name => self!canonical-name($/.Str), type => 'strategy' };
    }

    method rule-symbol($/) {
        make { name => self!canonical-name($/.Str), type => 'rule' };
    }

    method method-symbol($/) {
        make $<strategy-symbol> ?? $<strategy-symbol>.made !! $<decorator-strategy-symbol> ?? $<decorator-strategy-symbol>.made !! $<rule-symbol>.made;
    }

    method numeric-option-spec($/) {
        make self!canonical-name($<numeric-option-symbol>.Str) => $<number>.Int;
    }

    method singularity-handler-option($/) {
        make 'singularity-handler' => self!canonical-name($<singularity-handler-symbol>.Str);
    }

    method option($/) {
        my $value = ($<option-value> // $<number>).Str;
        make self!canonical-name($<symbol>.Str) => ($value ~~ /^\d+$/ ?? $value.Int !! self!canonical-name($value));
    }

    method method-option-value($/) {
        my %spec = $<method-symbol>.made.Hash; make self!merge-options(%spec, $/);
    }

    method method-option($/) {
        make 'method' => $<method-option-value>.made;
    }

    method strategy-spec($/) {
        my %spec = $<strategy-symbol>.made.Hash;
        make self!merge-options(%spec, $/);
    }

    method rule-component-spec($/) {
        my %spec = $<rule-symbol>.made.Hash; make self!merge-options(%spec, $/);
    }

    method rule-component-sequence($/) {
        make $<rule-component-spec>.map(*.made).Array;
    }

    method rule-component-option-value($/) {
        make $<rule-component-sequence> ?? $<rule-component-sequence>.made !! $<rule-component-spec>.made;
    }

    method cartesian-rule-method($/) {
        make 'method' => $<rule-component-option-value>.made;
    }

    method cartesian-rule-spec($/) {
        my %spec = name => self!canonical-name($<cartesian-rule>.Str), type => 'rule';
        make self!merge-options(%spec, $/);
    }

    method rule-spec($/) {
        make $<cartesian-rule-spec> ?? $<cartesian-rule-spec>.made !! $<rule-component-spec>.made;
    }

    method method-rule-spec($/) {
        make 'method' => $<rule-spec>.made;
    }

    method TOP($/) { make $<method-option>.made; }
}
