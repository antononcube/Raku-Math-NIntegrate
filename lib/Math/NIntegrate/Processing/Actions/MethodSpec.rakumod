use v6.d;

class Math::NIntegrate::Processing::Actions::MethodSpec {

    # Instead of having a dimension attribute using Whatever for the rule.
    # Then Builder's method make-rule() can apply its own default.
    has %!default =
            name => 'None',
            type => 'decorator',
            method => {name => 'GlobalAdaptive', type => 'strategy',
                       method => {name => Whatever, type => 'rule'}};

    has Bool:D $.full-spec is rw = False;

    method get-default() {
       return %!default.deepmap: *.clone;
    }

    method !canonical-name(Str:D $name --> Str:D) {
        my $compact = $name.trim.subst(/<[-_\s]>/, '', :g).subst(/ ^ <["\']> | <["\']> $/, :g).lc;

        my %names =
                globaladaptive => 'GlobalAdaptive', localadaptive => 'LocalAdaptive', doubleexponential => 'DoubleExponential',
                montecarlo => 'MonteCarlo', adaptivemontecarlo => 'AdaptiveMonteCarlo', quasimontecarlo => 'QuasiMonteCarlo',
                adaptivequasimontecarlo => 'AdaptiveQuasiMonteCarlo',
                symbolicpiecewisesubdivision => 'SymbolicPiecewiseSubdivision', evenoddsubdivision => 'EvenOddSubdivision',
                trapezoidalrule => 'TrapezoidalRule', newtoncotesrule => 'NewtonCotesRule',
                gausskronrodrule => 'GaussKronrodRule', lobattokronrodrule => 'LobattoKronrodRule',
                clenshawcurtisrule => 'ClenshawCurtisRule', montecarlorule => 'MonteCarloRule', cartesianrule => 'CartesianRule',
                multidimensionalrule => 'MultiDimensionalRule',
                maxpoints => 'max-points', maxrecursion => 'max-recursion', minrecursion => 'min-recursion',
                singularitydepth => 'singularity-depth', singularityhandler => 'singularity-handler',
                romberg => 'romberg-quadrature', rombergquadrature => 'romberg-quadrature',
                symbolicprocessing => 'symbolic-processing', nosingularityhandler => 'None',
                partitioning => 'Partitioning', initialestimaterelaxation => 'InitialSstimateRelaxation', randomseed => 'RandomSeed',
                duffycoordinates => 'DuffyCoordinates', imt => 'IMT', irimorigutitakesawa => 'IMT', doubleexponent => 'DoubleExponent',
                unitcuberescaling => 'UnitCubeRescaling',
                whatever => Whatever
                ;

        return %names{$compact} // $name.trim
    }

    method !merge-options(%spec is copy, $/ --> Hash:D) {
        for flat $<method-option>, $<numeric-option-spec>, $<boolean-option-spec>, $<singularity-handler-option>,
                $<cartesian-rule-method>, $<method-rule-spec>, $<method-strategy-spec>, $<partitioning-option>, $<option>
        -> $option {
            next unless $option.defined;
            my $pair = $option.made;
            %spec{$pair.key} = $pair.value;
        }
        %spec.Hash
    }

    method boolean($/) {
        make $/.lc eq 'true'
    }

    method number($/) {
        make $/.Int
    }

    method number-list($/) {
        make $<number>>>.Int.List
    }

    method strategy-symbol($/) {
        make { name => self!canonical-name($/.Str), type => 'strategy' };
    }

    method strategy-decorator-symbol($/) {
        make { name => self!canonical-name($/.Str), type => 'decorator' };
    }

    method rule-symbol($/) {
        make { name => self!canonical-name($/.Str), type => 'rule' };
    }

    method method-symbol($/) {
        make $<strategy-symbol> ?? $<strategy-symbol>.made !! $<strategy-decorator-symbol> ?? $<strategy-decorator-symbol>.made !! $<rule-symbol>.made;
    }

    method boolean-option-spec($/) {
        make self!canonical-name($<boolean-option-symbol>.Str) => $<boolean> ?? $<boolean>.made !! True;
    }

    method numeric-option-spec($/) {
        make self!canonical-name($<numeric-option-symbol>.Str) => $<number>.Int;
    }

    method singularity-handler-option($/) {
        make 'singularity-handler' => self!canonical-name($<singularity-handler-symbol>.Str);
    }

    method partitioning-option($/) {
        make 'partitioning' => $<number-list> ?? $<number-list>.made !! $<number>.made;
    }

    method option($/) {
        my $value = ($<option-value> // $<number-list> // $<number>).Str;
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
        if $!full-spec { %spec = %spec , {method => %!default<method><method>} }
        make self!merge-options(%spec, $/);
    }

    method rule-component-spec($/) {
        my %spec = $<rule-symbol>.made.Hash;
        make self!merge-options(%spec, $/);
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

    method method-strategy-spec($/) {
        make 'method' => $<strategy-spec>.made;
    }

    method strategy-decorator-spec($/) {
        my %spec = $<strategy-decorator-symbol>.made.Hash;
        if $!full-spec { %spec = %spec , {method => %!default<method>} }
        make self!merge-options(%spec, $/);
    }

    method TOP($/) {
        my %parsed = $/.values.head.made;

        if $!full-spec {
            my %res = self.get-default;
            if %parsed<type> eq 'rule' {
                %res<method><method> = %parsed
            } elsif %parsed<type> eq 'strategy' {
                %res<method> = %parsed
            } elsif %parsed<type> eq 'decorator' && %parsed<method><type> eq 'rule' {
                # There should be a more elegant way of doing this. I.e. via the class methods.
                my %s = %res<method>;
                %s<method> = %parsed<method>;
                %res = %parsed;
                %res<method> = %s;
            } else {
                %res = %parsed
            }
            make %res
        } else {
            make %parsed
        }
    }
}
