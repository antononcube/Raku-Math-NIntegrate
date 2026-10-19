use v6.d;

use Math::NIntegrate::Strategy::Decorator;

class Math::NIntegrate::Strategy::UnitCubeRescaling
        is Math::NIntegrate::Strategy::Decorator {

    # Whether the unit-cube rescaling is applied already or not.
    # I am not sure that it is needed since unit-cube rescaling an idempotent.
    has Bool:D $!rescaled = False;

    # This decorator class replaces the integrands and ranges of component's regions.
    # Since regions can have different integrands and/or ranges,
    # each region is treated separately.

    # In principle if at creation the component has regions, they can be modified at creation.
    # But it is more flexible to it as an added behavior.
    #submethod TWEAK(:$component, :@ranges) {
    #}

    #======================================================
    # Rescaling functions
    #======================================================

    sub fId($x) { $x }
    sub fIdJac($x) { 1 }

    sub fFin($x, $a, $b) { $a + ($b - $a) * $x }
    sub fFinJac($x, $a, $b) { $b - $a }

    sub fInf($x, $a) { $a + 1 / (1 - $x) - 1 }
    sub fInfJac($x, $a) { 1 / ($x - 1) ** 2}

    sub fNInf($x, $b) { 1 + $b - 1 / $x }
    sub fNInfJac($x, $b) { 1 / $x ** 2}

    sub fInfInf($x) { 1 / (1 - $x) - 1 / $x }
    sub fInfInfJac($x) {
        1 / (1 - $x) ** 2 + 1 / $x ** 2
    }

    #======================================================
    # Rescale by ranges
    #======================================================

    sub make-variable(Str:D $var-name, Str:D :$twigil='^') {
        with $var-name ~~ / ^ '$'? $<id>=(<alpha> [<alnum>|'-']*) $ / { "\${$twigil}{$<id>}" }
    }

    sub make-replacement(@range, Str:D :$twigil='^') {

        # Should numeric checks be done
        die 'The first argument is expected to be a list of three or more elements and the first element to be a string.'
        unless @range.elems ≥ 3 && @range.head ~~ Str:D;

        my $var = make-variable(@range.head, :$twigil);

        return do given @range {
            when $_[1] == 0 && $_.tail == 1 {
                { :$var, func => $var, jac => '1' }
            }

            when $_[1].abs < Inf && $_.tail.abs < Inf {
                { :$var, func => "fFin($var, {$_[1]}, {$_.tail})", jac => "fFinJac($var, {$_[1]}, {$_.tail})" }
            }

            when $_[1].abs < Inf && $_.tail eq Inf {
                { :$var, func => "fInf($var, {$_[1]})", jac => "fInfJac($var, {$_[1]})" }
            }

            when $_[1] eq Inf && $_.tail.abs < Inf {
                { :$var, func => "-fInf($var, {$_.tail})", jac => "(-fInfJac($var, {$_.tail}))" }
            }

            when $_[1] eq -Inf && $_.tail.abs < Inf {
                { :$var, func => "fNInf($var, {$_.tail})", jac => "fNInfJac($var, {$_.tail})" }
            }

            when $_[1].abs < Inf && $_.tail eq -Inf {
                { :$var, func => "-fNInf($var, {$_.tail})", jac => "(-fNInfJac($var, {$_.tail}))" }
            }

            when $_[1] eq -Inf && $_.tail eq Inf {
                { :$var, func => "fInfInf($var)", jac => "fInfInfJac($var)" }
            }

            when $_[1] eq Inf && $_.tail eq -Inf {
                { :$var, func => "-fInfInf($var)", jac => "(-fInfInfJac($var))" }
            }

            default {
                die 'The interation strategy decorator UnitCubeRescaling can handle only numeric ranges.'
            }
        }
    }

    # Function checks to be done: arity, variables to ranges correspondence
    sub make-rescaled-func(
            :&f,
            :@ranges where @ranges.all ~~ (Array:D | List:D | Seq:D),
            Bool:D :$eval = True
                           ) {
        my @mappings = @ranges.map(*.&make-replacement(twigil => ''));
        my $vars = @mappings.map(*<var>).join(', ');
        my $args = @mappings.map(*<func>).join(', ');
        my $jac = @mappings.grep(* ne '1').map(*<jac>).join(' * ');
        my $code = "sub ($vars) \{ \&f($args) * $jac \}";

        return do if $eval {
            use MONKEY-SEE-NO-EVAL;
            EVAL $code
        } else { $code }
    }

    #======================================================
    # Rescale integrand
    #======================================================

    # This can be an utility.
    multi method ranges-to-cube(&f, @ranges, Bool:D :$eval = True) {
        self.ranges-to-cube(:&f, :@ranges, :$eval)
    }

    multi method ranges-to-cube(:f(:&integrand), :@ranges, Bool:D :$eval = True) {
        my $newIntegrand = make-rescaled-func(f=> &integrand, :@ranges, :$eval);
        my @newRanges = [0, 1] xx @ranges.elems;
        return {integrand => $newIntegrand, ranges => @newRanges}
    }


    #======================================================
    # Algorithm
    #======================================================

    #| UnitCubeRescaling decorator strategy's algorithm
    method algorithm(*%args --> Map:D) {
        # Conditional rescaling should not be used since since the regions are rw in Math::NIntegrate::Strategy.
        # Unit cube rescaling, though, is idempotent, hence it can be bravely applied/invoked.

        # At this point rescaling to the unit cube is already done.
        # It it a interesting question should the original regions of the component be kept
        # and this object have the rescaled ones. The method rescale can substitute the former with the latter.
        # self.rescale

        #note 'Call in Math::NIntegrate::Strategy::UnitCubeRescaling::algorithm';
        self.Math::NIntegrate::Strategy::Decorator::algorithm(|%args)
    }
}
