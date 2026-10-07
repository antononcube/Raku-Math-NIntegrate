use v6.d;

use Math::NIntegrate::Rule;
use Math::NIntegrate::Strategy;

grammar Math::NIntegrate::Processing::Grammar {

    rule TOP { <strategy-decorator-spec> | <strategy-spec> | <rule-spec> }

    token symbol { <.quote>? <[\w \- _]>+ <.quote>? }

    token boolean {:i True | False | '!' }

    token number { \d+ | <[-+]>? Inf }

    rule number-list { <.lb> <number>+ % <.sep> <.rb>}

    # Conditional parsing with different predicates.
    # The regex predicates are probably redundant.
    # Probably it is a good idea to check similarity with known-symbols.

    rule top-level-strategy {
        (<symbol>)
        <?{
            my $t = $0.Str.subst(/ <["']> /):g;
            try ::($t);
            my $exists = $! ?? False !! True;
            given $t {
                when $_ ~~ / [ 'Rule' | '-rule' | 'Rescaling' | '-rescaling' ] <["']>? $/ { False }
                when !$exists { note "⎡$_⎦ does not exist."; False}
                when ::($_) ~~ Math::NIntegrate::Strategy { True }
                when ::($_) ~~ Math::NIntegrate::Rule { False }
                default {
                    note "If {$_} is a top level strategy then it is expected to be of type Math::NIntegrate::Strategy.";
                    False
                }
            }
        }> }

    rule top-level-rule {
        (<symbol>)
        <?{
            my $t = $0.Str.subst(/ <["']> /):g;
            try ::($t);
            my $exists = $! ?? False !! True;
            given $t {
                when $_ ~~ / [ 'Rescaling' | '-rescaling' ] <["']>? $/ { False }
                when !$exists { note "⎡$_⎦ does not exist."; False}
                when ::($_) ~~ Math::NIntegrate::Rule { True }
                when ::($_) ~~ Math::NIntegrate::Strategy { False }
                default {
                    note "If {$_} is a top level rule then it is expected to be of type Math::NIntegrate::Rule.";
                    False
                }
            }
        }> }

    rule option-value { <symbol> | <number> }

    regex option {
        # I do not like this conditional chack, but it was quick to do with the generic <method-option>
        # <!{ $<symbol>.Str.lc eq 'method' }>
        [ <symbol> \s* <.arrow> \s* <option-value> | ':' <number> <symbol> | ':' <symbol> <.lb> \s* <number> \s* <.rb>]
    }

    rule option-sequence { <option>+ % <.sep> }

    regex method-option {
     | <method> \s* <.arrow> \s* <method-option-value>
     | ':' <method> <.lb> \s* <method-option-value> \s* <.rb>
     | '"' <method> '"' \s* ':' \s* <method-option-value>
    }

    # This rule is too general: parses correct syntax, but does not impose only allowed sub-options.
    regex method-option-value {
        || <method-symbol>
        || <.lb> \s* <method-symbol> \s* <.rb>
        || <.lb> \s* <method-symbol> \s* <.sep> \s* [ <method-option> | <option> ]+ % [\s* <.sep> \s*] \s* <.rb>
    }

    rule method-symbol {
        <strategy-symbol> | <strategy-decorator-symbol> | <rule-symbol>
    }

    rule strategy-symbol-known {
        | <global-adaptive> | <local-adaptive> | <double-exponential>
        | <monte-carlo> | <adaptive-monte-carlo> | <quasi-monte-carlo>
        | <adaptive-quasi-monte-carlo>
    }

    rule strategy-symbol {
        || '"' ~ '"' <strategy-symbol-known>
        || <strategy-symbol-known>
        || <top-level-strategy>
    }

    rule strategy-decorator-known {
        <unit-cube-rescaling> | <symbolic-piecewise-subdivision> | <even-odd-subdivision>
    }

    rule strategy-decorator-symbol {
        || '"' ~ '"' <strategy-decorator-known>
        || <strategy-decorator-known>
    }

    rule rule-symbol-known {
        | <trapezoidal-rule> | <newton-cotes-rule>
        | <gauss-kronrod-rule> | <lobatto-kronrod-rule>
        | <clenshaw-curtis-rule> | <monte-carlo-rule>
        | <multidimensional-rule>
        | <cartesian-rule>
        | <whatever>
    }

    rule rule-symbol {
        || '"' ~ '"' <rule-symbol-known>
        || <rule-symbol-known>
        || <top-level-rule>
    }

    token numeric-option-symbol {
        <max-recursion> | <max-points> | <min-recursion> | <singularity-depth> | <points> | <generators> | <random-seed>
    }

    regex numeric-option-spec {
        | <numeric-option-symbol> \s* <.arrow> \s* <number>
        | ':' <number> <numeric-option-symbol>
        | ':' <numeric-option-symbol> <.lb> \s* <number> \s* <.rb>
        | '"' <numeric-option-symbol> '"' \s* ':' \s* <number>
    }

    token boolean-option-symbol {
        <romberg-quadrature>
    }

    regex boolean-option-spec {
        | <boolean-option-symbol> \s* <.arrow> \s* <boolean>
        | ':' <boolean-option-symbol>
        | ':' <boolean> <boolean-option-symbol>
        | ':' <numeric-option-symbol> <.lb> \s* <boolean> \s* <.rb>
        | '"' <numeric-option-symbol> '"' \s* ':' \s* <boolean>
    }

    token singularity-handler-symbol {
        | <imt>
        | <double-exponential>
        | <duffy-coordinates>
        | <no-singularity-handler>
        | <whatever>
    }

    regex singularity-handler-option {
        | <singularity-handler> \s* <.arrow> \s* <.quote>? <singularity-handler-symbol> <.quote>?
        | ':' <singularity-handler> <.lb> \s* <.quote>? <singularity-handler-symbol> <.quote>? \s* <.rb>
        | '"' <singularity-handler> '"' \s* ':' \s* <singularity-handler-symbol>
    }

    regex partitioning-option {
        | <partitioning> \s* <.arrow> \s* [ <number> | <number-list> ]
        | ':' <partitioning> <.lb> \s* [ <number> | <number-list> ] \s* <.rb>
        | ':' <number> <partitioning>
        | '"' <partitioning> '"' \s* ':' \s* [ <number> | <number-list> ]
    }

    rule strategy-spec {
        || <strategy-symbol>
        || <.lb> <strategy-symbol> [ <.sep> [ <method-rule-spec> | <partitioning-option> | <numeric-option-spec> | <singularity-handler-option> ]* % <.sep> ]? <.rb>
    }

    regex rule-component-spec {
        || <rule-symbol>
        || <.lb> \s* <rule-symbol> \s* [ <.sep> \s* [ <boolean-option-spec> | <numeric-option-spec> ]* % [ \s* <.sep> \s*] ]? <.rb>
    }

    regex rule-component-sequence {
        <rule-component-spec>+ % [\s* <.sep> \s*]
    }

    regex rule-component-option-value {
        <.lb> \s* <rule-component-sequence> \s* <.rb> || <rule-component-spec>
    }

    regex cartesian-rule-method {:i
        | <method> \s* <.arrow> \s* <rule-component-option-value>
        # It really does not want me to use <.lb> instead of <.punct> !!
        | ':' <method> <.punct> \s* <rule-component-option-value> \s* <.rb>
        | '"' <method> '"' \s* ':' \s* <.lb> \s* <rule-component-option-value> \s* <.rb>
    }

    # This might be not a sufficiently strong production rule.
    # Cartesian rule is not allowed to have method that is a Cartesian rule.
    # (Although, in principle, that can be processed and corresponding "flattened" Cartesian rule be created.)
    rule cartesian-rule-spec {
        || <cartesian-rule>
        || '"' ~ '"' <cartesian-rule>
        || <.lb> [ <cartesian-rule> | '"' ~ '"' <cartesian-rule> ] [ <.sep> [ <cartesian-rule-method> | <numeric-option-spec> ]* % <.sep> ]? <.rb>
    }

    rule rule-spec {
        || <cartesian-rule-spec>
        || <rule-component-spec>
    }

    # Method rule spec has to be separate -- a (leaf) strategy cannot have another strategy as method.
    regex method-rule-spec {
        | <method> \s* <.arrow> \s* <rule-spec>
        | ':' <method> <.lb> \s* <rule-spec> \s* <.rb>
        | '"' <method> '"' \s* ':' \s* <rule-spec>
    }

    regex method-strategy-spec {
        | <method> \s* <.arrow> \s* <strategy-spec>
        | ':' <method> <.lb> \s* <strategy-spec> || <rule-spec> \s* <.rb>
        | '"' <method> '"' \s* ':' \s* <strategy-spec>
    }

    rule strategy-decorator-spec {
        || <strategy-decorator-symbol>
        || <.lb> <strategy-decorator-symbol> [ <.sep> [ <method-strategy-spec> | <method-rule-spec> | <numeric-option-spec> ]* % <.sep> ]? <.rb>
    }

    # Brackets
    token lb { '(' }
    token rb { ')' }

    # Arrow
    token arrow { '=>' }
    
    # Separator
    token sep { ',' }

    # Quote
    token quote { <["']> }

    # Option names
    token generators {:i generators }
    token initial-estimate-relaxation { InitialEstimateRelaxation | initial <[\-_]> estimate <[\-_]> relaxation }
    token max-points { MaxPoints | max <[\-_]> points }
    token max-recursion { MaxRecursion | max <[\-_]> recursion }
    token method {:i method }
    token min-recursion { MinRecursion | min <[\-_]> recursion }
    token partitioning {:i partitioning }
    token points {:i points }
    token random-seed {:i RandomSeed | random <[\-_]> seed }
    token romberg-quadrature {:i RombergQuadrature | Romberg | romberg <[\-_]> quadrature }
    token singularity-depth { SingularityDepth | singularity <[\-_]> depth }
    token singularity-handler { SingularityHandler | singularity <[\-_]> handler }
    token symbolic-processing { SymbolicProcessing | symbolic <[\-_]> processing }

    # Automatic and Whatever
    token whatever { Whatever | WhateverCode }
    token automatic {:i 'automatic' }

    # Integration strategy names
    token adaptive-monte-carlo { AdaptiveMonteCarlo | adaptive <[_\-]> monte <[_\-]> carlo }
    token adaptive-quasi-monte-carlo { AdaptiveQuasiMonteCarlo | adaptive <[_\-]> quasi <[_\-]> monte <[_\-]> carlo }
    token double-exponential { DoubleExponential | double <[_\-]> exponential }
    token global-adaptive { GlobalAdaptive | global <[_\-]> adaptive }
    token local-adaptive { LocalAdaptive | local <[_\-]> adaptive }
    token monte-carlo { MonteCarlo | monte <[_\-]> carlo }
    token quasi-monte-carlo { QuasiMonteCarlo | quasi <[_\-]> monte <[_\-]> carlo }

    # Integration preprocessor names
    token even-odd-subdivision { EvenOddSubdivision | even <[_\-]> odd <[_\-]> subdivision }
    token symbolic-piecewise-subdivision { SymbolicPiecewiseSubdivision | symbolic <[_\-]> piecewise <[_\-]> subdivision }
    token unit-cube-rescaling { UnitCubeRescaling | unit <[_\-]> cube <[_\-]> rescaling }

    # Integration rule names
    token cartesian-rule { CartesianRule | cartesian <[_\-]> rule }
    token clenshaw-curtis-rule { ClenshawCurtisRule | clenshaw <[_\-]> curtis <[_\-]> rule }
    token gauss-kronrod-rule { GaussKronrodRule | gauss <[_\-]> kronrod <[_\-]> rule }
    token lobatto-kronrod-rule { LobattoKronrodRule | lobatto <[_\-]> kronrod <[_\-]> rule }
    token monte-carlo-rule { MonteCarloRule | monte <[_\-]> carlo <[_\-]> rule }
    token multidimensional-rule { Multi [d|D] imensionalRule | multidimensional <[_\-]> rule }
    token newton-cotes-rule { NewtonCotesRule | newton <[_\-]> cotes <[_\-]> rule }
    token trapezoidal-rule { TrapezoidalRule | trapezoidal <[_\-]> rule }

    # Singularity handler names
    token imt {:i imt | iri <[\-_]> moriguti <[\-_]> takesawa }
    token no-singularity-handler {:i none | no <[\-_]> singularity <[\-_]> handler }
    token duffy-coordinates {:i DuffyCoordinates | duffy <[\-_]> coordinates } # This is, actually, a preprocessor
}
