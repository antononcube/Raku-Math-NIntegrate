use v6.d;

grammar Math::NIntegrate::Processing::Grammar {

    rule TOP { <method-option> }

    token symbol { <[\w \- _]>+ }

    token number { \d+ }

    rule top-level-strategy { <symbol> }

    rule top-level-rule { <symbol> }

    rule option-value { <symbol> | <number> }

    regex option {
        # I do not like this conditional chack, but it was quick to do with the generic <method-option>
        [ <symbol> \s* <arrow> \s* <option-value> | ':' <number> <symbol> ] <!{ $<symbol>.Str.lc eq 'method' }>
    }

    rule option-sequence { <option>+ % <sep> }

    regex method-option {:i
     | <method> \s* <arrow> \s* <method-option-value>
     | ':' <method> <lb> \s* <method-option-value> \s* <rb>
     | '"' <method> '"' \s* ':' \s* <method-option-value>
    }

    # This rules is too general: parses correct syntax, but does not impose only allowed sub-options.
    rule method-option-value {
        | <method-symbol>
        | <lb> <method-symbol> <rb>
        | <lb> <method-symbol> <sep> [ <method-option> | <option> ]+ % <sep> <rb>
    }

    rule method-symbol {
        <strategy-symbol> | <decorator-strategy-symbol> | <rule-symbol>
    }

    rule strategy-symbol-known {
        | <global-adaptive> | <local-adaptive> | <double-exponential>
        | <monte-carlo> | <adaptive-monte-carlo> | <quasi-monte-carlo>
        | <adaptive-quasi-monte-carlo>
    }

    rule strategy-symbol {
        | <strategy-symbol-known>
        | <top-level-strategy>
    }

    rule decorator-strategy-symbol {
        <symbolic-piecewise-subdivision> | <even-odd-subdivision>
    }

    rule rule-symbol-known {
        | <trapezoidal-rule>
        | <gauss-kronrod-rule>  | <lobatto-kronrod-rule>
        | <clenshaw-curtis-rule>  | <monte-carlo-rule>
        | <cartesian-rule>
    }

    rule rule-symbol {
        | <rule-symbol-known>
        | <top-level-rule>
    }

    rule singularity-handler-symbol {
        | <imt>
        | <double-exponential>
        | <duffy-coordinates>
        | <no-singularity-handler>
    }

    token numeric-option-symbol { <max-recursion> | <max-points> | <min-recursion> | <singularity-depth> | <points> }

    regex numeric-option-spec {
        | <numeric-option-symbol> \s* <arrow> \s* <number>
        | ':' <number> <numeric-option-symbol>
        | ':' <numeric-option-symbol> <lb> \s* <number> \s* <rb>
        | '"' <numeric-option-symbpl> '"' \s* ':' \s* <number>
    }

    regex singularity-handler-option {
        | <singularity-handler> \s* <arrow> \s* <singularity-handler-symbol>
        | ':' <singularity-handler> <lb> \s* <singularity-handler-symbol> \s* <rb>
        | '"' <singularity-handler> '"' \s* ':' \s* <singularity-handler-symbol>
    }

    rule strategy-spec {
        || <strategy-symbol>
        || <lb> <strategy-symbol> [ <sep> [ <method-option> | <numeric-option-spec> | <singularity-handler-option> ]* % <sep> ]? <rb>
    }

    rule rule-component-spec {
        || <rule-symbol>
        || <lb> <rule-symbol> [ <sep> <numeric-option-spec>* % <sep> ]? <rb>
    }

    regex rule-component-sequence {
        <rule-component-spec>+ % [\s* <sep> \s*]
    }

    regex rule-component-option-value {
        <lb> \s* <rule-component-sequence> \s* <rb> || <rule-component-spec>
    }

    regex cartesian-rule-method {:i
        | <method> \s* <arrow> \s* <rule-component-option-value>
        | ':' <method> <lb> \s* <rule-component-option-value> \s* <rb>
        | '"' <method> '"' \s* ':' \s* <lb> \s* <rule-component-option-value> \s* <rb>
    }

    # This might be not a sufficiently strong production rule.
    # Cartesian rule is not allowed to have method that is a Cartesian rule.
    # (Although, in principle, that can be processed and corresponding "flattened" Cartesian rule be created.)
    rule cartesian-rule-spec {
        || <cartesian-rule>
        || <lb> <cartesian-rule> [ <sep> [ <cartesian-rule-method> | <numeric-option-spec> ]* % <sep> ]? <rb>
    }

    rule rule-spec {
        || <cartesian-rule-spec>
        || <rule-component-spec>
    }

    # Brackets
    token lb { '(' }
    token rb { ')' }

    # Arrow
    token arrow { '=>' }
    
    # Separator
    token sep { ',' }

    # Option names
    token max-points { MaxPoints | max <[\-_]> points }
    token max-recursion { MaxRecursion | max <[\-_]> recursion }
    token method {:i method }
    token min-recursion { MinRecursion | min <[\-_]> recursion }
    token points {:i points }
    token singularity-depth { SingularityDepth | singularity <[\-_]> depth }
    token singularity-handler { SingularityHandler | singularity <[\-_]> handler }
    token symbolic-processing { SymbolicProcessing | symbolic <[\-_]> processing }

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

    # Integration rule names
    token cartesian-rule { CartesianRule | cartesian <[_\-]> rule }
    token clenshaw-curtis-rule { ClenshawCurtisRule | clenshaw <[_\-]> curtis <[_\-]> rule }
    token gauss-kronrod-rule { GaussKronrodRule | gauss <[_\-]> kronrod <[_\-]> rule }
    token lobatto-kronrod-rule { LobattoKronrodRule | lobatto <[_\-]> kronrod <[_\-]> rule }
    token monte-carlo-rule { MonteCarloRule | monte <[_\-]> carlo <[_\-]> rule }
    token trapezoidal-rule { TrapezoidalRule | trapezoidal <[_\-]> rule }

    # Singularity handler names
    token imt {:i imt | iri <[\-_]> moriguti <[\-_]> takesawa }
    token no-singularity-handler {:i none | no <[\-_]> singularity <[\-_]> handler }
    token duffy-coordinates {:i DuffyCoordinates | duffy <[\-_]> coordinates } # This is, actually, a preprocessor
}
