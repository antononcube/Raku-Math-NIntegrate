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
        [ <symbol> \s* '=>' \s* <option-value> | ':' <number> <symbol> ] <!{ $<symbol>.Str.lc eq 'method' }>
    }

    rule option-sequence { <option>+ % ',' }

    rule method-option {:i
     | 'method' '=>' <method-option-value>
     | ':method(' <method-option-value> ')'
    }

    rule method-option-value {
        |  <method-symbol> 
        | '(' <method-symbol> ')'
        | '(' <method-symbol> ',' [ <method-option> | <option> ]+ % ',' ')'
    }

    rule method-symbol {
        <strategy-symbol> | <decorator-strategy-symbol> | <rule-symbol>
    }

    token global-adaptive { GlobalAdaptive | global <[_\-]> adaptive }
    token local-adaptive { LocalAdaptive | local <[_\-]> adaptive }
    token double-exponential { DoubleExponential | double <[_\-]> exponential }
    token monte-carlo { MonteCarlo | monte <[_\-]> carlo }
    token adaptive-monte-carlo { AdaptiveMonteCarlo | adaptive <[_\-]> monte <[_\-]> carlo }
    token quasi-monte-carlo { QuasiMonteCarlo | quasi <[_\-]> monte <[_\-]> carlo }
    token adaptive-quasi-monte-carlo { AdaptiveQuasiMonteCarlo | adaptive <[_\-]> quasi <[_\-]> monte <[_\-]> carlo }

    token symbolic-piecewise-subdivision { SymbolicPiecewiseSubdivision | symbolic <[_\-]> piecewise <[_\-]> subdivision }
    token even-odd-subdivision { EvenOddSubdivision | even <[_\-]> odd <[_\-]> subdivision }

    token trapezoidal-rule { TrapezoidalRule | trapezoidal <[_\-]> rule }
    token clenshaw-curtis-rule { ClenshawCurtisRule | clenshaw <[_\-]> curtis <[_\-]> rule }
    token gauss-kronrod-rule { GaussKronrodRule | gauss <[_\-]> kronrod <[_\-]> rule }
    token lobatto-kronrod-rule { LobattoKronrodRule | lobatto <[_\-]> kronrod <[_\-]> rule }
    token cartesian-rule { CartesianRule | cartesian <[_\-]> rule }
    token monte-carlo-rule { MonteCarloRule | monte <[_\-]> carlo <[_\-]> rule }

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
}
