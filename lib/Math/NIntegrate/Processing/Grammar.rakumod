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

    rule strategy-symbol {
        'GlobalAdaptive' | 'LocalAdaptive' | 'DoubleExponential'
      | 'MonteCarlo' | 'AdaptiveMonteCarlo' | 'QuasiMonteCarlo'
      | 'AdaptiveQuasiMonteCarlo' | <top-level-strategy>
    }

    rule decorator-strategy-symbol {
        'SymbolicPiecewiseSubdivision' | 'EvenOddSubdivision'
    }

    rule rule-symbol {
        'TrapezoidalRule' | 'GaussKronrodRule' | 'LobattoKronrodRule'
      | 'CartesianRule' | 'MonteCarloRule' | <top-level-rule>
    }
}
