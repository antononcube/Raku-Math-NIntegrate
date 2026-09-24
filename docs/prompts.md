
----

## Grammar actions

```text
Examine "./lib/Math/NIntegrate/Processing/Grammar.rakumod" and fill-in the corresponding actions in 
"./lib/Math/NIntegrate/Processing/Actions/MethodSpec.rakumod". 
I want those actions to produce a Raku nested hashmap. The integration rules should produce a hashmap like: 
`{name => 'GaussKronrodRule', type => 'rule', points => 4}`. 
Integration strategies like: 
`{name => 'GlobalAdaptive', type => 'strategy', method => <integration-rule>, max-recursion => 12}`.
```

Several changes and completions had to be done.