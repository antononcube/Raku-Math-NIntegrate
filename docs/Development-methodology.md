# Development methodology for "Math::NIntegrate"

## Introduction

This document outlines the development methodology for the Raku package "Math::NIntegrate".

----

## Motivation

In order to better understand some of the methodological points is helpful to list a few motivational reasons
for how the decision to develop "Math::NIntegrate" was made.

*TBF...*

----

## Methodology

### Done by *human mind* (or "manually")

- The Object-Oriented Design (OOP) of numerical integration framework
- The complete dynamics of the interaction between regions, rules, and variable transformers
- The method parsing grammar
  - The grammar for parsing the method option has several unique properties
  - A similar grammar in BNF was originally developed for WL's `NIntegrate` 
  - Framework's grammar is written to be general and executable by Raku
  - It handles several different style of method specification: JSON-strings, WL-code-string, Raku data structures
  - Tests fo"Sentences" par

### Done via LLM / AI-agent support

- Many/most of the integration rules are well known for decades and centuries
  - Hence, LLMs / AI-agents know them well, and can reliably generate computational code for them
- Reprogramming of Fortran codes from articles
  - Some numerical integration articles have corresponding Fortran code text files available on Web
- Grammar actions
  - For a well written grammar with tests simple interpreters are easy to generate via LLMs / AI-agents

----

## Concrete steps

Here are the macro/bird view steps:

- Read the WL `NIntegrate` function page and browse `NIntegrate`'s advanced documentation.
- Program the integration rules class hierarchy and program a few rules (Gauss-Kronrod and Clenshaw-Curtis).
- Come up with Region-Rule-VariableTransformer design.
- Program region splitting and test it.
- Implement the Builder class in order to simplify the tests.
- Program Cartesian Rule in order to verify that regions can handle nD integration.
- Read / review the WL implementations of GlobalAdpative and implement it in Raku.
- Program the method option grammar and parser.
  - Come up with a comprehensive set of tests.
- Implement in the Builder class methods for making integration strategies and rules.
- Test GlobalAaptive with Cartesian rule.
- Find, program, and use the 1D "good integrator" tests by Kahaner.
- Implement fully symmetric integration rules.
- Hook-up and test the MultiDimensional integration rule
