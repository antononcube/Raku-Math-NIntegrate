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