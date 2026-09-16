# Methodology and Design 

<p style="font-size:24px; font-weight:700;">...for the Raku package "Math::NIntegrate"</p>

---

## Introduction

The Raku numerical integration framework "Math::NIntegrate", [AAp1], is based on Object-Oriented Programming (OOP) design and implementation.
[Design Patterns (by GoF)](https://en.wikipedia.org/wiki/Design_Patterns) are extensively used.

This document outlines the OOP architecture and gives rationale for some of the less obvious design decisions.

---

## Overall process

Each integration strategy of "Math::NIntegrate" creates and manipulates a collection of integration regions. 
Each integration region can have its own integrand and/or integration rule. 
"Math::NIntegrate"'s main integration strategy "GlobalAdaptive" keeps the integration regions in a heap according to their error. 
The sum of the integral estimates of all regions make the global integral estimate. 
The sum of the integral errors make the global error. 
If the global error is larger than the desired tolerance "GlobalAdaptive" splits the region with the largest error estimate into two regions and applies the integration rule. 
If too many splittings have been done then a singularity handler is applied over the last region split.

At each step of an integration strategy the option "IntegrationMonitor" obtains as an argument the list of integration regions used in that step. 
Below is a table that shows methods that can be applied to each integration region in that list.

| Property               | Description                                                                                                          |
|------------------------|----------------------------------------------------------------------------------------------------------------------|
| `"Axis"`               | Returns the axis with maximum error.                                                                                 |
| `"Boundaries"`         | Returns the boundaries of the region.                                                                                |
| `"Dimension"`          | Returns the dimension of the region.                                                                                 |
| `"Error"`              | Returns the error estimate of the last integral estimation.                                                          |
| `"GetRule"`            | Returns the attached integration rule.                                                                               |
| `"GetSamplingPoints"`  | Returns the sampling points used for the integration estimation.                                                     |
| `"GetValues"`          | Returns the integrand values over the sampling points used for the integration estimation.                           |
| `"Integral"`           | Returns the integral estimate of the last integral estimation.                                                       |
| `"Integrand"`          | Returns the integrand (numerical function) attached to the region.                                                   |
| `"NumericalFunction"`  | Returns the integrand (numerical function) attached to the region.                                                   |
| `"OriginalBoundaries"` | Returns the original boundaries of the region (when singularity handling is applied it differs from `"Boundaries"`). |
| `"Properties"`         | Returns all methods applicable.                                                                                      |
| `"WorkingPrecision"`   | Returns the precision with which the function `N` is used during evaluation.                                         |

---

# Integration rule object

The Integration Rule (IRule) object has these three lists:

- Abscissas
- Weights
- Error weights

An IRule is attached to one or more Integration Region (IReg) objects.

IRule's method `integrate` uses the Integration Region (IReg) object method `eval-integrand`.

---

## Integration region object

As mentioned in the previous section, the Integration Region (IReg) objects have a central role in the framework.
IReg objects are instances of the class `Math::NIntegrate::Region`. Every IR has an Integration Rule (IRule) object and 
a Variable Transformer (VT) object. 

For a given IReg object, `R`:

- Its IRule object `R.rule` has a reference to `R` and uses the `R` method `R.eval-integrand` to do rule computations
- Its VT object `R.variable-transformer` has a reference to `R` and uses `R` boundaries.

---

## Variable transformer

The Variable Transformer (VT) hierarchy of classes implements the [Composite design pattern](https://en.wikipedia.org/wiki/Composite_pattern).

The fundamental classes `Math::NIntegrate::Infinity` and `Math::NIntegrate::Affine`. 
Some singularity handlers are implemented VT classes. Some integration rules must be paired with corresponding variable transformers.

```mermaid
%%{init: {"theme": "redux-dark-color"} }%%
classDiagram
    direction LR
    class Component["Math::NIntegrate::VariableTransformer"] {
        +context
        +region
        +transform()
    }
    class Composite["Math::NIntegrate::VariableTransformer::Composite"] {
        -stack: List~Component~
        -vtAffine
        +transform()
        +add(Component)
        +remove(Component)
    }
    class Affine["Math::NIntegrate::VariableTransformer::Affine"] {
        +transform()
    }
    class Infinity["Math::NIntegrate::VariableTransformer::Infinity"] {
        +transform()
    }
    class IMT["Math::NIntegrate::VariableTransformer::IMT"] {
        +transform()
    }
    class DoubleExponential["Math::NIntegrate::VariableTransformer::DoubleExponential"] {
        +transform()
    }
    class Reverse["Math::NIntegrate::VariableTransformer::Reverse"] {
        +transform()
    }

    Component <|-- Composite
    Component <|-- Affine
    Component <|-- Infinity
    Component <|-- IMT
    Component <|-- DoubleExponential
    Component <|-- Reverse
    Composite "1" *-- "0..*" Component
    Composite *-- Affine
```

---

## References

### Article, blog posts

[AAmse1] Anton Antonov,
["Determining which rule NIntegrate selects automatically", answer](https://mathematica.stackexchange.com/a/96663),
(2015)
[MathematicaStackExchange](https://mathematica.stackexchange.com).

### Wolfram Language documentation

### Raku packages

[AAp1] Anton Antonov,
[Math::NIntegrate, Raku package](https://github.com/antononcube/Raku-Data-NIntegrate),
(2021-2026),
[GitHub/antononcube](https://github.com/antononcube).

[AAp1] Anton Antonov, 
[Data::Transformers, Raku package](https://github.com/antononcube/Raku-Data-Transformers),
(2026),
[GitHub/antononcube](https://github.com/antononcube).

[AAp2] Anton Antonov,
[LeftistHeap, Raku package](https://github.com/antononcube/Raku-LeftistHeap),
(2026),
[GitHub/antononcube](https://github.com/antononcube).

### Repositories

[AAr1] Anton Antonov,
["NIntegrate, the missing manual"](https://github.com/antononcube/NIntegrateTheMissingManual-book),
(2019-2026),
[GitHub/antononcube](https://github.com/antononcube).