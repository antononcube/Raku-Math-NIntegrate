# Context use in variable transformers

## Introduction

This is a "low-level" design note.

----

## Problem with the initial design

Initially, `Math::NIntegrate::VariableTransformer` had this attribute:

```raku
    #| Contextual variable transformer (e.g. Composite)
    has Math::NIntegrate::VariableTransformer $.context is rw;
```

And methods like:

```raku
    #| Get original bounds
    method get-original-bounds(-->Map:D) {
        return do if $!context {
            $!context.get-original-bounds()
        } else {
            %(min => @!min-original-bounds, max => @!max-original-bounds)
        }
    }

    #| Get transformation bounds
    method get-transform-bounds(-->Map:D) {
        return do if $!context {
            $!context.get-transform-bounds()
        } else {
            %(min => @!min-transform-bounds, max => @!max-transform-bounds)
        }
    }

    #| Get region object
    method get-region() {
        return do if $!context {
            $!context.get-region()
        } else {
            self.region
        }
    }
```

Unfortunately, Raku hangs when initializations like this a made in `TWEAK`:

```raku
        for @!stack -> $vt { $vt.context = self }
```

Or when like this in a dedicated `Math::NIntegerate::Builder` method of `Math::NIntegerate::VariableTransformer::Composite`:

```raku
    method make-variable-transformer-composite(
            :$region
            --> Math::NIntegrate::VariableTransformer::Composite) {
    
        my $obj = Math::NIntegrate::VariableTransformer::Composite.new(:$region);

        my @min-original-bounds = $region.min;
        my @max-original-bounds = $region.max;
        my $vtInf = Math::NIntegrate::VariableTransformer::Infinity.new(:@min-original-bounds, :@max-original-bounds);
        $obj.add($vtInf);
        return $obj
    }
```

Where the method `add` has the line:

```raku
        @!stack.push($obj);
        $obj.context = self;
```

**Remark:** In C and C++ implementations of that design the stack would be just pointers, so, there would be no inherent problems.

----

## Remedy

Have the argument `:$context = Nil` of the `transform` methods of the `Math::NIntegerate::VariableTransformer` classes.

That would require a more explicit management of the use of the objects `Math::NIntegerate::VariableTransformer::Compoiste` and `Math::NIntegrate::Region`.
