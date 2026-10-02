use v6.d;

use Math::NIntegrate::Strategy;

class Math::NIntegrate::Strategy::Decorator
        is Math::NIntegrate::Strategy {

    # Should the component by rw?
    # It definitely should not be private since,
    # transparency and introspection of the integrator objects is one of framework's goals.
    # A strategy decorator cannot be an "empty" decorators -- always over a component.
    has Math::NIntegrate::Strategy $.component is required is rw;

    multi method new(Math::NIntegrate::Strategy:D $component) {
        self.bless(:$component)
    }

    method algorithm(*%args -->Map:D) {
        self.component.algorithm(|%args)
    }

    #======================================================
    # Representation
    #======================================================
    multi method Hash(::?CLASS:D:-->Map:D) {
        my %h =
                name => self.^name,
                component => self.component.Hash;
        return %h
    }

    multi method gist(::?CLASS:D:-->Str) {
        self.^name ~ '(component => ' ~ self.component.gist ~ ')'
    }
}
