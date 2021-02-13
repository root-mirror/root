""" Pythonization API.
"""

__all__ = [
    'add_pythonization',
    'remove_pythonization',
    'pin_type',
    ]

def _set_backend(backend):
    global _backend
    _backend = backend


# user-provided, general pythonizations
def add_pythonization(pythonizor, scope = ''):
    """<pythonizor> should be a callable taking two arguments: a class proxy,
    and its C++ name. It is called each time a named class from <scope> (the
    global one by default, but a relevant C++ namespace is recommended) is bound.
    """
    return _backend.add_pythonization(pythonizor, scope)

def remove_pythonization(pythonizor, scope = ''):
    """Remove previously registered <pythonizor> from <scope>.
    """
    return _backend.remove_pythonization(pythonizor, scope)


# prevent auto-casting (e.g. for interfaces)
def pin_type(klass):
    return _backend._pin_type(klass)

# exception pythonizations
def add_exception_mapping(cpp_exception, py_exception):
    _backend.UserExceptions[cpp_exception] = py_exception



#--- Pythonization factories --------------------------------------------

def set_gil_policy(match_class, match_method, release_gil=True):
    return set_method_property(match_class, match_method, '__release_gil__', int(release_gil))


def set_ownership_policy(match_class, match_method, python_owns_result):
    return set_method_property(match_class, match_method, 
                               '__creates__', int(python_owns_result))


# NB: Ideally, we'd use the version commented out below, but for now, we
#     make do with the hackier version here.
def rename_attribute(match_class, orig_attribute, new_attribute, keep_orig=False):
    class attribute_pythonizor(object):
        class getter(object):
            def __init__(self, attr):
                self.attr = attr
            def __call__(self, obj):
                return getattr(obj, self.attr)

        class setter(object):
            def __init__(self, attr):
                self.attr = attr
            def __call__(self, obj, value):
                return setattr(obj, self.attr, value)

        class deleter(object):
            def __init__(self, attr):
                self.attr = attr
            def __call__(self, obj):
                return delattr(obj, self.attr)

        def __init__(self, match_class, orig_attribute, new_attribute, keep_orig):
            import re
            self.match_class = re.compile(match_class)
            self.match_attr = re.compile(orig_attribute)
            self.new_attr = new_attribute
            self.keep_orig = keep_orig

        def __call__(self, obj, name):
            if not self.match_class.match(name):
                return
            for k in dir(obj): #.__dict__:
                if self.match_attr.match(k):
                    tmp = property(self.getter(k), self.setter(k), self.deleter(k))
                    setattr(obj, self.new_attr, tmp)
                    #if not self.keep_orig: delattr(obj, k)
    return attribute_pythonizor(match_class, orig_attribute, new_attribute, keep_orig)

# def rename_attribute(match_class, orig_attribute, new_attribute, keep_orig=False):
#     class method_pythonizor:
#         def __init__(self, match_class, orig_attribute, new_attribute, keep_orig):
#             import re
#             self.match_class = re.compile(match_class)
#             self.match_attr = re.compile(orig_attribute)
#             self.new_attr = new_attribute
#             self.keep_orig = keep_orig

#         def __call__(self, obj, name):
#             import sys
#             if not self.match_class.match(name):
#                 return
#             sys.stderr.write("%s %s %s %s" % ("!!!", obj, name, "\n"))
#             for k in dir(obj): #obj.__dict__:
#                 if not self.match_attr.match(k): continue
#                 try:
#                    tmp = getattr(obj, k)
#                 except Exception as e:
#                    continue
#                 setattr(obj, self.new_attr, tmp)
#                 if not self.keep_orig: delattr(obj, k)
#     return method_pythonizor(match_class, orig_attribute, new_attribute, keep_orig)


# Shared with PyPy:

def add_overload(match_class, match_method, overload):
    class method_pythonizor(object):
        def __init__(self, match_class, match_method, overload):
            import re
            self.match_class = re.compile(match_class)
            self.match_method = re.compile(match_method)
            self.overload = overload

        def __call__(self, obj, name):
            if not self.match_class.match(name):
                return
            for k in dir(obj): #.__dict__:
               try:
                   tmp = getattr(obj, k)
               except:
                   continue
               if self.match_method.match(k):
                   try:
                       tmp.__add_overload__(overload)
                   except AttributeError: pass
    return method_pythonizor(match_class, match_method, overload)


def compose_method(match_class, match_method, g):
    class composition_pythonizor(object):
        def __init__(self, match_class, match_method, g):
            import re
            self.match_class = re.compile(match_class)
            self.match_method = re.compile(match_method)
            self.g = g

        def __call__(self, obj, name):
            if not self.match_class.match(name):
                return
            g = self.g
            for k in obj.__dict__:
                if not self.match_method.match(k):
                    continue
                try:
                    f = getattr(obj, k)
                except:
                    continue
                def make_fun(f, g):
                    def h(self, *args, **kwargs):
                        return g(self, f(self, *args, **kwargs))
                    return h
                h = make_fun(f, g)
                setattr(obj, k, h)
    return composition_pythonizor(match_class, match_method, g)


def set_method_property(match_class, match_method, prop, value):
    class method_pythonizor(object):
        def __init__(self, match_class, match_method, prop, value):
            import re
            self.match_class = re.compile(match_class)
            self.match_method = re.compile(match_method)
            self.prop = prop
            self.value = value

        def __call__(self, obj, name):
            if not self.match_class.match(name):
                return
            for k in dir(obj): #.__dict__:
                try:
                    tmp = getattr(obj, k)
                except:
                    continue
                if self.match_method.match(k):
                    setattr(tmp, self.prop, self.value)
    return method_pythonizor(match_class, match_method, prop, value)


def make_property(match_class, match_get, match_set=None, match_del=None, prop_name=None):
    class property_pythonizor(object):
        def __init__(self, match_class, match_get, match_set, match_del, prop_name):
            import re
            self.match_class = re.compile(match_class)

            self.match_get = re.compile(match_get)
            match_many_getters = self.match_get.groups == 1

            if match_set:
                self.match_set = re.compile(match_set)
                match_many_setters = self.match_set.groups == 1
                if match_many_getters ^ match_many_setters:
                    raise ValueError('Must match getters and setters equally')
            else:
                self.match_set = None

            if match_del:
                self.match_del = re.compile(match_del)
                match_many_deleters = self.match_del.groups == 1
                if match_many_getters ^ match_many_deleters:
                    raise ValueError('Must match getters and deleters equally')
            else:
                self.match_del = None

            self.match_many = match_many_getters
            if not (self.match_many or prop_name):
                raise ValueError("If not matching properties by regex, need a property name with exactly one substitution field")
            if self.match_many and prop_name:
                if prop_name.format(').!:(') == prop_name:
                    raise ValueError("If matching properties by regex and providing a property name, the name needs exactly one substitution field")

            self.prop_name = prop_name

        def make_get_del_proxy(self, getter):
            class proxy(object):
                def __init__(self, getter):
                    self.getter = getter

                def __call__(self, obj):
                    return getattr(obj, self.getter)()
            return proxy(getter)

        def make_set_proxy(self, setter):
            class proxy(object):
                def __init__(self, setter):
                    self.setter = setter

                def __call__(self, obj, arg):
                    return getattr(obj, self.setter)(arg)
            return proxy(setter)

        def __call__(self, obj, name):
            if not self.match_class.match(name):
                return

            names = []
            named_getters = {}
            named_setters = {}
            named_deleters = {}

            if not self.match_many:
                fget, fset, fdel = None, None, None

            for k in dir(obj): #.__dict__:
                match = self.match_get.match(k)
                try:
                    tmp = getattr(obj, k)
                except:
                    continue
                if match and hasattr(tmp, '__call__'):
                    if self.match_many:
                        name = match.group(1)
                        named_getters[name] = k
                    else:
                        fget = self.make_get_del_proxy(k)
                        break

            if self.match_set:
                for k in dir(obj): #.__dict__:
                    match = self.match_set.match(k)
                    try:
                        tmp = getattr(obj, k)
                    except:
                        continue
                    if match and hasattr(tmp, '__call__'):
                        if self.match_many:
                            name = match.group(1)
                            named_setters[name] = k
                        else:
                            fset = self.make_set_proxy(k)
                        break

            if self.match_del:
                for k in dir(obj): #.__dict__:
                    match = self.match_del.match(k)
                    try:
                        tmp = getattr(obj, k)
                    except:
                        continue
                    if match and hasattr(tmp, '__call__'):
                        if self.match_many:
                            name = match.group(1)
                            named_deleters[name] = k
                        else:
                            fdel = self.make_get_del_proxy(k)
                            break

            if not self.match_many:
                new_prop = property(fget, fset, fdel)
                setattr(obj, self.prop_name, new_prop)
                return

            names += list(named_getters.keys())
            names += list(named_setters.keys())
            names += list(named_deleters.keys())
            names = set(names)

            properties = []
            for name in names:
                if name in named_getters:
                    fget = self.make_get_del_proxy(named_getters[name])
                else:
                    fget = None

                if name in named_setters:
                    fset = self.make_set_proxy(named_setters[name])
                else:
                    fset = None

                if name in named_deleters:
                    fdel = self.make_get_del_proxy(named_deleters[name])
                else:
                    fdel = None

                new_prop = property(fget, fset, fdel)
                if self.prop_name:
                    prop_name = self.prop_name.format(name)
                else:
                    prop_name = name

                setattr(obj, prop_name, new_prop)

    return property_pythonizor(match_class, match_get, match_set, match_del, prop_name)

