# coding: utf-8
from __future__ import absolute_import
from __future__ import unicode_literals

import functools
import inspect
import warnings


class classproperty(object):
    """A class `property` decorator.
    """

    def __init__(self, getter):
        self.getter = getter

    def __get__(self, instance, owner):
        return self.getter(owner)


def isabstract(cls):
    return inspect.isabstract(cls) or any(
        getattr(cls, attr, None) is NotImplemented for attr in dir(cls)
    )


def catch_warnings(action, category=Warning, lineno=0, append=False):
    """Wrap the function in a `warnings.catch_warnings` context.

    It can be used to silence some particular specific warnings, or instead
    to treat them as errors within the function body.

    Example:
        >>> import warnings
        >>> from moclo.utils import catch_warnings
        >>> @catch_warnings('ignore')
        ... def are_you_scared():
        ...     warnings.warn("I'm warning you !")
        ...     return False
        >>> are_you_scared()
        False

    """

    def decorator(func):
        @functools.wraps(func)
        def newfunc(*args, **kwargs):
            with warnings.catch_warnings():
                warnings.simplefilter(action, category, lineno, append)
                return func(*args, **kwargs)

        return newfunc

    return decorator
