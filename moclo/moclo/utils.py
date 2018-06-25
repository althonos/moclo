# coding: utf-8
from __future__ import absolute_import
from __future__ import unicode_literals

import functools
import warnings


def catch_warnings(action, category=Warning, lineno=0, append=False):
    def decorator(func):
        @functools.wraps(func)
        def newfunc(*args, **kwargs):
            with warnings.catch_warnings():
                warnings.simplefilter(action, category, lineno, append)
                return func(*args, **kwargs)
        return newfunc
    return decorator
