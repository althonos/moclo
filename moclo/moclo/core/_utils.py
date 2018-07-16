# coding: utf-8
from __future__ import absolute_import
from __future__ import unicode_literals


def cutter_check(cutter, name):
    if cutter is NotImplemented:
        raise NotImplementedError('{} does not declare a cutter'.format(name))
    elif cutter.is_blunt():
        raise ValueError('cannot use a blunt cutter for Golden Gate')
    elif cutter.is_unknown():
        raise ValueError('cannot use an unknown cutter for Golden Gate')
