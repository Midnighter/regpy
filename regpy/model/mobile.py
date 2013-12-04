#!/usr/bin/env python
# -*- coding: utf-8 -*-


"""
===============
Mobile Elements
===============

:Authors:
    Moritz Emanuel Beber
:Date:
    2011-10-05
:Copyright:
    Copyright(c) 2011 Jacobs University of Bremen. All rights reserved.
:File:
    mobile.py
"""


import numpy
#import random
import logging

from .misc import NullHandler, ModelParameters

logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)
logger.addHandler(NullHandler())


parameters = ModelParameters()


class BaseProduct(object):
    """
    """

    _counter = 1
    _memory = dict()

    def __new__(cls, name=u"", *args, **kw_args):
        """
        Ensures the unique instance policy of all gene products.
        """
        return cls._memory.get((cls, name), object.__new__(cls))

    def __init__(self, name=u"", *args, **kw_args):
        """
        Parameters
        ----------
        name: str (optional)
            A string uniquely identifying this element among its class.
        """
        if self.__class__._memory.has_key((self.__class__, name)):
            return
        object.__init__(self)
        self._index = self.__class__._counter
        self.__class__._counter += 1
        if name:
            self._name = name
        else:
            self._name = u"%s_%d" % (self.__class__.__name__, self._index)
        self.location = 0
        self.degradation_constant = 0
        self.__class__._memory[(self.__class__, self._name)] = self

    def __str__(self):
        """
        """
        return self._name

    def __repr__(self):
        return u"<%s.%s, %d>" % (self.__module__, self.__class__.__name__, self._index)

    def degrade(self, concentration):
        """
        """
        return numpy.random.binomial(concentration, self.degradation_constant)\
                if concentration > 0.0 else 0.0


class Enzyme(BaseProduct):

    def __init__(self, name=u"", *args, **kw_args):
        if self.__class__._memory.has_key((self.__class__, name)):
            return
        BaseProduct.__init__(self, name, *args, **kw_args)
        self.diffusion_constant = parameters.mobile.tf.diffusion()
        self.degradation_constant = parameters.mobile.tf.degradation()


class TranscriptionFactor(BaseProduct):

    def __init__(self, name=u"", *args, **kw_args):
        if self.__class__._memory.has_key((self.__class__, name)):
            return
        BaseProduct.__init__(self, name, *args, **kw_args)
        self.association_constant = parameters.mobile.tf.association()
        self.dissociation_constant = parameters.mobile.tf.dissociation()
        self.diffusion_constant = parameters.mobile.tf.diffusion()
        self.degradation_constant = parameters.mobile.tf.degradation()


class NucleoidAssociatedProtein(BaseProduct):

    def __init__(self, name=u"", *args, **kw_args):
        if self.__class__._memory.has_key((self.__class__, name)):
            return
        BaseProduct.__init__(self, name, *args, **kw_args)
        self.association_constant = parameters.mobile.nap.association()
        self.dissociation_constant = parameters.mobile.nap.dissociation()
        self.diffusion_constant = parameters.mobile.nap.diffusion()
        self.degradation_constant = parameters.mobile.nap.degradation()


class RNAPolymerase(BaseProduct):

    def __init__(self, name=u"", *args, **kw_args):
        if self.__class__._memory.has_key((self.__class__, name)):
            return
        BaseProduct.__init__(self, name, *args, **kw_args)
        self.association_constant = parameters.mobile.rnap.association()
        self.dissociation_constant = parameters.mobile.rnap.dissociation()
        self.diffusion_constant = parameters.mobile.rnap.diffusion()
        self.degradation_constant = parameters.mobile.rnap.degradation()
        self.bound = False

