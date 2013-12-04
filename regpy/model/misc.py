#!/usr/bin/env python
# -*- coding: utf-8 -*-


"""
=============
Miscellaneous
=============

:Authors:
    Moritz Emanuel Beber
:Date:
    2011-10-05
:Copyright:
    Copyright(c) 2011 Jacobs University of Bremen. All rights reserved.
:File:
    misc.py
"""


import logging
import numpy

from meb.utils.classes import BasicOptionsManager


class NullHandler(logging.Handler):
    """
    A stub logging handler that swallows all messages. The default for all
    loggers in this package.
    """

    def emit(self, record):
        """
        Swallows messages intended for logging.
        """
        pass


class ModelParameters(BasicOptionsManager):
    """
    Singleton container class for some global model parameters.
    """

    def __init__(self, *args, **kw_args):
        """
        """
        BasicOptionsManager.__init__(self, *args, **kw_args)
        self.sequence = SequenceManager()
        self.mobile = MobileManager()
        # rng
        self.rnd_float = numpy.random.random_sample
        self.rnd_int = numpy.random.random_integers


class SequenceManager(BasicOptionsManager):
    """
    Singleton class for all sequence elements related parameters.
    """

    def __init__(self, *args, **kw_args):
        """
        """
        BasicOptionsManager.__init__(self, *args, **kw_args)
        self.gene = GeneSequenceManager()
        self.nap = NAPSequenceManager()
        self.tf = TFSequenceManager()
        self.empty = EmptySequenceManager()


class DefaultSequenceManager(BasicOptionsManager):
    """
    Singleton class for all default sequence elements related parameters.
    """

    def __init__(self, *args, **kw_args):
        """
        Attributes contain default values for all sequence elements.
        """
        BasicOptionsManager.__init__(self, *args, **kw_args)
        self.length = lambda : int(numpy.floor((numpy.random.exponential(1.0) + 1.0) * 1000.0))


class GeneSequenceManager(DefaultSequenceManager):
    """
    Singleton class for all gene sequence elements related parameters.
    """

    def __init__(self, *args, **kw_args):
        """
        Attributes inherited from DefaultSequenceManager can be replaced here or after
        instanciation to model behaviour of Genes differently.
        """
        DefaultSequenceManager.__init__(self, *args, **kw_args)
        self.production = lambda : 1.0
        self.leakage = lambda : 0.0


class NAPSequenceManager(DefaultSequenceManager):
    """
    Singleton class for all NAP binding site sequence elements related parameters.
    """

    def __init__(self, *args, **kw_args):
        """
        Attributes inherited from DefaultSequenceManager can be replaced here or after
        instanciation to model behaviour of NAP binding sites differently.

        Attributes
        ----------
        num: method
            Call to this function returns the number of NAP binding sites per
            NAP. Currently follows a binomial distribution taking into account
            the mean and prob.
        """
        DefaultSequenceManager.__init__(self, *args, **kw_args)
        self.length = lambda : 10
        self.mean = 0
        self.prob = 0.2
        self.num = lambda : numpy.random.binomial(self.mean / self.prob, self.prob)
        self.states = lambda : 5


class TFSequenceManager(DefaultSequenceManager):
    """
    Singleton class for all TF binding site sequence elements related parameters.
    """

    def __init__(self, *args, **kw_args):
        """
        Attributes inherited from DefaultSequenceManager can be replaced here or after
        instanciation to model behaviour of TF binding sites differently.
        """
        DefaultSequenceManager.__init__(self, *args, **kw_args)
        self.length = lambda : 10
        self.threshold = 1.0


class EmptySequenceManager(DefaultSequenceManager):
    """
    Singleton class for all empty sequence elements related parameters.
    """

    def __init__(self, *args, **kw_args):
        """
        Attributes inherited from DefaultSequenceManager can be replaced here or after
        instanciation to model behaviour of empty sites differently.
        """
        DefaultSequenceManager.__init__(self, *args, **kw_args)


class MobileManager(BasicOptionsManager):
    """
    Singleton class for all mobile elements related parameters.
    """

    def __init__(self, *args, **kw_args):
        """
        """
        BasicOptionsManager.__init__(self, *args, **kw_args)
        self.enzyme = EnzymeMobileManager()
        self.nap = NAPMobileManager()
        self.tf = TFMobileManager()
        self.rnap = RNAPMobileManager()


class DefaultMobileManager(BasicOptionsManager):
    """
    Default singleton class for all mobile elements related parameters.
    """

    def __init__(self, *args, **kw_args):
        """
        Attributes reflect default state for all mobile elements.
        """
        BasicOptionsManager.__init__(self, *args, **kw_args)
        self.diffusion = lambda : float(1E+04)
        self.association = lambda : 0.5
        self.dissociation = lambda : 0.2
        self.degradation = lambda : 1.0


class EnzymeMobileManager(DefaultMobileManager):
    """
    Singleton class for all enzyme products related parameters.
    """

    def __init__(self, *args, **kw_args):
        """
        Attributes inherited from DefaultMobileManager can be replaced here or after
        instanciation to model behaviour of NAPs differently.
        """
        DefaultMobileManager.__init__(self, *args, **kw_args)


class NAPMobileManager(DefaultMobileManager):
    """
    Singleton class for all mobile NAP elements related parameters.
    """

    def __init__(self, *args, **kw_args):
        """
        Attributes inherited from DefaultMobileManager can be replaced here or after
        instanciation to model behaviour of NAPs differently.
        """
        DefaultMobileManager.__init__(self, *args, **kw_args)


class TFMobileManager(DefaultMobileManager):
    """
    Singleton class for all mobile TF elements related parameters.
    """

    def __init__(self, *args, **kw_args):
        """
        Attributes inherited from DefaultMobileManager can be replaced here or after
        instanciation to model behaviour of TFs differently.
        """
        DefaultMobileManager.__init__(self, *args, **kw_args)


class RNAPMobileManager(DefaultMobileManager):
    """
    Singleton class for all mobile RNAP elements related parameters.
    """

    def __init__(self, *args, **kw_args):
        """
        Attributes inherited from DefaultMobileManager can be replaced here or after
        instanciation to model behaviour of RNAPs differently.
        """
        DefaultMobileManager.__init__(self, *args, **kw_args)

