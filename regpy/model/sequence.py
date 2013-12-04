#!/usr/bin/env python
# -*- coding: utf-8 -*-


"""
=================
Sequence Elements
=================

:Authors:
    Moritz Emanuel Beber
:Date:
    2011-10-05
:Copyright:
    Copyright(c) 2011 Jacobs University of Bremen. All rights reserved.
:File:
    sequence.py
"""


import copy
import operator
import numpy
import logging
import networkx as nx

from . import mobile
from .misc import NullHandler, ModelParameters


logger = logging.getLogger(__name__)
logger.addHandler(NullHandler())


parameters = ModelParameters()


class SequenceElement(object):
    """
    """

    _counter = 1
    _memory = dict()

    def __new__(cls, name=u"", *args, **kw_args):
        """
        Ensures the unique instance policy of all sequence sites.
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
        self._length = 0
        self._symbol = u""
        self.occupied = False
        self.__class__._memory[(self.__class__, self._name)] = self

    def __len__(self):
        """
        """
        return self._length

    def __str__(self):
        """
        """
        return u"%s-%d" % (self._symbol, self._index)

    def __repr__(self):
        return u"<%s.%s, %d>" % (self.__module__, self.__class__.__name__, self._index)

    def reset(self):
        """
        """
        self.occupied = False


class EmptySite(SequenceElement):
    """
    """

    def __init__(self, name="", *args, **kw_args):
        """
        """
        if self.__class__._memory.has_key((self.__class__, name)):
            return
        SequenceElement.__init__(self, name, *args, **kw_args)
        self._length = parameters.sequence.empty.length()
        self._symbol = u"E"


class GeneSite(SequenceElement):
    """
    """

    def __init__(self, name=u"", promoters=False, product=None, *args, **kw_args):
        """
        """
        if self.__class__._memory.has_key((self.__class__, name)):
            return
        SequenceElement.__init__(self, name, *args, **kw_args)
        self._length = parameters.sequence.gene.length()
        self._symbol = u"G"
        self.promoters = promoters if promoters else list()
        self.product = product
        self.rate = parameters.sequence.gene.production()
        self._active = False

    def __len__(self):
        """
        """
        return self._length + sum(len(elem) for elem in self.promoters)

    def __str__(self):
        """
        """
        promoter = u"|".join(str(r) for r in self.promoters)
        ss = u"%s%s%s" % (promoter, u"|", super(GeneSite, self).__str__()) if promoter else super(GeneSite, self).__str__()
        return ss

    def is_active(self):
        """
        """
        states = [tf_site.regulation if tf_site.bound else 0\
                for tf_site in self.promoters]
        if not states:
            self.rate = parameters.sequence.gene.leakage()
            return True
#        state = reduce(operator.add, states)
        state = reduce(operator.mul, states)
        if state > 0:
            self.rate = parameters.sequence.gene.production()
            self._active = True
        elif state < 0:
            self._active = False
        else:
            self.rate = parameters.sequence.gene.leakage()
            self._active = True
        return self._active

    def reset(self):
        """
        """
        super(GeneSite, self).reset()
        self._active = False
        for tf_site in self.promoters:
            tf_site.reset()


class BindingSite(SequenceElement):
    """
    """

    def __new__(cls, ligand, regulation, name=u"", *args, **kw_args):
        """
        Ensures the unique instance policy of all ligand binding sites.
        """
        return cls._memory.get((cls, name), SequenceElement.__new__(cls,
                name, *args, **kw_args))

    def __init__(self, ligand, regulation, name=u"", *args, **kw_args):
        """
        Parameters
        ----------
        regulation: int
            type of regulation how a bound ligand affects the site
        """
        if self.__class__._memory.has_key((self.__class__, name)):
            return
        SequenceElement.__init__(self, name, *args, **kw_args)
        self.regulation = int(regulation)
        self.ligand = ligand
        self.distance = 0
        self.factor = 0.0
        self.bound = False

    def update_distance(self, index, sequence):
        if index > self.ligand.location:
            sequence = sequence[self.ligand.location:index]
        else:
            sequence = sequence[index:self.ligand.location]
        self.distance = sum(len(site) for site in sequence)

    def reset(self):
        """
        """
        super(BindingSite, self).reset()
        self.bound = False


class TFBindingSite(BindingSite):
    """
    """

    def __init__(self, ligand, regulation, name=u"", *args, **kw_args):
        """
        Parameters
        ----------
        regulation: int
            -1 inhibitory
             0 neutral
             1 activating

        """
        if self.__class__._memory.has_key((self.__class__, name)):
            return
        BindingSite.__init__(self, ligand, regulation, name, *args, **kw_args)
        self._length = parameters.sequence.tf.length()
        self._symbol = u"T"


class NAPBindingSite(BindingSite):
    """
    """

    def __init__(self, ligand, regulation, name=u"", *args, **kw_args):
        """
        Parameters
        ----------
        regulation: int
            In combination with the super-coiling state of the sequence bound
            NAPBingindSites modify the sequence accessibility.
            -2 strongly restrictive
            -1 restrictive
             0 neutral
             1 enhancing
             2 strongly enhancing
        """
        if self.__class__._memory.has_key((self.__class__, name)):
            return
        BindingSite.__init__(self, ligand, regulation, name, *args, **kw_args)
        self._length = parameters.sequence.nap.length()
        self._symbol = u"N"


class Sequence(list):
    """
    """

    def __init__(self, *args, **kw_args):
        """
        Parameters
        ----------
        """
        list.__init__(self, *args, **kw_args)
        self.concentrations = dict()
        self.polymerases = dict()

    def __str__(self):
        return u"|".join(str(item) for item in self)

    def linearise_trn(self, trn):
        for gene in trn:
            logger.debug("%s product: %s", repr(gene), repr(gene.product))
            logger.debug("\t%s", trn.pred[gene])
            for (regulator, data) in trn.pred[gene].iteritems():
                tf_site = TFBindingSite(ligand=regulator.product,
                        regulation=data.get("regulation", 0))
                logger.debug("\t%s added", repr(tf_site))
                gene.promoters.append(tf_site)
            self.append(gene)
            logger.debug("%s", str(self))

    def initialise_promoters(self, genes):
        """
        """
        for site in genes:
            for tf_site in site.promoters:
                while not tf_site.regulation:
                    if parameters.rnd_float() < 0.5:
                        tf_site.regulation = 1
                    else:
                        tf_site.regulation = -1
            # binding site location plays a small role in diffusion
#            random.shuffle(site.promoters) # want something deterministic here

    def initialise_naps(self, genes):
        """
        """
        num_states = parameters.sequence.nap.states()
        nap_pdf = [float(i) / num_states for i in range(1, num_states + 1)]
        for site in genes:
            nap = mobile.NucleoidAssociatedProtein()
            site.product = nap
            num_sites = parameters.sequence.nap.num()
            while not num_sites:
                num_sites = parameters.sequence.nap.num()
            for i in range(num_sites):
                regulation = 0
                while not regulation:
                    prob = parameters.rnd_float()
                    for (i, val) in enumerate(nap_pdf):
                        if prob <= val:
                            break
                    regulation = i - 2
                nap_site = NAPBindingSite(ligand=nap, regulation=regulation)
                self.append(nap_site)

    def initialise(self):
        """
        """
        for (i, site) in enumerate(self):
            if isinstance(site, BindingSite):
                site.update_distance(i, self)
                # diffusion factor
                site.factor = site.ligand.association_constant *\
                        numpy.exp(-site.distance/site.ligand.diffusion_constant)
            elif isinstance(site, GeneSite):
                for tf_site in site.promoters:
                    tf_site.update_distance(i, self)
                    # diffusion factor
                    tf_site.factor = tf_site.ligand.association_constant *\
                            numpy.exp(-tf_site.distance / tf_site.ligand.diffusion_constant)
                    logger.debug("%s constant binding factor = %f",
                            str(tf_site), tf_site.factor)

    def next(self):
        """
        """
        old = copy.copy(self.concentrations)
        # update sequence elements
        for site in self:
#            logger.debug(str(site))
            # update TFs in promoter regions
            if isinstance(site, GeneSite):
                for tf_site in site.promoters:
                    conc = old.get(tf_site.ligand, 0.0)
                    if tf_site.factor * conc >= parameters.sequence.tf.threshold:
                        tf_site.bound = True
                    else:
                        tf_site.bound = False
##                    logger.debug("[%s] = %f @ %s factor = %f:",
##                            str(tf_site.ligand), conc, str(tf_site), tf_site.factor)
#                    if tf_site.bound and\
#                            parameters.rnd_float() < tf_site.ligand.dissociation_constant:
##                        logger.debug("\treleased")
#                        tf_site.bound = False
#                    elif conc > 0.0 and parameters.rnd_float() < tf_site.factor * conc:
##                        logger.debug("\tbound")
#                        tf_site.bound = True
            # update NAP binding sites
            elif isinstance(site, NAPBindingSite):
                pass
        # update polymerases
#        rm = set()
#        last = len(self) - 1
#        for (rnap, pos) in self.polymerases.iteritems():
#            if pos >= len(self):
#                rm.add(rnap)
#                continue
#            site = self[pos]
#            logger.debug("%s @ %s:", str(rnap), str(site))
#            if rnap.bound:
#                # bound polymerase can dissociate or express
#                if parameters.rnd_float() < rnap.dissociation_constant:
#                    logger.debug("\treleased")
#                    rnap.bound = False
#                else:
#                    logger.debug("\ttranscribed %s", str(site.product))
#                    self.concentrations[site.product] =\
#                            self.concentrations.get(site.product, 0.0) + site.rate
#            elif isinstance(site, GeneSite) and site.is_active() and\
#                    parameters.rnd_float() < rnap.association_constant:
#                # unbound polymerase can bind
#                logger.debug("\tbound")
#                rnap.bound = True
#            elif pos == last:
#                logger.debug("\tleft")
#                self.polymerases[rnap] += 1
#                self[pos].occupied = False
#            elif not self[pos + 1].occupied:
#                # unbound polymerase moves on
#                logger.debug("\tmoved on")
#                self.polymerases[rnap] += 1
#                self[pos].occupied = False
#                self[pos + 1].occupied = True
#        for rnap in rm:
#            del self.polymerases[rnap]
        rm = set()
        last = len(self) - 1
        for (rnap, pos) in self.polymerases.iteritems():
            if pos >= len(self):
                rm.add(rnap)
                continue
            site = self[pos]
            logger.debug("%s @ %s:", str(rnap), str(site))
            if rnap.bound:
                logger.debug("\ttranscribed %s", str(site.product))
                logger.debug(str(site.rate))
                self.concentrations[site.product] =\
                        self.concentrations.get(site.product, 0.0) + site.rate
                logger.debug("\treleased")
                rnap.bound = False
                rnap.was_bound = True
            elif isinstance(site, GeneSite) and site.is_active() and\
                    not rnap.was_bound:
                logger.debug("\tbound")
                rnap.bound = True
            elif pos == last:
                logger.debug("\tleft")
                self.polymerases[rnap] += 1
                self[pos].occupied = False
            elif not self[pos + 1].occupied:
                # unbound polymerase moves on
                logger.debug("\tmoved on")
                self.polymerases[rnap] += 1
                self[pos].occupied = False
                self[pos + 1].occupied = True
                rnap.was_bound = False
        for rnap in rm:
            del self.polymerases[rnap]
        # update concentrations
        logger.debug(str(self.concentrations))
        for (mol, conc) in self.concentrations.iteritems():
#            deg = mol.degrade(conc)
            deg = numpy.ceil(mol.degradation_constant * conc)
            logger.debug(str(deg))
            self.concentrations[mol] = max(conc - deg, 0.0)
        logger.debug(str(self.concentrations))
        del old

    def introduce_polymerase(self):
        """
        In this function and in 'next', there is no collision check between
        polymerases yet.
        """
        if not self:
            return False
        if self[0].occupied:
            return False
        rnap = mobile.RNAPolymerase()
        rnap.was_bound = False
        self.polymerases[rnap] = 0
        self[0].occupied = True
        return True

    def reset(self):
        self.polymerases = dict()
        self.concentrations = dict()
        for site in self:
            site.reset()


def network2trn(network):
    """
    """
    trn = nx.DiGraph(name="TRN")
    mapping = dict()
    for node in network:
        if network.out_degree(node) > 0:
            tf = mobile.TranscriptionFactor()
        else:
            tf = None
        gene = GeneSite(product=tf)
        trn.add_node(gene)
        mapping[node] = gene
    for (u, v, data) in network.edges_iter(data=True):
        trn.add_edge(mapping[u], mapping[v], regulation=data.get("regulation", 0))
    return trn

