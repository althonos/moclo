# coding: utf-8
"""An implementation of the Golden Braid 3.0 for the Python MoClo library.
"""
from __future__ import absolute_import
from __future__ import unicode_literals

from Bio.Restriction import BsaI, BsmBI

from ..core import parts, modules, vectors

__author__ = 'Martin Larralde <martin.larralde@ens-paris-saclay.fr>'
__version__ = (
    __import__('pkg_resources')
        .resource_string(__name__, 'gb3.version')
        .strip()
        .decode('ascii')
)


# VECTORS ####################################################################

class GB3OmegaVector(vectors.CassetteVector):
    """A GoldenBraid 3.0 Level Omega vector.
    """

    cutter = BsaI

    @classmethod
    def structure(cls):  # noqa: D105
        return (
            'CGTCTC'  # BsmBI
            'NN'
            '(NNNN)'  # Vector downstream overhang
            '(N'
            'GAGACC'  # BsaI
            'N*'      # Placeholder sequence
            'GGTCTC'  # BsaI
            'N)'
            '(NNNN)'  # Vector upstream overhang
            'NN'
            'GAGACG'  # BsmBI
        )


class GB3AlphaVector(vectors.DeviceVector):
    """A GoldenBraid 3.0 Level Alpha vector.
    """

    cutter = BsmBI

    @classmethod
    def structure(cls):  # noqa: D105
        return (
            'GGTCTC'  # BsaI
            'N'
            '(NNNN)'  # Vector downstream overhang
            '(N'
            'GAGACG'  # BsmBI
            'N*'      # Placeholder sequence
            'CGTCTC'  # BsmBI
            'N)'
            '(NNNN)'  # Vector upstream overhang
            'N'
            'GAGACC'  # BsaI
        )


# MODULES ####################################################################

class GB3OmegaModule(modules.Entry):
    """A GoldenBraid 3.0 Level Omega module.
    """

    cutter = BsaI

class GB3AlphaModule(modules.Cassette):
    """A GoldenBraid 3.0 Level Alpha module.
    """

    cutter = BsmBI

# PARTS ######################################################################

class GB3Part(parts.AbstractPart):
    """An GoldenBraid 3.0 part, also refered as *element* in the GB lexicon.

    A part is a plasmid with standardized flanking overhang sequences
    that allows immediate type recognition.
    """

    cutter = BsaI
    signature = NotImplemented

# Transcriptional Units Parts ################################################

class GB3Dist(GB3Part, GB3OmegaModule):
    signature = ('GGAG', 'TGAC')

class GB3Prox(GB3Part, GB3OmegaModule):
    signature = ('TGAC', 'TCCC')

class GB3Core(GB3Part, GB3OmegaModule):
    signature = ('TCCC', 'TACT')

class GB3UpstreamUTR(GB3Part, GB3OmegaModule):
    signature = ('TACT', 'CCAT')

class GB3NTag(GB3Part, GB3OmegaModule):
    signature = ('CCAT', 'AATG')

class GB3CodingSequence1(GB3Part, GB3OmegaModule):
    signature = ('AATG', 'AGCC')

class GB3CodingSequence2(GB3Part, GB3OmegaModule):
    signature = ('AGCC', 'TTCG')

class GB3CTag(GB3Part, GB3OmegaModule):
    signature = ('TTCG', 'GCTT')

class GB3DownstreamUTR(GB3Part, GB3OmegaModule):
    signature = ('GCTT', 'GGTA')

class GB3Terminator(GB3Part, GB3OmegaModule):
    signature = ('GGTA', 'CGCT')

# Crisprs Parts ##############################################################

class GB3PolIIIMonocotPromoter(GB3Part, GB3OmegaModule):
    signature = ('GGAG', 'GGCA')

class GB3MTarget(GB3Part, GB3OmegaModule):
    signature = ('GGCA', 'GTTT')

class GB3PolIIIDicotPromoter(GB3Part, GB3OmegaModule):
    signature = ('GGAG', 'ATTG')

class GB3DTarget(GB3Part, GB3OmegaModule):
    signature = ('ATTG', 'GTTT')

class GB3sgRNA(GB3Part, GB3OmegaModule):
    signature = ('GTTT', 'CGCT')

# Composite Parts ############################################################

class GB3InteractionAdaptor(GB3Part, GP3OmegaModule):
    signature = ('GGAG', 'AATG')
