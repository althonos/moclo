# coding: utf-8
"""An implementation of the Yeast ToolKit for the Python MoClo library.

This module is tested against the officials parts available in the Yeast
ToolKit (YTK), and also against the Pichia ToolKit (PTK) parts since they were
designed to be compatible with each other.

The documentation of this module is mostly adapted from the *Lee et al.*
supplementary data. Each item also has specific sections that are organized
as follow:

**Note**:
    this section describes a behaviour that is not part of the YTK
    standard, but that is implemnted in all YTK official parts, and encouraged
    to follow by the YTK authors.
**Caution**
    this section describes a behaviour that goes against the MoClo standard,
    but which you are entitled to follow for your parts to be valid YTK parts.
**Danger**
    this section describes a quirk specific to the ``moclo-ytk`` library.

References:
    1. `Lee, M. E., DeLoache, W. C., Cervantes, B., Dueber, J. E. (2015).
       A Highly Characterized Yeast Toolkit for Modular, Multipart Assembly.
       ACS Synthetic Biology, 4(9), 975–986.
       <https://doi.org/10.1021/sb500366v>`__
    2. `Obst, U., Lu, T. K., Sieber, V. (2017).
       A Modular Toolkit for Generating Pichia pastoris Secretion Libraries.
       ACS Synthetic Biology, 6(6), 1016–1025
       <https://doi.org/10.1021/acssynbio.6b00337>`__
    3. `Weber, E., Engler, C., Gruetzner, R., Werner, S., Marillonnet, S. (2011).
       A Modular Cloning System for Standardized Assembly of Multigene Constructs.
       PLOS ONE, 6(2), e16765.
       <https://doi.org/10.1371/journal.pone.0016765>`__

"""

import abc

import six
from Bio.Restriction import BsaI, BsmBI

from ..core import modules, vectors, parts

__author__ = 'Martin Larralde <martin.larralde@ens-paris-saclay.fr>'
__version__ = (
    __import__('pkg_resources')
        .resource_string(__name__, 'ytk.version')
        .strip()
        .decode('ascii')
)


### VECTORS

class YTKEntryVector(vectors.EntryVector):
    """A MoClo Yeast ToolKit entry vector.

    Any plasmid with two *BsmBI* restriction sites can be used to create a YTK
    entry, although the toolkit-provided entry vector (*pYTK001*) is probably
    the most appropriate plasmid to use.

    Caution:
        To the contrary of the usual MoClo entry vectors described in the
        *Weber et al.* paper, the YTK entry vectors do not provide another
        *BsaI* restriction site enclosing the placeholder sequence. As such,
        YTK Level -1 modules must embed the *BsaI* binding site.

    """

    cutter = BsmBI


class YTKCassetteVector(vectors.CassetteVector):
    """A MoClo Yeast ToolKit cassette vector.

    The YTK provides a canonical integration plasmid, preassembled from
    several other parts, that can be used as a cassette vector for an assembly
    of Type 2, 3 and 4 parts. Type 8, 8a and 678 parts are also considered as
    cassette vectors.

    References:
        *Lee et al.*, Figure 2.

    """

    cutter = BsaI


class YTKDeviceVector(vectors.DeviceVector):
    """A MoClo Yeast ToolKit multigene vector.

    Parts of Type 1 and 5 are used to order the cassette plasmids within
    the multigene assembly. The vector always contains a `ConLS` and `ConRE`
    parts.

    References:
        *Lee et al.*, Supplementary Figure S21.

    """

    cutter = BsmBI


### MODULES

class YTKProduct(modules.Product):  # FIXME ?
    """A MoClo Yeast ToolKit product.

    As the YTK entry vector does not contain the required *BsaI* restriction
    site, the site must be contained in the product sequence.

    Caution:
        The standard construction describe in the *Lee et al.* paper directly
        inserts the beginning of the *BsaI* recognition site inside of the two
        *BsmBI* overhangs at both ends of the product. Other valid constructs
        that do not proceed like so won't be considered a valid product,
        although they contain the required *BsaI* site.

    References:
        *Lee et al.*, Supplementary Figure S19.

    """

    cutter = BsmBI

    @staticmethod
    def structure():  # noqa: D105
        return (
            'CGTCTC'  # BsmBI
            'N'
            '(NNGG)'  # Product overhangs (start)
            '(TCTC'   # BsaI (first 2 nucleotides in the overhang)
            'N'
            'NNNN'    # Type specific overhang (start)
            'N*?'     # Template
            'NNNN'    # Type specific overhang (end)
            'N'
            'GA)'     # BsaI (last 4 nucleotides in the overhang)
            '(GACC)'  # Entry overhangs (end)
            'N'
            'GAGACG'  # BsmBI
        )


class YTKEntry(modules.Entry):
    """A MoClo Yeast ToolKit entry.

    YTK entries are stored and shared as plasmids flanked by *BsaI* binding
    sites at both ends of the target sequence.

    Danger:
        Although the *BsaI* binding sites is not located within the target
        sequence for almost all the standard toolkit parts, special Type 234r
        parts have these sites reversed, because these parts are used to
        assemble cassette vectors and require the final construct to contain
        a *BsaI* site to allow assembly with other parts. **Those parts will
        not match** the default `~moclo.kits.ytk.YTKEntry`, and must be
        used as `~moclo.kits.ytk.YTKPart234r` instances for the assembly
        logic to work as expected.

    """

    cutter = BsaI


class YTKCassette(modules.Cassette):
    """A MoClo Yeast ToolKit cassette.
    """

    cutter = BsmBI


### PARTS

@six.add_metaclass(abc.ABCMeta)
class YTKPart(parts.AbstractPart):
    """A Yeast MoClo ToolKit standard part.

    A part is a plasmid with standardized flanking overhangs sequences
    that allows immediate type recognition.

    Note:
        The base MoClo framework defines entries and vectors for each level
        of assembly, but the Yeast ToolKit authors instead identify their
        sequences as "parts" of different, strongly-defined types. In order
        to follow the base MoClo frameworks, parts of Type 8, 8a, and 678 are
        handled as cassette vectors, while other parts are handled as entries.
    """

    cutter = BsaI
    signature = NotImplemented


class YTKPart1(YTKPart, YTKEntry):
    """A YTK Type 1 part (**Upstream assembly connector**).

    .. image:: type1.svg
       :align: center

    Parts of this type contain non-coding and non-regulatory sequences that
    are used to direct assembly of multigene plasmids, such as ligation sites
    for other Type IIS endonucleases (e.g. *BsmBI*).

    Note:
        Official toolkit Type 1 parts also include a *EcoRI* and *XbaI* site
        just after the upstream overhang for BioBrick compatibility of the
        assembled cassettes and multi-gene plasmids.
    """

    signature = ('CCCT', 'AACG')


class YTKPart2(YTKPart, YTKEntry):
    """A YTK Type 2 part (**Promoter**).

    .. image:: type2.svg
       :align: center

    Parts of this type contain a promoter. The downstream overhang doubles as
    the start codon for the subsequent Type 3 or Type 3a coding sequence.

    Note:
        Official toolkit Type 2 parts also include a *BglII* site immediately
        preceding the start codon (overlapping the downstream overhang) for
        BglBrick compatibility.
    """

    signature = ('AACG', 'TATG')


class YTKPart3(YTKPart, YTKEntry):
    """A YTK Type 3 part (**Coding sequence**).

    .. image:: type3.svg
       :align: center

    Parts of this type contain a coding sequence, with the start codon located
    on the upstream overhang. If a stop codon is omitted from the part, and
    two bases are added before the downstream overhang, the resulting site can
    be used as a two amino acid linker to a Type 4 or 4a C-terminal fusion.

    Note:
        Official toolkit Type 3 parts also include a *BamHI* recognition site
        at the end of the included CDS (overlapping the downstream overhang)
        for BglBrick compatibility.
    """

    signature = ('TATG', 'ATCC')


class YTKPart3a(YTKPart, YTKEntry):
    """A YTK Type 3a part (**N-terminal coding sequence**).

    .. image:: type3a.svg
       :align: center

    """

    signature = ('TATG', 'TTCT')


class YTKPart3b(YTKPart, YTKEntry):
    """A YTK Type 3b part (**C-terminal coding sequence**).

    .. image:: type3b.svg
       :align: center

    Note:
        As with Type 3 parts, official toolkits Type 3b parts also include a
        *BamHI* recognition site at the end of the included CDS (overlapping
        the downstream overhang) for BglBrick compatibility.
    """

    signature = ('TTCT', 'ATCC')


class YTKPart4(YTKPart, YTKEntry):
    """A YTK Type 4 part (**Transcriptional terminator**).

    .. image:: type4.svg
       :align: center

    As Type 3 parts do not include a stop codon, parts of this type should
    encode an in-frame stop codon before the transcriptional terminator.
    Commonly used C-terminal fusions, such as purification or epitope tags,
    but it is recommended to use `~moclo.kits.ytk.YTKPart4a` and
    `~moclo.kits.ytk.YTKPart4b` subtypes instead.

    Note:
        Official toolkit Type 4 parts all start by a stop codon directly after
        the upstream overhang, followed by a *XhoI* recognition site which
        enables BglBrick compatibility, then followed by the terminator
        sequence itself.
    """

    signature = ('ATCC', 'GCTG')


class YTKPart4a(YTKPart, YTKEntry):
    """A YTK Type 4a part (**C-terminal tag sequence**).

    .. image:: type4a.svg
       :align: center

    Type 4a parts contain additional coding sequences that will be fused to
    the C-terminal extremity of the protein. These parts include, but are not
    limited to: localisation tags, purification tags, fluorescent proteins.

    Caution:
        In contrast to the Type 3 and 3b parts, the convention for 4a parts
        is to include the stop codon rather than enable read-through of the
        downstream overhang, although that convention it is not enforced.

    Note:
        Official toolkit Type 4a parts contain a stop codon after the CDS,
        itself immediately followed by a *XhoI* recognition site just before
        the downstream overhang, for BglBrick compatibility.
    """

    signature = ('ATCC', 'TGGC')


class YTKPart4b(YTKPart, YTKEntry):
    """A YTK Type 4b part (**Terminator sequence**).

    .. image:: type4b.svg
       :align: center

    Type 4b contain transcriptional terminators, but are not required to
    encode an in-frame start codon, as it should be located in the Type 4a
    part that precedes it.
    """

    signature = ('TGGC', 'GCTG')


class YTKPart234(YTKPart, YTKEntry):
    """A YTK Type 234 part (**Composite 2, 3, 4**).

    .. image:: type234.svg
       :align: center

    Type 234 parts are composed of a complete expression cassette (promoter,
    coding sequence, and terminator) fused into a single part, instead of
    separate Type 2, 3 and 4 parts.
    """

    signature = ('AACG', 'GCTG')


class YTKPart234r(YTKPart, YTKEntry):
    """A YTK Type 234 part (**Composite 2, 3, 4**) with reversed BsaI sites.

    .. image:: type234r.svg
       :align: center

    Type 234r parts are designed so that the BsaI sites are kept within the
    final cassette. They are used to assemble canonical integration vectors,
    where the Type 234 part acts as a placeholder until replaced by actual
    Type 2, 3 and 4 parts in the final construct.
    """

    signature = ('AACG', 'GCTG')

    @staticmethod
    def structure():  # noqa: D105
        return '(AACG)(NGAGACCN*?GGTCTCN)(GCTG)'


class YTKPart5(YTKPart, YTKEntry):
    """A YTK Type 5 part (**Downstream assembly connector**).

    .. image:: type5.svg
       :align: center

    As with Type 1 parts, parts of this type provide sequences such as
    restriction enzymes recognition sites, for instance in order to direct
    multigene expression plasmids.

    Note:
        Official toolkit parts also include a *SpeI* and *PstI* site at the
        end of the part sequence for BioBrick compatibility of the assembled
        cassettes and multi-gene plasmids.

    """

    signature = ('GCTG', 'TACA')


class YTKPart6(YTKPart, YTKEntry):
    """A YTK Type 6 part (**Yeast marker**).

    .. image:: type6.svg
       :align: center

    Parts of this type contain a selectable marker for *S. cerevisiae*, as
    a full expression cassette (promoter, ORF, and terminal) for conferring
    the selectable phenotype (such as drug-resistance or bioluminescence).
    """

    signature = ('TACA', 'GAGT')


class YTKPart7(YTKPart, YTKEntry):
    """A YTK Part Type 7 part (**Yeast origin** / **3' homology**).

    .. image:: type7.svg
       :align: center

    Depending on the expression organism (*E.coli* or *S. ceverisiae*), this
    sequence will either hold a yeast origin of replication, or a 3' homology
    sequence for integration in the bacterial genome.
    """

    signature = ('GAGT', 'CCGA')


class YTKPart8(YTKPart, YTKCassetteVector):
    """A YTK Type 8 part (**Bacterial origin & marker**).

    .. image:: type8.svg
       :align: center

    Parts of this type contain a bacterial origin of replication, as well
    as an antibiotic resistance marker. They act as the Golden Gate Assembly
    vector when assembling a cassette, and as such should also embbed a
    dropout sequence, such as a fluorescent protein expression cassette.

    Note:
        Official toolkit parts use an mRFP coding sequence as the dropout,
        and also include *NotI* restriction site at each end of the part to
        allow the verification of new assemblies.
    """

    signature = ('CCGA', 'CCCT')


class YTKPart8a(YTKPart, YTKCassetteVector):
    """A YTK Part 8a part (**Bacterial origin & marker**).

    .. image:: type8a.svg
       :align: center

    Parts of this type, like Type 8 parts, include a bacterial origin of
    replication and an antibiotic resistance marker, and act as Assembly
    vectors.

    Note:
        Official toolkit parts use an mRFP coding sequence as the dropout,
        and also include *NotI* restriction site at each end of the part so
        the integration plasmid can be linearized prior to transformation
        into yeast.
    """

    signature = ('CCGA', 'CAAT')


class YTKPart8b(YTKPart, YTKEntry):
    """A YTK Type 8b part (**5' homology**).

    .. image:: type8b.svg
       :align: center

    As with certain Type 7 parts, parts of this type contain long sequences
    of homology to the genome that is upstream of the target locus.
    """

    signature = ('CAAT', 'CCCT')


class YTKPart678(YTKPart, YTKCassetteVector):
    """A YTK Type 678 part (**Composite 6, 7, 8**).

    .. image:: type678.svg
       :align: center

    Type 678 parts are used when there is no requirement for yeast markers
    and origins to be included in the final assembly, for instance when
    assembling an intermediary plasmid acting as a vector for a multi-gene
    construct.
    """

    signature = ('TACA', 'CCCT')
