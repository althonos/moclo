Concepts
========

Introduction
------------

The MoClo standard was first presented in the *Weber et al., 2011* [21364738]_
paper, as an attempt to standardize the process of assembling complex DNA
molecules from smaller genetic elements. It is inspired by two previous
standards:

- NOMAD [8855278]_, which proposed generic notions of *modules* and *vectors*,
  as well as assembly using Type IIS enzymes. Modules can be combined in any
  order, but are clone sequentially one module at a time.
- BioBrick [18410688]_, which defines *parts* with a stable structure: assembling
  two parts together always gives a part with the same flanking restriction
  sites.

The MoClo standard enhances both of these assembly standards by relying on the
Golden Gate Assembly, which allows single-step assembly of an arbitrary number
of modules into a vector. Furthermore, MoClo parts are flanked by stereotypical
overhangs, enforcing a particular assembly order, therefore allowing only the
desired contruct to be obtained.


Type II-S enzymes
-----------------

Restriction enzymes are enzymes that are able to cut DNA at or near specific
recognition sites. Among those enzymes, Type IIS enzymes cut DNA out of the
sequence they recognize, at a defined distance. The cut can produce *cohesive ends*,
which can then recombine with other sequences sharing the complementary cohesive
ends, or *blunt ends*, which cannot recombine. The design of the cohesive ends
is of great importance when using Type II-S enzymes to do molecular cloning.


Golden Gate Assembly
--------------------

The Golden Gate Assembly relies on Type II-S enzymes to assemble several DNA
sequences. The sequences are first cut by restriction enzymes, and then
assembled together using a T4 DNA ligase. These two steps can be repeated in
a single reaction tube using a *thermo cycler*, as the two enzymes typically do
not work at the same temperature. As standard Type II-S enzymes, such as *BsaI*
or *BsmBI*, create a 4-base-long cohesive end when cutting the DNA, there can be
as much as 256 fragments combined together in a deterministic way in a single
assembly, although *in vivo* the chemical properties of the nucleotides will
most likely prevent assemblies that large to succeed.

.. figure:: assembly.svg
    :align: center

    *Example GoldenGate assembly of two modules in a vector using*
    `BsaI <https://international.neb.com/products/r0535-bsai#Product%20Information>`_.



The MoClo system
----------------

The MoClo system combines the idea of a standard *part* format from the BioBrick
standard, with the Golden Gate assembly protocol, allowing several modules to be
assembled in a vector at the same time.

Hierarchy
'''''''''

MoClo modules and vectors are divided into several levels, describing their
structural and transcriptional features:

- Level -1 modules are sequences that are not yet in a standardized backbone,
  but can be assembled in a dedicated vector to form a level 0 module. They are
  most of the time obtained via oligonucelotide synthesis, or PCR.
- Level 0 modules are standardized genetic elements: promoter, 5' UTR,
  signal sequence, CDS, terminator.
- Level 1 modules are transcription units, formed by a combination of Level 0
  modules, and are able to express proteins
- Level 2 modules are multigenic units, containing several transcription units,
  and are able to express many genes at onces.

Furthermore, the enzyme used during the Golden Gate Assembly depends on the
assembly level. Alternating between the two enzymes makes it possible for an
infinite number of genes to be inserted in the same plasmid, although biological
limits are reached *in vivo*.


Types definition
''''''''''''''''

Although transcription units can be assembled in any possible order in their
destination vectors, level 0 modules must be assembled in a specific order to
obtain a functional genetic construct. In order to enforce the assembly order,
parts are flanked by fusion sites with standard sequences, which are unique to
the *type* of the part. A valid level 1 module is obtained by assembling a part
of each type into the destination vector.


Assembly markers
''''''''''''''''

Once the Golden Gate Assembly is finished, the obtained constructs can be
amplified using a bacterial host. After transformation, bacteria are selected
using two different factors:

- An antibiotic for which a resistance cassette is only availble on the vector,
  but not on any module: this allows selecting all the bacterias that received
  the vector plasmid
- A marker for a dropout reporter gene that can only be found in the vector but
  not in the final construct (such as the *gfp* or *lacZ* genes).

This double screening makes it possible to select only the bacterias that
contain the expected construct, discarding the others, and retrieving the
assembled plasmid through a miniprep.


References
----------


.. [8855278] Rebatchouk, D, N Daraselia, and J O Narita.
             ‘NOMAD: A Versatile Strategy for in Vitro DNA Manipulation Applied to Promoter Analysis and Vector Design.’
             *Proceedings of the National Academy of Sciences of the United States of America* 93, no. 20 (1 October 1996): 10891–96.
             `pmid:8855278 <https://www.ncbi.nlm.nih.gov/pubmed/8855278>`_

.. [18410688] Shetty, Reshma P, Drew Endy, and Thomas F Knight.
              ‘Engineering BioBrick Vectors from BioBrick Parts’.
              *Journal of Biological Engineering* 2 (14 April 2008): 5.
              `doi:10.1186/1754-1611-2-5 <https://doi.org/10.1186/1754-1611-2-5>`_

.. [21364738] Weber, Ernst, Carola Engler, Ramona Gruetzner, Stefan Werner, and Sylvestre Marillonnet.
              ‘A Modular Cloning System for Standardized Assembly of Multigene Constructs’.
              *PLOS ONE* 6, no. 2 (18 February 2011): e16765.
              `doi:10.1371/journal.pone.0016765 <https://doi.org/10.1371/journal.pone.0016765>`_
