Preliminary Definitions
=======================

.. toctree::
  :maxdepth: 1



Genetic Alphabet
----------------

.. admonition:: Definition
    :class: math-definition

    A *genetic alphabet* :math:`\langle \Sigma,\sim \rangle` is an algebraic structure on an alphabet :math:`\Sigma`
    with a unary operation :math:`\sim` verifying the following properties:

    * :math:`\sim: \Sigma^\star \to \Sigma^\star` is a bijection
    * :math:`\forall x \in \Sigma^\star, \lvert \widetilde{x} \rvert = \lvert x \rvert`
    * :math:`\forall (x, y) \in (\Sigma^\star)^2, \quad \widetilde{x \cdot y} = \widetilde{\,y\,} \cdot \widetilde{\,x\,}`


.. note::

    To stay consistent with the biology lexicon, we will be referring to a word
    over a genetic alphabet as a *sequence*, only explicitly naming a mathematical
    sequence when needed to.


.. admonition:: Examples
    :class: math-example

    * :math:`(\{A, T, G, C\}, \sim)` is the standard genetic alphabet, with
      :math:`\sim` defined as:

      .. math::

         \forall x \in \{A, T, G, C\}, \widetilde{x} = \overline{\mu(x)}

      where:

      * :math:`\overline{y}` is the *reversal* of :math:`y`
      * :math:`\mu` is the *automorphism* defined as :math:`\mu(A) = T`, :math:`\mu(C) = G`

    * :math:`(\{A, T, G, C, d5SICS, dNaM\}, \sim)` is the genetic alphabet using the
      unnatural base pairs from
      `Malyshev et al., Nature 2014 <https://www.nature.com/articles/nature13314>`_
      defined as:

      .. math::

          \forall x \in \{A, T, G, C, d5ICS, dNaM\}, \widetilde{x} = \overline{\mu\prime(x)}

      where:

      * :math:`\mu\prime` is the automorphism defined as :math:`\mu\prime(A) = T`,
        :math:`\mu\prime(C) = G`, :math:`\mu\prime(d5SICS) = dNaM`



Circular Sequences
------------------


.. admonition:: Definition
   :class: math-definition

   A *circular word* over an alphabet :math:`\Sigma` is a finite word with no end. It
   can be noted :math:`w^{(c)}`, where :math:`w` is a finite word of :math:`\Sigma^\star`.


.. admonition:: Definition: *Cardinality*
   :class: math-definition

   Given a circular sequence :math:`s^{(c)}`, the cardinal of :math:`s^{(c)}`,
   noted :math:`\lvert s^{(c)} \rvert`, is defined as:

   .. math:: \lvert s^{(c)} \rvert = \lvert s \rvert


.. admonition:: Definition: *Equality*
    :class: math-definition

    Given two sequences :math:`a^{(c)}` and :math:`b^{(c)}` with

    .. math::

        \begin{array}{lllll}
        a &=& a_0 \cdot a_1 \cdot \, \dots \, \cdot a_m & \in \Sigma^{(m)}, & m \in \mathbb{N} \\
        b &=& b_0 \cdot b_1 \cdot \, \dots \, \cdot b_n & \in \Sigma^{(n)}, & n \in \mathbb{N}
        \end{array}

    let the :math:`=` relation be defined as:

    .. math:: a = b \iff \exists k \in \mathbb{N}, a = \sigma^{k}(b)

    where :math:`\sigma` is the circular shift defined as:

    .. math::

        \begin{array}l
        \forall u = u_1 \cdot u_2 \cdot\,\dots\,\cdot u_k \in \Sigma^k, \\
        \quad \quad \sigma(u_1 \cdot u_2 \cdot\,\dots\,\cdot u_k) =
        u_k \cdot u_1 \cdot u_2 \cdot \, \dots \, \cdot u_{k-1}
        \end{array}



.. admonition:: Property
    :class: math-property

    :math:`=` is a relation of equivalence over :math:`\Sigma^{(c)}`


.. admonition:: Demonstration
    :class: math-demo

    Given the set of circular sequences :math:`\Sigma^{(c)}` using an alphabet
    :math:`\Sigma`:

    * **Reflexivity**:

      .. math::

          s^{(c)} \in \Sigma^{(c)} \implies s = Id(s) = \sigma^{0}(s) \implies s^{(c)} = s^{(c)}

    * **Symetry**:
      :math:`\forall s_1^{(c)}, s_2^{(c)} \in \Sigma^{(c)} \times \Sigma^{(c)}`:

      .. math::

          \begin{array}{lll}
          s_1^{(c)} = s_2^{(c)}
          &\iff& \exists k \in \mathbb{N}, s_1 = \sigma^k(s_2) \\
          &\iff& \exists k \in \mathbb{N}, s_2 = \sigma^{-k}(s_1) \\
          &\iff& \exists k \in \mathbb{N}, s_2 = \sigma^{\lvert s_1 \rvert - k}(s_1) \\
          &\iff& s_2^{(c)} = s_1^{(c)}
          \end{array}

    * **Transitivity**:
      :math:`\forall s_1, s_2, s_3 \in \Sigma^{(c)} \times \Sigma^{(c)} \times \Sigma^{(c)}`

      .. math::

          \begin{array}{lll}
          \begin{cases}
          s_1^{(c)} = s_2^{(c)} \\
          s_2^{(c)} = s_3^{(c)}
          \end{cases}
          &\implies& \begin{cases}
          \exists k_1 \in \mathbb{N}, s_1 = \sigma^{k_1}(s_2) \\
          \exists k_2 \in \mathbb{N}, s_2 = \sigma^{k_2}(s_3)
          \end{cases} \\
          &\implies&
          \exists k_1, k_2 \in \mathbb{N}^2, s_1 = \sigma^{k_1} \circ \sigma^{k_2}(s_3) \\
          &\implies&
          \exists k_1, k_2 \in \mathbb{N}^2, s_1 = \sigma^{k_1 + k_2}(s_3) \\
          &\implies& s_1^{(c)} = s_3^{(c)}
          \end{array}



.. admonition:: Definition: *Automaton acception*
    :class: math-definition

    Given a finite automate :math:`A` over an alphabet :math:`\Sigma`, and
    :math:`u^{(c)}` a sequence of :math:`\Sigma^{(c)}`, :math:`A` *accepts*
    :math:`u^{(c)}` iff there exist a sequence :math:`v` of :math:`\Sigma^\star`
    such that:

    * :math:`v^{(c)} = u^{(c)}`
    * :math:`A` accepts :math:`v`



Restriction Enzymes
-------------------

.. admonition:: Definition
    :class: math-definition

    Given a genetic alphabet :math:`\langle \Sigma, \sim \rangle`, a restriction enzyme :math:`e` can
    be defined as a tuple :math:`(S, n, k)` where:

    * :math:`S \subseteq \Sigma^\star` is the finite set of *recognition sites*
      that :math:`e` binds to
    * :math:`\forall (s, s\prime) \in S^2, \lvert s \rvert = \lvert s\prime \rvert`
    * :math:`n \in \mathbb{Z}` is the *cutting offset* between the last nucleotides
      of the site and the first nucleotide of the restriction cut
    * :math:`n \ge \lvert -s \rvert, s \in S`
    * :math:`k \in \mathbb{Z}` is the *overhang length*:

      * :math:`k = 0` if the enzyme produces blunt cuts
      * :math:`k > 0` if the enzyme produces :math:`5\prime` overhangs
      * :math:`k < 0` if the enzyme produce :math:`3\prime` overhangs


.. note::

    This definition only covers single-cut restriction enzymes found *in vivo*,
    but we don’t need to cover the case of double-cut restriction enzymes since
    they are not used in modular cloning.

.. admonition:: Definition: *Enzyme types*
    :class: math-definition

    A restriction enzyme :math:`(S, n, k)` is:

    * a *blunt cutter* is :math:`k = 0`
    * an *asymmetric cutter* if :math:`k \ne 0`
    * a *Type IIS* enzyme if:

      * :math:`n \ge 0`
      * :math:`\forall s \in S, s \ne \overline{s}`



Golden Gate Assembly
--------------------

.. admonition:: Definition
    :class: math-definition

    An *assembly* is a function of :math:`\mathcal{P}(\Sigma^\star \cup \Sigma^{(c)}) \times \mathcal{P}(E)`
    to :math:`\mathcal{P}(\Sigma^\star \cup \Sigma^{(c)})`, which to a set
    of distinct sequences :math:`\{d_1, \dots, d_m\}` and a set of restriction
    enzymes :math:`\{e_1, \dots, e_n\}` associates the set of digested/ligated sequences
    :math:`A = \{a_1, \dots a_k\}`.

    The notation for an assembly is:

    .. math:: d_1 + \dots d_m \xrightarrow{\quad e_1, \dots, e_n \quad} A
