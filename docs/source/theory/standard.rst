Standard Modular Cloning System
===============================

.. toctree::
  :maxdepth: 1



System Definition
-----------------


.. admonition:: Definition
    :class: math-definition

    Given a genetic alphabet :math:`\langle \Sigma, \sim \rangle`, a Modular
    Cloning System :math:`S` is defined as a mathematical sequence

    .. math:: (M_l,\ V_l,\ e_l)_ {\ l\ \ge -1}

    where:

    * :math:`M_l \subseteq \Sigma^\star \cup \Sigma^{(c)}` is the set of modules
      of level :math:`l`
    * :math:`V_l \subseteq \Sigma^{(c)}` is the set of vectors of level :math:`l`
    * :math:`e_l \subseteq E` is the finite, non-empty set of *asymmetric*,
      *Type IIS* restriction enzymes of level :math:`l`


.. admonition:: Definition: :math:`k`-*cyclicity*
    :class: math-definition

    A Modular Cloning System :math:`(M_l, V_l, e_l)_ {l \ge -1}` is said to
    be :math:`k`-*cyclic after a level* :math:`\lambda` if:

    .. math::

        \begin{array}{ll}
        \exists k \in N^\star, & \\
        \forall l \ge \lambda, & \\
        & \begin{cases}
        M_{l+k} \subseteq M_l \\
        V_{l+k} \subseteq V_l \\
        e_{l+k} \subseteq e_l
        \end{cases}
        \end{array}


.. admonition:: Definition: :math:`\lambda`-limit
    :class: math-definition

    A Modular Cloning System :math:`(M_l, V_l, e_l)_ {l \ge -1}` is said to
    be :math:`\lambda`-*limited* if:

    .. math::

        \forall l \ge \lambda,
        M_l = \emptyset,
        V_l = \emptyset,
        e_l = \emptyset



Modules
-------


.. admonition:: Definition
    :class: math-definition

    For a given level :math:`l`, :math:`M_l` is defined as the set of modules :math:`m \in \Sigma^\star \cup \Sigma^{(c)}`
    for which:

    .. math::

        \begin{array}{l}
        \exists ! (S, n, k) \in e_l, \\
        \exists ! (S^\prime, n^\prime, k^\prime) \in e_l,  \\
        \exists ! (s, s^\prime) \in S \times S^\prime, \\
        \exists ! (x, y, o_5, o_3) \in (\Sigma^\star)^4, \\
        \\
        \quad \exists ! t \in \Sigma^\star,
        \left\{ \begin{array}{lll}
        \exists ! b \in \Sigma^\star,\ & m = (s \cdot x \cdot o_5 \cdot t \cdot o_3 \cdot y \cdot \widetilde{s^\prime} \cdot b)^{(c)}, & \text{ if } m \in \Sigma^{(c)}\\
        \exists ! u, v \in (\Sigma^\star)^2, & m = u \cdot s \cdot x \cdot o_5 \cdot t \cdot o_3 \cdot y \cdot \widetilde{s^\prime} \cdot v, & \text{ if } m \not \in \Sigma^{(c)}
        \end{array} \right.
        \end{array}

    with:

    * :math:`|x| = n`
    * :math:`|y| = n^\prime`
    * :math:`|o_5| = abs(k)`
    * :math:`|o_3| = abs(k^\prime)`


.. note::

    This decomposition is called the *canonic module decomposition*, where:

    * :math:`t` is the *target sequence* of the module :math:`m`
    * :math:`b` is the *backbone* of the module :math:`m` (if :math:`m` is circular)
    * :math:`u` and :math:`v` are called the *prefix* and *suffix* of the module :math:`m` (if :math:`m` is not circular)
    * :math:`o_5` and :math:`o_3` are the *upstream* and *downstream overhangs* respectively.


.. admonition:: Property
    :class: math-property

    :math:`\forall \langle \Sigma, \sim \rangle`, :math:`\forall l \ge -1`,
    :math:`\forall e_l \subset E`:

    .. math::

        M_l \text{ is a rational language }


.. admonition:: Demonstration
    :class: math-demo

    Let there be a genetic alphabet :math:`\langle \Sigma, \sim \rangle`
    and a Modular Cloning System :math:`(M_l, V_l, e_l)_ {l \ge -1}` over
    it.


    :math:`\forall l \ge -1`, the regular expression:

    .. math::

        \begin{array}l
        \bigcup_{\begin{array}l(S, n, k) \in e_l \\ (S\prime, n\prime, k\prime) \in e_l\end{array}}
        \Sigma^\star \cdot S \cdot \Sigma^n \cdot \Sigma^{abs(k)} \cdot \Sigma^\star \cdot \overline{(S | S^\prime)} \cdot \Sigma^\star \cdot \Sigma^{abs(k\prime)} \cdot \Sigma^{n\prime} \cdot \widetilde{\,S\prime\,} \cdot \Sigma^\star \\
        \end{array}

    where:

    * :math:`\star` is the `Kleene star <https://en.wikipedia.org/wiki/Kleene_star>`_.
    * :math:`\widetilde{S} = \{\widetilde{s}, s \in S\}` (:ref:`reverse complementation <reverse-complement>` operator).
    * :math:`\overline{S} = \{w \in \Sigma^\star, w \not \in S\}` (`complement <https://en.wikipedia.org/wiki/Complement_(set_theory)>`_ operator).
    * :math:`S | S^\prime = S \cup S^\prime` (`alternation <https://en.wikipedia.org/wiki/Alternation_(formal_language_theory)>`_ operator).


    matches a sequence :math:`m \in \Sigma^\star \cup \Sigma^{(c)}` if and only if
    :math:`m \in M_l`.

    :math:`M_l` is regular, so given Kleene's Theorem, :math:`M_l` is rational.



Vectors
-------

.. admonition:: Definition
    :class: math-definition

    For a given level :math:`l`, :math:`V_l` is defined as the set of vectors :math:`v \in \Sigma^{(c)}`
    for which:

    .. math::

        \begin{array}{l}
        \exists ! (S, n, k) \in e_l, \\
        \exists ! (S^\prime, n^\prime, k^\prime) \in e_l,  \\
        \exists ! (s, s^\prime) \in S \times S^\prime, \\
        \exists ! (x, y, o_5, o_3) \in (\Sigma^\star)^4, \\
        \\
        \quad \exists ! (b, p) \in (\Sigma^\star)^2,
        \exists ! b \in \Sigma^\star,\ v = (o_3 \cdot b \cdot o_5 \cdot y \cdot \widetilde{s} \cdot p \cdot s\prime \cdot x)^{(c)} \\
        \end{array}

    with:

    * :math:`|x| = n`
    * :math:`|y| = n^\prime`
    * :math:`|o_5| = abs(k)`
    * :math:`|o_3| = abs(k^\prime)`
    * :math:`o_3 \ne o_5`


.. note::

    This decomposition is called the *canonic vector decomposition*, where:

    * :math:`p` is the *placeholder sequence* of the vector :math:`v`
    * :math:`b` is the *backbone* of the vector :math:`v`
    * :math:`o_3` and :math:`o_5` are the *upstream* and *downstream overhangs* respectively.



Overhangs
---------

By definition, every valid level :math:`l` module and vector only have a single canonic
decomposition where they have unique :math:`o_5` and :math:`o_3` overhangs. As such,
let the function :math:`up` (resp. :math:`down`) be defined as the function which:

* to a module :math:`m` associates the word :math:`o_5` (resp. :math:`o_3`) from its
  canonic module decomposition
* to a vector :math:`v` associates the word :math:`o_3` (resp. :math:`o_5`) from its
  canonic vector decomposition.



Standard Assembly
-----------------


.. admonition:: Definition: *Standard MoClo Assembly*
    :class: math-definition

    Given an assembly of level :math:`l`, where :math:`m_1, \dots, m_k \in M_l^k, v \in V_l`:

    .. math::

        a:\quad m_1 + \dots + m_k \xrightarrow{\quad e_l \quad} A \subset (\Sigma^\star \cup \Sigma^{(c)})

    and the partial order :math:`le` over :math:`S = \{m_1, \dots, m_k\}` defined as:

    .. math::

        \begin{array}{l}
        \forall x, y \in S^2, \\
        \quad x \le y \iff \begin{cases}
        x = y & \\
        down(x) = up(y) & \text{ if } x \ne y\\
        \exists z \in S \backslash \{x, y\}, down(x) = up(z), \ z \le y & \text{ if } x \ne y \text{ and } down(x) \ne up(y)
        \end{cases}
        \end{array}

    then a chain :math:`\langle S\prime, \le \rangle \subset \langle S, \le \rangle` is
    an *insert* if:

    .. math::

        \begin{cases}
        v \le min(S^\prime) \\
        max(S^\prime) \le v
        \end{cases}
        \iff
        \begin{cases}
        down(v) = up(min(S^\prime)) \\
        up(v) = down(max(S^\prime))
        \end{cases}

    :math:`a` is:

    * *invalid* if  :math:`\langle S, \le \rangle` is an antichain or :math:`\langle S, \ge \rangle`
      has no insert.
    * *valid* if :math:`\langle S, \le \rangle` has at least one insert.
    * *ambiguous* if :math:`\langle S, \le \rangle` has more than one insert.
    * *unambiguous* if :math:`\langle S, \le \rangle` has exactly one insert.
    * *complete* if :math:`\langle S, \le \rangle` is an insert.


.. admonition:: Corollary
    :class: math-property

    If an assembly :math:`a` is complete, then there exist a permutation
    :math:`\pi` of :math:`[\![1, k]\!]` such that:

    .. math:: m_{\pi(1)} \le m_{\pi(2)} \le \dots \le m_{\pi(k-1)} \le m_{\pi(k)}

    and:

    .. math::

        \begin{array}{lll}
        up(m_{\pi(1)}) &=& down(v) \\
        down(m_{\pi(k)}) &=& up(v)
        \end{array}


.. admonition:: Property: *Uniqueness of the cohesive ends*
   :class: math-property

   If an assembly

   .. math:: m_1 + \dots + m_k \xrightarrow{\quad e_l \quad} A \subset (\Sigma^\star \cup \Sigma^{(c)})

   is unambiguous and complete, then :math:`\forall i \in [\![1, k]\!]`,

   .. math::

       \left\{
       \begin{array}{llll}
       up(m_i) &\ne& down(m_i)& \\
       up(m_i) &\ne& up(m_j),     & j \in [\![1, k]\!] \backslash \{i\} \\
       down(m_i) &\ne& down(m_j), & j \in [\![1, k]\!] \backslash \{i\} \\
       \end{array}
       \right .



.. admonition:: Demonstration
    :class: math-demo

    Let there be an unambiguous complete assembly

    .. math:: a:\quad m_1 + \dots + m_k \xrightarrow{\quad e_l \quad} A



    .. rubric:: :math:`up(m_i) \ne down(m_i)`

    Let's suppose that :math:`\exists i \in [\![1, k]\!]` such that

    .. math:: up(m_i) = down(m_i)

    then :math:`\langle \{m_1, \dots, m_k\} \backslash \{m_i\}, \le \rangle`
    is also an insert, which cannot be since :math:`a` is complete.



    .. rubric:: :math:`up(m_i) \ne up(m_j)`


    Let's suppose that :math:`\exists (i, j) \in [\![1, k]\!]^2` such that

    .. math:: up(m_i) = up(m_j)

    Since the :math:`a` is complete, there exists :math:`pi` such that

    .. math:: m_{\pi(1)} \le m_{\pi(2)} \le \dots \le m_{\pi(k-1)} \le m_{\pi(k)}

    and since :math:`a` is unambiguous, :math:`\langle \{m_1, \dots, m_k\}, \le \rangle`
    is the only insert.

    Or if :math:`up(m_i) = up(m_j)`, then

    .. TODO

    .. .. math:: \langle \{m_1, \dots, m_k\} \backslash \{m_\pi(i), \le \rangle

    is also an insert, which cannot be since :math:`a` is unambiguous.



    .. rubric:: :math:`down(m_i) \ne down(m_j)`

    TODO


.. admonition:: Property: *Uniqueness of the assembled plasmid*
    :class: math-property

    If an assembly

    .. math:: m_1 + \dots + m_k \xrightarrow{\quad e_l \quad} A \subset (\Sigma^\star \cup \Sigma^{(c)})

    is unambiguous, then

    .. math:: A \cap \Sigma^{(c)} = \{p\}

    with

    .. math:: p = \left( up(v) \cdot b \cdot up(m_{\pi(1)}) \cdot t_{\pi(1)} \cdot \, \dots \, \cdot up(m_{\pi(n)}) \cdot t_{\pi(n)} \right) ^{(c)}

    (:math:`n \le k`, :math:`n = k` if :math:`a` is complete).




.. admonition:: Demonstration
    :class: math-demo

    TODO
