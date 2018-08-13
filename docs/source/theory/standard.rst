Standard Modular Cloning system
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
    :class: math-demonstration

    Let there be a genetic alphabet :math:`\langle \Sigma, \sim \rangle`
    and a Modular Cloning System :math:`(M_l, V_l, e_l)_ {l \ge -1}` over
    it.

    :math:`\forall l \ge -1`, the regular expression:

    .. math::

        \begin{array}l
        \bigcup_{\begin{array}l(S, n, k) \in e_l \\ (S\prime, n\prime, k\prime) \in e_l\end{array}}
        \bigcup_{\begin{array}l s \in S \\ s\prime \in S\prime\end{array}}
        x^\star \cdot  s \cdot x^n \cdot x^{abs(k)} \cdot x^\star \cdot x^{abs(k\prime)} \cdot x^{n\prime} \cdot \widetilde{\,s\prime\,} \cdot x^\star \\
        \end{array}

    where:

    * :math:`x` is the regular expression accepting any single letter
      of :math:`\Sigma`

    matches a sequence :math:`m \in \Sigma^\star \cup \Sigma^{(c)}` if and only if
    :math:`m \in M_l`.

    :math:`M_l` is regular, so given Kleene's Theorem, :math:`M_l` is rational.
