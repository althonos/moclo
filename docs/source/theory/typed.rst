Typed Modular Cloning System
============================

.. toctree::
  :maxdepth: 1



System Definition
-----------------


.. admonition:: Definition
    :class: math-definition

    Given a genetic alphabet :math:`\langle \Sigma, \sim \rangle`, a
    Typed Modular Cloning System :math:`S` is defined as a mathematical
    sequence

    .. math:: (M_l,\ V_l,\ \mathcal{M}_l,\ \mathcal{V}_l,\ e_l)_ {\ l\ \ge -1}

    where:

    * :math:`(M_l, V_l, e_l)_{l \ge -1}` is a standard Modular Cloning System
    * :math:`\mathcal{M}_l \subseteq \mathcal{P}(M_l) \to \mathcal{P}(M_l)`
      is the set of *module types* of level :math:`l`
    * :math:`\mathcal{V}_l \subseteq \mathcal{P}(V_l) \to \mathcal{P}(V_l)`
      is the set of *vector types* of level :math:`l`



Types
-----

.. admonition:: Definition
    :class: math-definition

    :math:`\forall l \ge -1`, we define types using their signatures (*i.e.* the
    sets of upstream and downstream overhangs of elements using this type):

    .. math::

        \begin{array}{ll}
        \forall t \in \mathcal{M}_l,& \begin{cases}
        Up(t) &= \bigcup_{m \in t(M_l)} \{ up(m) \} \\
        Down(t) &= \bigcup_{m \in t(M_l)} \{ down(m) \}
        \end{cases} \\
        \forall t \in \mathcal{V}_l,& \begin{cases}
        Up(t) &= \bigcup_{v \in t(V_l)} \{ up(v) \} \\
        Down(t) &= \bigcup_{v \in t(V_l)} \{ down(v) \}
        \end{cases}
        \end{array}


.. admonition:: Corollary
    :class: math-property

    :math:`\forall l \ge -1`,

    .. math::

        \begin{array}{lll}
        \forall t \in \mathcal{M}_l,&\ t(M_l) &= \{ m \in M_l\ |\ up(m) \in Up(t),\ down(m) \in Down(t) \} \\
        \forall t \in \mathcal{V}_l,&\ t(V_l) &= \{ v \in V_l\ |\ up(v) \in Up(t),\ down(v) \in Down(t) \}
        \end{array}


.. admonition:: Property: *Structural equivalence of module types*
    :class: math-property

    Given a valid (*resp.* unambiguous) (*resp.* complete) assembly

    .. math::

        m_1 + \dots + m_k + v \xrightarrow{e_l} A \subset (\Sigma^\star \cup \Sigma^{(c)})

    then if there exist :math:`t \in \mathcal{M}_l` such that

    .. math::

        \begin{cases}
        \lvert Up(t) \rvert = \lvert Down(t) \rvert = 1 \\
        m_1 \in t(M_l)
        \end{cases}

    then :math:`\forall m_1\prime \in t(M_l)`,

    .. math::

        m_1\prime + \dots + m_k + v \xrightarrow{e_l} A \subset (\Sigma^\star \cup \Sigma^{(c)})

    is valid (*resp.* unambiguous) (*resp.* complete).
