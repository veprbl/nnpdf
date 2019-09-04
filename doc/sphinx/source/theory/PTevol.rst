| **Notes on Perturbative Evolution and PDF flavor decomposition**
| AG & MU (with the help of Jacob Haddo, summer student)

Perturbative PDF evolution
==========================

Notation
--------

* **The strong coupling constant**

Define the coupling

.. math:: a_{s} = \frac{\alpha_{s}(Q^{2})}{4\pi}

.. math:: a_{0} = a_{s}(Q_{0}^{2})

which satisfies the Renormalisation Group Equation

.. math:: \frac{da_{s}}{d\ln\mu^{2}} = \beta(a_{s}) = - \sum_{n = 0}^{\infty}\beta_{n}a_{s}^{n + 2}\,,

where

.. math:: \beta_0 = \frac{11}{3}C_A - \frac{4}{3}T_FN_f

.. math:: \beta_1 = \frac{34}{3}C^2_A - 4C_FT_FN_f - \frac{20}{3}C_AT_FN_f

.. math:: \beta_2 = \frac{2857}{54}C^3_A + 2C^2_FT_FN_f - \frac{205}{9}C_FC_AT_FN_f - \frac{1415}{27}C^2_AT_FN_f - \frac{44}{9}C_FT^2_FN^2_f - \frac{158}{27}C_AT^2_FN^2_f.


* **Mellin transform**

The Mellin transform of a function is defined as

.. math:: f(N,Q^{2}) = \int_{0}^{1}dx\, x^{N - 1}f(x,Q^{2})\,,

and we can get back the x-space distribution as

.. math:: f(x,Q^{2}) = \int_{c - i\infty}^{c + i\infty}\mspace{6mu}\frac{dN}{2\pi i}\, x^{- N}f(N,Q^{2})\,,

where the intercept c of integration contour is chosen to be to the
right of all singularities of f(N,Q2) in the complex N plane.

Parton evolution
--------------------

The scale dependence of the parton distribution functions is described
by the renormalisation group equations for mass factorisation (DGLAP)

.. math:: \mu^{2}\frac{\partial}{\partial\mu^{2}}f_{i}(x,\mu^{2}) = P_{ij}(x,\mu^{2}) \otimes f(x,\mu^{2})\,

where f\ :sub:`i` is the generic parton distribution function, P\ :sub:`ij` are the
Altarelli-Parisi kernels and :math:`\otimes` denotes the Mellin convolution

.. math:: f(x) \otimes g(x) \equiv \int_{x}^{1}dyf(y)g\left( \frac{x}{y} \right)

We have a system of (2n\ :sub:`f` + 1) coupled
integro-differential equations, where the summation over the parton
species j is understood.

The N\ :sup:`m`\ LO approximation for the splitting functions :math:`P_{ij}(x,\mu^2)`

.. math:: P_{ij}^{N^{m}LO}(x,\mu^{2}) = \sum_{k = 0}^{m}a_{s}^{k + 1}(\mu^{2})P_{ij}^{(k)}(x)

where we note that the only dependence on the scale :math:`\mu^2`
is through the coupling constant :math:`a_s(\mu^2)`. The splitting
functions in the case of unpolarised partons are known up to NNLO and,
in the notation we adopt, their explicit expressions are found in .

In the following, to describe the solution to the DGLAP evolution
equations we will be working in Mellin space where, as we have seen,
convolutions are turned into products.

Flavour decomposition
---------------------

The primary quantities are the :math:`2n_f` quark and antiquark
distributions qi(x,Q2), Qi(x,Q2) and the gluon distribution g(x,Q2).

From considerations based on charge conjugation and flavour symmetry it
is possible to rewrite the system of equations as :math:`(2N_f - 1)` equations
describing the
independent evolution of the non-singlet quark asymmetries and

.. math:: q_{NS,ij}^\pm = q_i \pm Q_i - (q_j \pm Q_j)

.. math:: q_{NS}^v = \sum_{i = 1}^{N_f}(q_i - Q_i)

and a system of 2 equations describing the coupled evolution of the
singlet and gluon parton distributions.

.. math::

   \begin{matrix}
   \mu^{2}\frac{\partial}{\partial\mu^{2}}q_{NS}^{\pm ,v}(x,\mu^{2}) & = & P_{NS}^{\pm ,v} \otimes q_{NS}^{\pm ,v}(x,\mu^{2}) \\
   \mu^{2}\frac{\partial}{\partial\mu^{2}}\begin{pmatrix}
   \Sigma \\
   g \\
   \end{pmatrix}(x,\mu^{2}) & = & \begin{pmatrix}
   P_{qq} & P_{qg} \\
   P_{gq} & P_{gg} \\
   \end{pmatrix} \otimes \begin{pmatrix}
   \Sigma \\
   g \\
   \end{pmatrix}(x,\mu^{2}) \\
   \end{matrix}

where the singlet combination, :math:`\Sigma`, is defined as

.. math:: \Sigma = \sum_{i = 1}^{N_{f}}(q_{i} + {\overline{q}}_{i})\,,

where :math:`N_{f}` is the number of *light flavors*, *i.e.* the number
of flavors with :math:`m_{q}^{2} < Q^{2}`.

At LO
:math:`P_{NS}^{(0), +} = P_{NS}^{(0), -} = P_{NS}^{(0),v} = P_{qq}^{(0)}`.
At NLO :math:`P_{NS}^{(0), -} = P_{NS}^{(0),v}` while all the other
splitting functions are different. Starting form :math:`\mathcal{O}(\alpha_s^2)`
all splitting functions are different from each other.

The evolution of the individual quark distributions with the scale can
be computed by introducing the following set of non-singlet
distributions:

.. math:: \begin{matrix} V & = & u^{-} + d^{-} + s^{-} + c^{-} + b^{-} + t^{-} \\ \end{matrix}

.. math:: \begin{matrix} V_{3} & = & u^{-} - d^{-} \\ \end{matrix}

.. math:: \begin{matrix} V_{8} & = & u^{-} + d^{-} - 2s^{-} \\ \end{matrix}

.. math:: \begin{matrix} V_{15} & = & u^{-} + d^{-} + s^{-} - 3c^{-} \\ \end{matrix}

.. math:: \begin{matrix} V_{24} & = & u^{-} + d^{-} + s^{-} + c^{-} - 4b^{-} \\ \end{matrix}

.. math:: \begin{matrix} V_{35} & = & u^{-} + d^{-} + s^{-} + c^{-} + b^{-} - 5t^{-} \\ \end{matrix}

.. math:: \begin{matrix} T_{3} & = & u^{+} - d^{+} \\ \end{matrix}

.. math:: \begin{matrix} T_{8} & = & u^{+} + d^{+} - 2s^{+} \\ \end{matrix}

.. math:: \begin{matrix} T_{15} & = & u^{+} + d^{+} + s^{+} - 3c^{+} \\ \end{matrix}

.. math:: \begin{matrix} T_{24} & = & u^{+} + d^{+} + s^{+} + c^{+} - 4b^{+} \\ \end{matrix}

.. math:: \begin{matrix} T_{35} & = & u^{+} + d^{+} + s^{+} + c^{+} + b^{+} - 5t^{+} \\ \end{matrix}

where :math:`q_{i}^{\pm} = q_{i} \pm {\overline{q}}_{i}`, and
:math:`u,d,s,c,b,t` are the various flavour distributions.

The combinations :math:`V_{j}` and :math:`T_{j}` evolve according to eq.
(`[eq:DGLAPdecomp] <#eq:DGLAPdecomp>`__) with :math:`P_{NS}^{-}` and
:math:`P_{NS}^{+}` respectively, while the total valence :math:`V`
evolves with the :math:`P_{NS}^{v}` kernel. Inverting the linear system
Eq.\ `[eq:lincomb] <#eq:lincomb>`__ we obtain the individual pdf’s as a
function of the evolved non-singlet and singlet distributions:

.. math::

   \begin{matrix}
   u & = & (10\Sigma + 30T_{3} + 10T_{8} + 5T_{15} + 3T_{24} + 2T_{35} + 10V + 30V_{3} + 10V_{8} + 5V_{15} + 3V_{24} + 2V_{35})/120 \\
   \overline{u} & = & (10\Sigma + 30T_{3} + 10T_{8} + 5T_{15} + 3T_{24} + 2T_{35} - 10V - 30V_{3} - 10V_{8} - 5V_{15} - 3V_{24} - 2V_{35})/120 \\
   d & = & (10\Sigma - 30T_{3} + 10T_{8} + 5T_{15} + 3T_{24} + 2T_{35} + 10V - 30V_{3} + 10V_{8} + 5V_{15} + 3V_{24} + 2V_{35})/120 \\
   \overline{d} & = & (10\Sigma - 30T_{3} + 10T_{8} + 5T_{15} + 3T_{24} + 2T_{35} - 10V + 30V_{3} - 10V_{8} - 5V_{15} - 3V_{24} - 2V_{35})/120 \\
   s & = & (10\Sigma - 20T_{8} + 5T_{15} + 3T_{24} + 2T_{35} + 10V - 20V_{8} + 5V_{15} + 3V_{24} + 2V_{35})/120 \\
   \overline{s} & = & (10\Sigma - 20T_{8} + 5T_{15} + 3T_{24} + 2T_{35} - 10V + 20V_{8} - 5V_{15} - 3V_{24} - 2V_{35})/120 \\
   c & = & (10\Sigma - 15T_{15} + 3T_{24} + 2T_{35} + 10V - 15V_{15} + 3V_{24} + 2V_{35})/120 \\
   \overline{c} & = & (10\Sigma - 15T_{15} + 3T_{24} + 2T_{35} - 10V + 15V_{15} - 3V_{24} - 2V_{35})/120 \\
   b & = & (5\Sigma - 6T_{24} + T_{35} + 5V - 6V_{24} + V_{35})/60 \\
   \overline{b} & = & (5\Sigma - 6T_{24} + T_{35} - 5V + 6V_{24} - V_{35})/60 \\
   t & = & (\Sigma - T_{35} + V - V_{35})/12 \\
   \overline{t} & = & (\Sigma - T_{35} - V + V_{35})/12 \\
   \end{matrix}

N space solutions to DGLAP evolution equations
----------------------------------------------

We pointed out before that the splitting functions (and therefore the
anomalous dimensions) depend on the scale only through the coupling
constant. We can then choose :math:`a_{s}` as evolution variable and
rewrite the DGLAP evolution equation for the quark-singlet and gluon
distributions, in Mellin-\ :math:`N` space, as.

   .. math:: a_s\frac{\partial}{\partial a_s} \binom{\Sigma}{g}(N, a_s) = -\mathbf{R} \cdot \binom{\Sigma}{g}(N, a_s),

where the matrix **R** has the following perturbative expansion

.. math:: \mathbf{R} = \mathbf{R}_0+a_s\mathbf{R}_1+a_s\mathbf{R}_2 + \dots

with

.. math:: \mathbf{R}_0 \equiv \frac{\boldsymbol{\gamma}^{(0)}}{\beta_0}

.. math:: \mathbf{R}_k \equiv \frac{\boldsymbol{\gamma}^{(k)}}{\beta_0} - \sum_{i=1}^k \frac{\beta_i}{\beta_0}R_{k-i}

where the :math:`\mathbf{\gamma}` stands for the matrix of anomalous dimensions.

The solution of the singlet evolution equation at leading order is:

.. math:: \mathbf{q}_{LO}(x,Q^2) = \mathbf{L}(a_s,a_0,N)\mathbf{q}_{LO}(x,Q_0^2).

The leading order evolution operator :math:`\mathbf{L}` is written, in terms
of the eigenvalues of the leading order anomalous dimension matrix

.. math:: \lambda_{\pm} = \frac{1}{2\beta_{0}}\left\lbrack \gamma_{qq}^{0} + \gamma_{gg}^{0} \pm \sqrt{\left( \gamma_{qq}^{0} - \gamma_{gg}^{0} \right)^{2} + 4\gamma_{qg}^{0}\gamma_{gq}^{0}} \right\rbrack

and the corresponding projector matrices

.. math:: \mathbf{e}_\pm=\frac{\pm 1}{\lambda_+ - \lambda_-}(R^{(0)}-\lambda_\mp\mathbb{I}),

in the following form

.. math:: \mathbf{L}(a_s,a_0,N)= \mathbf{e}_-(\frac{a_s}{a_0})^{-\lambda_{-(N)}} + \mathbf{e}_+(\frac{a_s}{a_0})^{-\lambda_{+(N)}}.

We express the solution of the evolution equation as a perturbative expansion
around the LO solution :math:`\mathbf{L}(a_s,a_0,N)`

.. math:: \binom{\Sigma}{g}(N,a_s) = \bigg[\mathbb{I}+\sum_{k=1}^{\infty}a_s^kU_k(N)\bigg] \mathbf{L}(a_s,a_0,N)\bigg[\mathbb{I}+\sum_{k=1}^{\infty}a_0^kU_k(N)\bigg]^{-1}\binom{\Sigma}{g}(N,a_0)\equiv \mathbf{\Gamma}_S(N,a_s,a_0)\binom{\Sigma}{g}(N,a_0)

The *fully truncated*\  expression of the matrix evolution kernel up to NNLO reads

.. math:: \mathbf{\Gamma}_S(N) = \big[\mathbf{L} + a_s\mathbf{U}_1\mathbf{L} - a_0\mathbf{LU}_1 + a_s^2 \mathbf{U}_2\mathbf{L} - a_sa_0 \mathbf{U}_1\mathbf{LU}_1 + a_0^2\mathbf{L}(\mathbf{U}_1^2 - \mathbf{U}_2)\big].

The :math:`U` matrices introduced in the previous equation are defined by this commutation relations

.. math:: \big[ \mathbf{U}_1, \mathbf{R}_0 \big] =  \mathbf{R}_1 + \mathbf{R}_1

.. math:: \big[ \mathbf{R}_2, \mathbf{R}_0 \big]  =  \mathbf{R}_2 +\mathbf{R}_1 \mathbf{U}_1 + 2 \mathbf{U}_2

.. math:: \vdots

.. math:: \big[ \mathbf{U}_k, \mathbf{R}_0 \big] =  \mathbf{R}_k + \sum_{i=1}^{k-1} \mathbf{R}_{k-i} \mathbf{U}_i + k \mathbf{U}_k \equiv\ \widetilde{\mathbf{R}}_k + k \mathbf{U}_k.

as

.. math:: \mathbf{U}_k=-\frac{1}{k}[e_+\widetilde{\mathbf{R}}_ke_+ + e_-\widetilde{\mathbf{R}}_ke_-] + \frac{e_+ \widetilde{\mathbf{R}}_k e_-}{\lambda_- -\lambda_+ - k} + \frac{e_-\widetilde{\mathbf{R}}_ke_+}{\lambda_+ -\lambda_- - k}

where

.. math:: \widetilde{\mathbf{R}}_k = \mathbf{R}_k+\sum_{i=1}^{k-1}\mathbf{R}_{k-i}\mathbf{U}_i.

By solving recursively the above equations

.. math:: \mathbf{R}_0 \equiv \frac{\boldsymbol{\gamma}^{(0)}}{\beta_0}

.. math:: \mathbf{R}_k\equiv - b_1 \mathbf{R}_{k-1} + \mathcal{O}(\textrm{NNLO})

the NLO full solution (corresponding to IMODEV=1 in ref.) can be easily implemented into the code.
Practically the sum in the above equation is stopped to a sufficiently high order such as k=20.

Getting back the x-space PDF’s
------------------------------

The :math:`x` space parton distributions are obtained by taking the
inverse Mellin transforms of the solutions obtained in eq.
(`[eq:solutionexpand] <#eq:solutionexpand>`__) which, making use of the
convolution theorem, can be written as

.. math::

   \begin{matrix}
   q_{NS}^{\pm ,v}(x,Q^{2}) & = & \int_{x}^{1}\frac{dy}{y}\Gamma_{qq}(y,a_{s},a_{0})\, q_{NS}^{\pm ,v}\left( \frac{x}{y},Q_{0}^{2} \right) \\
   \begin{pmatrix}
   \Sigma \\
   g \\
   \end{pmatrix}(x,Q^{2}) & = & \int_{x}^{1}\frac{dy}{y}\Gamma_{S}(y,a_{s},a_{0})\begin{pmatrix}
   \Sigma \\
   g \\
   \end{pmatrix}\left( \frac{x}{y},Q_{0}^{2} \right) \\
   \end{matrix}

The evolution kernels :math:`\Gamma(x)` are defined as the inverse
Mellin transforms of the evolution factors introduced in eqs.
(`[eq:solutionexpand] <#eq:solutionexpand>`__)

.. math:: \Gamma_{S}(x,a_{s},a_{0}) = \int_{c - i\infty}^{c_{+}i\infty}\frac{dN}{2\pi i}x^{- N}\Gamma_{S}(N,a_{s},a_{0})

Note however that all splitting functions, except the off-diagonal
entries of the singlet matrix, diverge when :math:`x = 1`, this implies
that the evolution kernels :math:`\Gamma(x)` will likewise be divergent
in :math:`x = 1`.

We now show that, like the splitting functions, the evolution factors
can be defined as distributions. To this purpose consider the generic
evolution factor :math:`\Gamma` such that (omitting the explicit
dependence of :math:`\Gamma` on the coupling :math:`a_{s}`)

.. math:: f(x,Q^{2}) = \int_{x}^{1}\frac{dy}{y}\Gamma(y)f\left( \frac{x}{y},Q_{0}^{2} \right)\,.

Defining the distribution

.. math:: \Gamma_{+}(x) = \Gamma(x) - \gamma\delta(1 - x)\,,\text{\quad\quad}where\quad\gamma = \int_{0}^{1}dx\Gamma(x)\,.

Equation (`[eq:gengamma] <#eq:gengamma>`__) can then be rewritten as

.. math::  f(x,Q^{2}) & = \gamma f(x,Q_{0}^{2}) + \int_{x}^{1}\frac{dy}{y}\Gamma_{+}(y)f\left( \frac{x}{y},Q_{0}^{2} \right) \\
    & = \gamma f(x,Q_{0}^{2}) + \int_{x}^{1}\frac{dy}{y}\Gamma(y)\left\lbrack f\left( \frac{x}{y},Q_{0}^{2} \right) - yf\left( x,Q_{0}^{2} \right) \right\rbrack - f(x,Q_{0}^{2})\int_{0}^{x}dy\Gamma(y)\,. \\

Due to the subtraction eq. `[eq:gammadist] <#eq:gammadist>`__, all
integrals on the r.h.s of eq. `[eq:genexp] <#eq:genexp>`__ converge and
can be evaluated numerically. We can then use this expression to compute
the parton distribution functions in :math:`x` space, determining
:math:`\Gamma` numerically from eq.\ `[eq:xkernels] <#eq:xkernels>`__
and :math:`\gamma` as

.. math:: \gamma = \int_{0}^{1}dx\int_{c - i\infty}^{c + i\infty}\frac{dN}{2\pi i}x^{- N}\Gamma(N) = \int_{c - i\infty}^{c + i\infty}\frac{dN}{2\pi i}\frac{\Gamma(N)}{1 - N}\,.

In this singlet case, however this prescription has been slightly
modified because :math:`\Gamma(N)|_{N = 1}` is indeed infinite. So
eq.\ `[eq:genexp] <#eq:genexp>`__ is rewritten in another equivalent
form. Let us define

.. math:: f^{(1)}(x,Q^{2}) = x f(x,Q^{2}) .

.. math:: \Gamma^{(1)}(x,Q_{0}^{2},Q^{2}) = x\Gamma(x,Q_{0}^{2},Q^{2}) .

Thus

.. math::

   \begin{matrix}
   f^{(1)}(x,Q^{2}) & = & x\, f(x,Q^{2}) = \int_{x}^{1}\,\frac{dy}{y}\,\Gamma(y,Q_{0}^{2},Q^{2})\, x\, f\left( \frac{x}{y},Q_{0}^{2} \right) \\
    & = & \int_{x}^{1}\,\frac{dy}{y}\,\Gamma^{(1)}(y,Q_{0}^{2},Q^{2})\, f^{(1)}\left( \frac{x}{y},Q_{0}^{2} \right) \\
    & = & \int_{x}^{1}\,\frac{dy}{y}\,\Gamma^{(1)}(y,Q_{0}^{2},Q^{2})\,\left( f^{(1)}\left( \frac{x}{y},Q_{0}^{2} \right) - yf^{(1)}(x,Q_{0}^{2}) \right) \\
    & + & \int_{x}^{1}\,\frac{dy}{y}\, y\Gamma^{(1)}(y,Q_{0}^{2},Q^{2})\, f^{(1)}(x,Q_{0}^{2}) \\
    & = & \int_{x}^{1}\,\frac{dy}{y}\,\Gamma^{(1)}(y,Q_{0}^{2},Q^{2})\,\left( f^{(1)}\left( \frac{x}{y},Q_{0}^{2} \right) - yf^{(1)}(x,Q_{0}^{2}) \right) \\
    & + & f^{(1)}(x,Q_{0}^{2})\left\lbrack \int_{0}^{1}\, dy\, y\Gamma(y,Q_{0}^{2},Q^{2}) - \int_{0}^{x}\, y\Gamma(y) \right\rbrack \\
    \Rightarrow f(x,Q^{2}) & = & \int_{x}^{1}\,\frac{dy}{y}\, y\Gamma(y,Q_{0}^{2},Q^{2})\,\left( \frac{1}{y}f\left( \frac{x}{y},Q_{0}^{2} \right) - yf(x,Q_{0}^{2}) \right) \\
    & + & f(x,Q_{0}^{2})\left\lbrack \Gamma(N,Q_{0}^{2},Q^{2})|_{N = 2} - \int_{0}^{x}\, y\Gamma(y,Q_{0}^{2},Q^{2}) \right\rbrack \\
   \end{matrix} .
