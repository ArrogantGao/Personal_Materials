#import "@preview/touying-flow:1.1.0":*
#show: flow-theme.with(
  aspect-ratio: "4-3",
  footer: self => self.info.title,
  footer-alt: self => self.info.subtitle,
  navigation: "none",
  primary:rgb("#004098"),//rgb(0,108,57),//rgb("#006c39"),
  secondary:rgb("#004098"),//rgb(161,63,61),//rgb("#a13f3d"),
  text-font: ("Libertinus Serif"),
  text-size: 20pt,
  // code-font: ("Jetbrains Mono NL","PingFang SC"),
  // code-size: 16pt,

  config-info(
    title: [Numerical Methods for Modeling and Simulation],
    subtitle: [Fast Summation Algorithms and Tensor Network Methods],
    author: text(23pt)[Xuanzhao Gao],
    institution: text(20pt)[Hong Kong University of Science and Technology],
    date: text(23pt)[2024-12-26],
  ),
  // config-common(show-notes-on-second-screen: right),
)

#title-slide()

// == About Me

// My name is Xuanzhao Gao, a fourth year PhD student at HKUST, supervised by Prof. Zecheng Gan.

= Fast Summation Algorithms

== Background

=== Molecular dynamics simulation

Molecular dynamics simulation is a widely used method to study the motion of atoms and molecules.

#align(center,
  grid(columns: 3, 
    image("figs/MD.png", width: 400pt),
    h(50pt),
    image("figs/MD_workflow.png", width: 270pt),
  )
)

#pagebreak()

=== Doubly periodic systems

Quasi-2D systems are at the macroscopic scale in $x y$, but microscopic in $z$, so that are always modeled as doubly periodic in simulations.
Q2D systems are widely exist in nature and engineering, for example, cell membranes and electrolyte near surfaces.

#figure(
  image("figs/Q2D.png", width: 600pt),
  caption: [Illustration of a Quasi-2D system.],
)

#pagebreak()
=== Long range Coulomb interaction and Ewald summation

Coulomb interaction plays a key role in Quasi-2D systems, leading to effect such as ion transport and self-assembly.

#figure(
  image("figs/self-assembly.png", width: 450pt),
  caption: [Self-assembly of nanoparticles in a Quasi-2D system #footnote(text(12pt)[Barros, Kipton, and Erik Luijten. Phys. Rev. Lett. 113, 017801],).
],)

However, the Coulomb interaction decays as $r^(-1)$ in 3D, so that is long ranged and singular at $r=0$, which make such simulation computationally expensive.
Complexity of a direct sum of the Coulomb interaction in doubly periodic systems is about $O(N^2 epsilon^(-1/3))$.

#pagebreak()

=== Long range Coulomb interaction and Ewald summation

A traditional method to deal with the long range Coulomb interaction is Ewald summation, where the Coulomb kernel is split into short-range and long-range parts:
$
  1 / r = (#text[erfc($alpha r$)]) / r + (#text[erf($alpha r$)]) / r
$
where $#text[erfc($x$)]$ is the complementary error function and $#text[erf($x$)]$ is the error function.

#align(center,
  image("figs/erfc.png", width: 350pt),
)


#pagebreak()

=== Long range Coulomb interaction and Ewald summation

The short-range part then is truncated in real space, and the long-range part is truncated in Fourier space.

Due to the periodicity in $x$ and $y$, the long-range interaction energy can be given by a 
$
  E_("long") = 1/(L_x L_y) sum_(i, j) q_i q_j sum_(k_x, k_y) integral_(- infinity)^(infinity) e^(-i arrow(k) dot arrow(r)_(i j))/(k^2) e^(-k^2/(4 alpha^2)) d k_z - alpha / sqrt(pi) sum_i q_i^2
$
where the infinite integral can be calculated analytically:
$
  integral_(- infinity)^(infinity) e^(-i arrow(k) dot arrow(r))/(k^2) e^(-k^2/(4 alpha^2)) d k_z = pi / (h) [e^(k_z z_(i j))"erfc"(h/(2 alpha) + alpha z_(i j)) + e^(-k_z z_(i j))"erfc"(h/(2 alpha) - alpha z_(i j))]
$
where $h = sqrt(k_x^2 + k_y^2)$. This function decay exponentially fast in the reciprocal space.

The resulting algorithm is called Ewald2D#footnote(text(12pt)[D. Parry, Surf. Sci. 49 (2) (1975) 433–440.],), with a complexity of $O(N^2 log epsilon)$, which prohibits the simulation of large systems.


== Sum-of-Exponential Ewald2D Method

Inspired by the fast Gaussian transform via sum-of-exponential (SOE) approximation#footnote(text(12pt)[S. Jiang, L. Greengard, Commun. Comput. Phys. 31 (1) (2022) 1–26],), we propose a sum-of-exponential Ewald2D method utilizing the state-of-the-art SOE approximation#footnote(text(12pt)[Z. Gao, J. Liang, Z. Xu, J. Sci. Comput. 93 (2) (2022) 40],).


== Fast Spectral Sum-of-Gaussian Method

Introduce the fast spectral sum-of-Gaussian method.

= Tensor Network Based Algorithms

== Background

=== Tensor Networks and Contraction Order Optimization

Introduce the tensor networks and contraction order optimization.

== Probabilistic Inference via Tensor Networks and Differential Programming

Tensor inference.

== Automatic Discovery of the Optimal Branching Rules

Optimal branching for combinatorial optimization.

= Future Research Plans

== Fast Algorithms Based on the DMK Framework

Fast algorithms on doubly periodic Coulomb systems based on the DMK framework.

== Boundary Integral Equations

For sharp boundaries.

== Quantum Many-Body Systems

Using general tensor networks as ansatz to represent quantum many-body systems.

= Appendix

== Details of the SoEwald2D method

