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

#pagebreak()

=== Doubly periodic Coulomb systems

Introduce the doubly periodic Coulomb systems.

== Sum-of-Exponential Ewald2D Method

Introduce the sum-of-exponential Ewald2D method.

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
