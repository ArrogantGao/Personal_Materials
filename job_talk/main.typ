#import "@preview/cetz:0.2.2": canvas, draw, tree, plot
#import "graph.typ": *
#import "@preview/subpar:0.2.0"

#set page(height: auto)
#set par(justify: true)

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

#show link: underline

#title-slide()

// == About Me

// My name is Xuanzhao Gao, a fourth year PhD student at HKUST, supervised by Prof. Zecheng Gan.

= Fast Summation Algorithms

== Quasi-2D Systems and Long Range Coulomb Interaction

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
  caption: [Illustration of a doubly periodic system.],
)

#pagebreak()
=== Coulomb interaction

Coulomb interaction plays a key role in nature, leading to effect such as ion transport and self-assembly.

#figure(
  image("figs/self-assembly.png", width: 450pt),
  caption: [Self-assembly of nanoparticles due to Coulomb interaction#footnote(text(12pt)[Barros, Kipton, and Erik Luijten. Phys. Rev. Lett. 113, 017801],).
],)

However, the Coulomb interaction decays as $r^(-1)$ in 3D, so that is long ranged and singular at $r=0$, which make such simulation computationally expensive.
Complexity of a direct sum of the Coulomb interaction in doubly periodic systems is about $O(N^2 epsilon^(-1/3))$.

#pagebreak()

=== Ewald summation

A traditional method to deal with the long range Coulomb interaction is Ewald summation, where the Coulomb kernel is split into short-range and long-range parts:
$
  1 / r = (#text[erfc($alpha r$)]) / r + (#text[erf($alpha r$)]) / r
$
where $#text[erfc($x$)]$ is the complementary error function and $#text[erf($x$)]$ is the error function.

#align(center,
  image("figs/erfc.png", width: 350pt),
)


#pagebreak()

=== Ewald summation

The short-range part then is truncated in real space, and the long-range part is truncated in Fourier space.

Due to the periodicity in $x$ and $y$, the long-range interaction energy can be given by a 
$
  E_("long") = 1/(L_x L_y) sum_(i, j) q_i q_j sum_(k_x, k_y) integral_(- infinity)^(infinity) e^(-i arrow(k) dot arrow(r)_(i j))/(k^2) e^(-k^2/(4 alpha^2)) d k_z - alpha / sqrt(pi) sum_i q_i^2
$
where the infinite integral can be calculated analytically:
$
  h / pi integral_(- infinity)^(infinity) e^(-i k_z z_(i j))/(k^2) e^(-k_z^2/(4 alpha^2)) d k_z = underbrace(e^(h abs(z_(i j)))"erfc"(h/(2 alpha) + alpha abs(z_(i j))), xi^+(z_(i j))) + underbrace(e^(-h abs(z_(i j)))"erfc"(h/(2 alpha) - alpha abs(z_(i j))), xi^-(z_(i j)))
$
where $h = sqrt(k_x^2 + k_y^2)$. This function decay exponentially fast in the reciprocal space.

The resulting algorithm is called Ewald2D#footnote(text(12pt)[D. Parry, Surf. Sci. 49 (2) (1975) 433–440.],), with a complexity of $O(N^2 log epsilon)$.
// , which prohibits the simulation of large systems.


== Sum-of-Exponential Ewald2D Method

=== The SOE approximation

Inspired by the fast Gaussian transform via sum-of-exponential (SOE) approximation#footnote(text(12pt)[S. Jiang, L. Greengard, Commun. Comput. Phys. 31 (1) (2022) 1–26],), we propose a sum-of-exponential Ewald2D method#footnote(text(12pt)[Z. Gan, X. Gao, J. Liang, and Z. Xu. arXiv preprint arXiv:2405.06333, 2024; also see #link("https://github.com/HPMolSim/SoEwald2D.jl")],) utilizing VPMR#footnote(text(12pt)[Z. Gao, J. Liang, Z. Xu, J. Sci. Comput. 93 (2) (2022) 40],) method.

The Guassian kernel is approximated by a sum of exponential functions as follows:
$
  e^(- alpha^2 t^2) approx sum_(l = 1)^M  w_l e^(- s_l alpha abs(t))
$
where $w_l$ and $s_l$ are complex parameters of the SOE approximation, and the following error bound is satisfied:
$
  abs( e^(- alpha^2 t^2) - sum_(l = 1)^M  w_l e^(- s_l alpha abs(t)) ) < epsilon.
$
In VPMR, $M = 4$, $8$ and $16$ gives the error bound $epsilon = 10^(-4)$, $10^(-8)$ and $10^(-14)$, respectively.

#pagebreak()

=== Approximation of the Ewald2D kernel

In Ewald2D summation, we observed the following identity:
$
  xi^+ (z) = e^( h abs(z))"erfc"(h/(2 alpha) + alpha abs(z)) & = (2 alpha) / sqrt(pi) e^(-h^2/(4 alpha^2)) e^(h abs(z)) integral_(abs(z))^(infinity) e^(-alpha^2 t^2 - h t) d t\
$
and the SOE approximation is applied to get:
$
   xi^+ (z) approx (2 alpha) / sqrt(pi) e^(-h^2/(4 alpha^2)) e^(h abs(z)) integral_(abs(z))^(infinity) sum_(l = 1)^M  w_l e^(- s_l alpha abs(t)  - h t) d t = (2 alpha) / sqrt(pi) e^(-h^2/(4 alpha^2)) sum_(l = 1)^M  w_l e^(- s_l alpha abs(z)) / (alpha s_l + h)
$
For $xi^-(z)$, the same approximation is applied.

Now the double summation over $i$ and $j$ in energy is approximated by
$
  sum_(i, j) q_i q_j xi (z_(i j)) approx sum_(l = 1)^M C_l sum_(i, j) q_i q_j e^(- s_l alpha abs(z_(i j)))
$
where $C_l$ is a constant irrelevant to $i$ and $j$.

#pagebreak()

=== Summing up the exponentials via sorting

The double summation can be simplified as 
$
  S = sum_(i, j) q_i q_j e^(- abs(z_i - z_j))
$
which still takes $O(N^2)$ operations due to the absolute value.

For further acceleration, we reorder the indices via sorting (at most $O(N log N)$ operations):
$
  z_1 < z_2 < ... < z_N
$
so that the absolute value can be removed: 
$
  S = sum_(i = 1)^N q_i e^(-z_i) underbrace(sum_(j = 1)^i q_j e^(z_j), A_i) = sum_(i = 1)^N q_i e^(-z_i) A_i
$
where the array $A$ can be computed iteratively in $O(N)$ operations, then calculating $S$ takes another $O(N)$ operations.

#pagebreak()

=== Complexity of the SOEwald2D method

Utilizing the SOE approximation, the double summation can be calculated in $O(N)$ operations.
Calculation of the energy can be simplified as
$
  sum_(k_x, k_y) e^(-h^2 / (4 alpha^2)) sum_(i, j) S(k_x, k_y, r_i, r_j) = underbrace(sum_(k_x, k_y), O(N^0.4)) e^(-h^2 / (4 alpha^2)) underbrace(f(k_x, k_y), O(N))
$
However, a summation in the reciprocal space is still required, and the resulting total complexity is $O(N^(1.4))$ since number of the Fourier mode to be summed is of $O(N^(0.4))$.

#pagebreak()

=== Random batch sampling

To further reduce the complexity, we use the the random batch sampling technique#footnote(text(12pt)[S. Jin, L. Li, Z. Xu, Y. Zhao, SIAM J. Sci. Comput. 43 (4) (2021) B937–B960.],), where a importance sampling is applied:
$
  sum_(k_x, k_y) e^(-h^2 / (4 alpha^2)) f(k_x, k_y) approx H / P sum_((k_x, k_y) in cal(K)_P) f(k_x, k_y)
$
where $P$ Fourier modes are selected using the Gaussian as distribution, $H$ is the summation of Guassians.

It has been shown that for a system with fixed density, $P ~ O(1)$ is sufficient to achieve the correct equilibrium state, and gives an accurate ensemble average.

The resulting algorithm is called RBSE2D, with a complexity of $O(N)$.


#pagebreak()

=== Numerical Results

The accuracy of the SOEwald2D method is determined by the error bound of the SOE approximation, and coverage exponentially fast as the Ewald2D method.

#figure(
  image("figs/soewald2d_error.png", width: 600pt),
  caption: [Error of the SOEwald2D method (a) with different number of term $M$ used in the SOE approximation (b) with different system sizes ($L_z$ is fixed).],
)

#pagebreak()

=== Numerical Results

We use the RBSE2D method to simulate a typical Q2D electrolyte system as shown below:

#figure(
  image("figs/soewald2d_md.png", width: 500pt),
  caption: [The charge density of the Q2D electrolyte system in $z$ direction with different $P$, inset shows the relative error of the average energy.],
)

#pagebreak()

=== Numerical Results

Time complexity of the method is verified as below:

#figure(
  image("figs/soewald2d_complexity.png", width: 500pt),
  caption: [Complexity of the Ewald2D method, SOEwald2D method and RBSE2D method.],
)

#pagebreak()

=== Summary

The SoEwald2D utilized the SOE approximation and spirit of the FGT to accelerate the Ewald2D method.

It has the following advantages:
- Reducing the complexity from $O(N^2)$ to $O(N^1.4)$ with a controlled error bound.
- With the random batch sampling technique, the complexity can be further reduced to $O(N)$.
- Not sensitive to the aspect ratio of the system.

However, it also has some disadvantages:
- A sorting step is required, which is not efficient for parallel computing.
- Importance sampling is necessary to achieve linear complexity, which may not be efficient for cases requiring exact results.

== Fast Spectral Sum-of-Gaussian Method

=== Sum-of-Gaussian approximation of Coulomb kernel



= Tensor Network Based Algorithms

== Tensor Networks and Their Contraction Order

=== Tensor Networks

Tensor network (TN) is a powerful tool to represent and manipulate high-dimensional data, and have been used in various fields such as quantum physics, machine learning and combinatorial optimization.

#subpar.grid(
  figure(image("figs/tn_manybody.png"), caption: [
    Tensor networks for quantum many-body systems#footnote(text(12pt)[#link("https://tensornetwork.org/")],).
  ]), <a>,
  figure(image("figs/tn_inference.png"), caption: [
    Tensor network for Bayesian inference.
  ]), <b>,
  columns: (350pt, 240pt),
  label: <full>,
  caption: [Examples of tensor networks.],
)

#pagebreak()

=== Einsum notation

Tensor networks can be represented as the so called Einsum notation:
$
  Y_(i_y ...) = sum_(i in.not {i_y ...}) A_(i_a ...) B_(i_b ...) C_(i_c ...) ...
$

It also has a hyper-graph representation, where each node is a tensor and each edge is an index.
For example, the following Einsum expression can be represented as the following hyper-graph:
#align(center,
// grid(columns: 3,
// align(horizon,
// $
//   
// $),
// h(50pt),
align(center, canvas(length:40pt, {
  import draw: *
  content((-7, 1.2), text(20pt)[$s = sum_(i, j, k, l) A_(i j) B_(j k) C_(k l) D_(l i) E_(i) F_(j) G_(k) H_(l)$])

  content((-2, 1.3), text(30pt)[$arrow$])

  // petersen graph
  let vertices1 = ((0, 0), (5, 0), (5, 5), (0, 5))
  let edges = ((0, 1), (1, 2), (2, 3), (3, 0))

  show-graph((vertices1).map(v=>(0.5 * (v.at(0)), 0.5 * (v.at(1)))), edges, radius:0.0)

  let vertices3 = ((0, 0), (5, 0), (5, 5), (0, 5), (-1.5, -1.5), (6.5, -1.5), (6.5, 6.5), (-1.5, 6.5))
  let edges3 = ((0, 4), (1, 5), (2, 6), (3, 7))
  show-graph((vertices3).map(v=>(0.5 * (v.at(0)), 0.5 * (v.at(1)))), edges3, radius:0.0)
  // 
  let indices = ("i", "l", "k", "j")
  for (i, v) in (vertices1.map(v=>(0.5 * (v.at(0)), 0.5 * (v.at(1))))).enumerate() {
    circle(v, radius:0.3, fill:white, stroke:none)
    content(v, text(20pt, black)[#(indices.at(i))])
  }

  let vertices2 = ((2.5, 0), (5, 2.5), (2.5, 5), (0, 2.5), (-1.5, -1.5), (6.5, -1.5), (6.5, 6.5), (-1.5, 6.5))
  let edges2 = ()
  show-graph-content((vertices2).map(v=>(0.5 * (v.at(0)), 0.5 * (v.at(1)))), edges2, ((0, "D"), (1, "C"), (2, "B"), (3, "A"), (4, "E"), (5, "H"), (6, "G"), (7, "F")), radius:0.4, fontsize:20pt)
}))
)

#pagebreak()

=== Contraction order

To extract information from the tensor network, we need to contract them, i.e. sum over the indices.
A naive way is to loop over all indices directly. However, in this way the cost grows exponentially with the number of indices.

In practice, we prefer to do binary contraction, i.e. contracting two tensors at a time so that BLAS can be utilized. The resulting order can be represented as a binary tree.
The contraction order is only related to the structure of the tensor network, and is independent of the data.

#figure(
align(center, canvas(length:40pt, {
  import draw: *
  set-origin((4, 0.35))
  let DY = 1
  let DX1 = 2.5
  let DX2 = 1.2
  let DX3 = 0.7
  let root = (0, DY)
  let left = (-DX1, 0)
  let right = (DX1, 0)
  let left_left = (-DX1 - DX2, -DY)
  let left_right = (-DX1 + DX2, -DY)
  let right_left = (DX1 - DX2, -DY)
  let right_right = (DX1 + DX2, -DY)

  let left_left_left = (-DX1 - DX2 - DX3, -2*DY)
  let left_left_right = (-DX1 - DX2 + DX3, -2*DY)
  let left_right_left = (-DX1 + DX2 - DX3, -2*DY)
  let left_right_right = (-DX1 + DX2 + DX3, -2*DY)
  let right_left_left = (DX1 - DX2 - DX3, -2*DY)
  let right_left_right = (DX1 - DX2 + DX3, -2*DY)
  let right_right_left = (DX1 + DX2 - DX3, -2*DY)
  let right_right_right = (DX1 + DX2 + DX3, -2*DY)

  for (a, b) in ((root, left), (root, right), (left, left_left), (left, left_right), (right, right_left), (right, right_right), (left_left, left_left_left), (left_left, left_left_right), (left_right, left_right_left), (left_right, left_right_right), (right_left, right_left_left), (right_left, right_left_right), (right_right, right_right_left), (right_right, right_right_right)){
    line(a, b)
  }

  for (l, t) in ((left_left_left, [$A_(i j)$]), (left_left_right, [$E_i$]), (left_right_left, [$B_(j k)$]), (left_right_right, [$F_k$]), (right_left_left, [$C_(k l)$]), (right_left_right, [$G_k$]), (right_right_left, [$D_(l i)$]), (right_right_right, [$H_l$])){
    circle(l, radius:0.4, fill: black)
    content(l, text(15pt, white)[#t])
  }

})),
caption: [An example of binary contraction tree.],
)

Different contraction orders can lead to different complexities, the order with the minimum complexity is called the optimal contraction order.

#pagebreak()

=== Optimizing the contraction order

Finding the optimal contraction order is a NP-hard problem.
However, it is lucky that we are always satisfied with a good enough solution, thus heuristic methods are sufficient in practice.

In the past few years, tools have been developed to optimize the contraction order, such as the OMEinsumContractionOrder.jl#footnote(text(12pt)[#link("https://github.com/TensorBFS/OMEinsumContractionOrder.jl")],) in Julia and cotengra#footnote(text(12pt)[#link("https://github.com/jcmgray/cotengra")],) in Python.

#figure(
  image("figs/tn_order.png", width: 280pt),
  caption: [Comparison of different contraction orders.],
)

== Probabilistic Inference via Tensor Network Contraction

=== Probabilistic graphical models

The probabilistic graphical models (PGM) is a probabilistic model for which a graph expresses the conditional dependence structure between random variables.

A undirected graph can be used to factorize the joint distribution of the random variables, and can be represented as a tensor network.

For example, given the following joint distribution:
$
  P(i_0, j_0, k_0, l_0) = P_(i j) (i_0, j_0) P_(j k) (j_0, k_0) P_(k l) (k_0, l_0) P_(l i) (l_0, i_0)
$
it can be represented as a tensor network with open edges:

#align(center, canvas(length:40pt, {
  import draw: *
  // petersen graph
  let vertices1 = ((0, 0), (5, 0), (5, 5), (0, 5))
  let edges = ((0, 1), (1, 2), (2, 3), (3, 0))

  show-graph((vertices1).map(v=>(0.5 * (v.at(0)), 0.5 * (v.at(1)))), edges, radius:0.0)

  let vertices3 = ((0, 0), (5, 0), (5, 5), (0, 5), (-1.5, -1.5), (6.5, -1.5), (6.5, 6.5), (-1.5, 6.5))
  let edges3 = ((0, 4), (1, 5), (2, 6), (3, 7))
  show-graph((vertices3).map(v=>(0.5 * (v.at(0)), 0.5 * (v.at(1)))), edges3, radius:0.0)
  // 
  let indices = ("i", "l", "k", "j")
  for (i, v) in (vertices1.map(v=>(0.5 * (v.at(0)), 0.5 * (v.at(1))))).enumerate() {
    circle(v, radius:0.3, fill:white, stroke:none)
    content(v, text(20pt, black)[#(indices.at(i))])
  }

  let vertices2 = ((2.5, 0), (5, 2.5), (2.5, 5), (0, 2.5), )
  let edges2 = ()
  show-graph-content((vertices2).map(v=>(0.5 * (v.at(0)), 0.5 * (v.at(1)))), edges2, ((0, $P_(i l)$), (1, $P_(l k)$), (2, $P_(j k)$), (3, $P_(i j)$)), radius:0.4, fontsize:15pt)
}))

#pagebreak()

=== Tasks for inference

The following tasks are commonly required in probabilistic inference:
- PR: the partition function
- MAR: the marginal probability distribution over sets of variables
- MPE: the most probable explanation (most likely assignment) to all variables of the model
- MMAP: the maximum marginal a posteriori, i.e., the most likely assignment to a set of variables after marginalizing a different set

In the following, we will show how to perform these tasks via the modern tensor network techniques in our work#footnote(text(12pt)[M. Roa-Villescas, X. Gao, S. Stuijk, H. Corporaal, and J.-G. Liu, Phys. Rev. Research 6, 033261 (2024). Also see #link("https://github.com/TensorBFS/TensorInference.jl")],).

#pagebreak()

=== Partition function

Among these tasks, the PR task can be directly computed by contracting the network with an all-one tensor attached to each indices.

#align(center, canvas(length:40pt, {
  import draw: *
  // petersen graph
  let vertices1 = ((0, 0), (5, 0), (5, 5), (0, 5))
  let edges = ((0, 1), (1, 2), (2, 3), (3, 0))

  show-graph((vertices1).map(v=>(0.5 * (v.at(0)), 0.5 * (v.at(1)))), edges, radius:0.0)

  let vertices3 = ((0, 0), (5, 0), (5, 5), (0, 5), (-1.5, -1.5), (6.5, -1.5), (6.5, 6.5), (-1.5, 6.5))
  let edges3 = ((0, 4), (1, 5), (2, 6), (3, 7))
  show-graph((vertices3).map(v=>(0.5 * (v.at(0)), 0.5 * (v.at(1)))), edges3, radius:0.0)
  // 
  let indices = ("i", "l", "k", "j")
  for (i, v) in (vertices1.map(v=>(0.5 * (v.at(0)), 0.5 * (v.at(1))))).enumerate() {
    circle(v, radius:0.3, fill:white, stroke:none)
    content(v, text(15pt, black)[#(indices.at(i))])
  }

  let vertices2 = ((2.5, 0), (5, 2.5), (2.5, 5), (0, 2.5), (-1.5, -1.5), (6.5, -1.5), (6.5, 6.5), (-1.5, 6.5))
  let edges2 = ()
  show-graph-content((vertices2).map(v=>(0.5 * (v.at(0)), 0.5 * (v.at(1)))), edges2, ((0, $P_(i l)$), (1, $P_(l k)$), (2, $P_(j k)$), (3, $P_(i j)$), (4, $bb(1)$), (5, $bb(1)$), (6, $bb(1)$), (7, $bb(1)$)), radius:0.4, fontsize:15pt)
}))

For the complex network (more than 1000 random variables), the contraction order is optimized automatically by the OMEinsumContractionOrder.jl package via the heuristic local search solvers.

#pagebreak()

=== Marginal probability by differential programming

A direct forward way to compute the marginal probability, for example, $P_i$, is to leave the $i$ index open and contract the rest of the network. However, one have to repeat $N$ times to get the marginal probability for all variables.

Instead, we use the backward mode automatic differentiation, since:
$
  Z = sum_(i_0) Z_(i = i_0) times bb(1)_(i_0), #text("then") frac(partial Z, partial bb(1)_(i_0)) = Z(i = i_0)
$
which gives the marginal probability together with the partition: $P_(i = i_0) = Z(i = i_0) / Z$.

With the backward mode automatic differentiation, the marginal probability of all variables can be computed in a single contraction and a backward pass.

#pagebreak()

=== Most probable explanation

The MPE task is to find a most likely assignment to all variables so that the probability is maximized:
$
  arg max_(i_0, j_0, ...) P(i = i_0, j = j_0, ...)
$
Solving the maximum $P$ is equivalent to contracting the logarithm of the network under a tropical semiring, where
$
  a times.circle b &= a + b \ a plus.circle b &= max(a, b)
$
The contraction order is exactly the same as the one in MAR.

Utilizing a similar backpropagation method#footnote(text(12pt)[J. G. Liu, X. Gao, M. Cain, M. D. Lukin, and S. T. Wang, SIAM J. Sci. Comput. 45, A1239 (2023).],) as employed in MAR, one can derive the gradient with respect to the all-one tensor. 
In each gradient computation, there exists precisely one non-zero element, representing the most likely assignment.

#pagebreak()

== Automatic Discovery of the Optimal Branching Rules

Optimal branching for combinatorial optimization.

= Future Research Plans

== Fast Algorithms Based on the DMK Framework

Fast algorithms on doubly periodic Coulomb systems based on the DMK framework.

== Boundary Integral Equations

For sharp boundaries.

== Quantum Many-Body Systems

Using general tensor networks as ansatz to represent quantum many-body systems.

= Summary

#pagebreak()

A summary of the talk.

= Acknowledgements

#pagebreak()

#align(center,
grid(columns: 5, 
  grid(rows : 3, image("photos/zechenggan.jpeg", width : 150pt), "" ,  text[Prof. Zecheng Gan \ HKUST(GZ)]),
  h(50pt), 
  grid(rows : 3, image("photos/zhenlixu.jpg", width : 150pt), "" ,  text[Prof. Zhenli Xu \ SJTU]),
  h(50pt), 
  grid(rows : 3, image("photos/shidongjiang.jpeg", width : 150pt), "" ,  text[Prof. Shidong Jiang \ Flatiron Institute]),
  )
)

// #pagebreak()

// \

#align(center,
grid(columns: 3, 
  grid(rows : 3, image("photos/jinguoliu.jpeg", width : 150pt), "" ,  text[Prof. Jinguo Liu, \ HKUST(GZ)]),
  h(60pt), 
  grid(rows : 3, image("photos/panzhang.jpeg", width : 150pt), "" ,  text[Prof. Pan Zhang, \ ITP, CAS]),
  )
)

= Appendix

== Details of the SoEwald2D method

