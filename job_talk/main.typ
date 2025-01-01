#import "@preview/cetz:0.2.2": canvas, draw, tree, plot
#import "graph.typ": *
#import "@preview/subpar:0.2.0"

#import "@preview/pinit:0.1.3": *
#import "@preview/ctheorems:1.1.3": *
#import "diagbox.typ": *

// Theorems configuration by ctheorems
#show: thmrules.with(qed-symbol: $square$)

#let theorem = thmbox("theorem", "Theorem", fill: rgb("#eeffee")).with(numbering: none)
#let corollary = thmplain(
  "corollary",
  "Corollary",
  base: "theorem",
  titlefmt: strong
)

#let definition = thmbox("definition", "Definition", inset: (x: 1.2em, top: 1em)).with(numbering: none)
#let example = thmplain("example", "Example").with(numbering: none)
#let proof = thmproof("proof", "Proof")
#let hd(name) = table.cell(text(12pt)[#name], fill: green.lighten(50%))
#let s(name) = table.cell(text(12pt)[#name])

#set page(height: auto)
#set par(justify: true)

#let globalvars = state("t", 0)
#let timecounter(minutes) = [
  #globalvars.update(t => t + minutes)
  #place(dx: 100%, dy: 0%, align(right, text(16pt, red)[#context globalvars.get()min]))
]
#let clip(image, top: 0pt, bottom: 0pt, left: 0pt, right: 0pt) = {
  box(clip: true, image, inset: (top: -top, right: -right, left: -left, bottom: -bottom))
}

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
    title: [Modeling and Simulating Scientific Problems],
    subtitle: [Tensor Network Methods and Fast Summation Algorithms],
    author: text(23pt)[Xuanzhao Gao],
    institution: text(20pt)[Hong Kong University of Science and Technology],
    date: text(23pt)[2024-12-26],
  ),
  // config-common(show-notes-on-second-screen: right),
)

#show link: underline

#title-slide()

= Tensor Network Based Algorithms

== Tensor Networks and Exact Contraction

=== Tensor Networks

Tensor network (TN) is a powerful tool to represent and manipulate high-dimensional data, and have been used in various fields such as quantum physics and machine learning.

#subpar.grid(
  figure(image("figs/tn_manybody.png"), caption: [
    Tensor networks for quantum\ many-body systems#footnote(text(12pt)[#link("https://tensornetwork.org/")],).
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

It also has a hyper-graph representation, where each node is a tensor and each edge is an index:

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

// To extract information from the tensor network, we need to contract them, i.e. sum over the indices.
// A naive way is to loop over all indices directly. However, in this way the cost grows exponentially with the number of indices.

// In practice, we prefer to do binary contraction, i.e. contracting two tensors at a time so that BLAS can be utilized. The resulting order can be represented as a binary tree.
// The contraction order is only related to the structure of the tensor network, and is independent of the data.

A contraction order can be represented as a rooted (binary) tree:

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

Finding the optimal contraction order is a NP-hard#footnote(text(12pt)[I.L. Markov, Y. Shi, SIAM J. Comput. 38, 963–981 (2008).],) problem! 
// However, it is lucky that we are always satisfied with a good enough solution, thus heuristic methods are sufficient in practice.

In the past few years, tools have been developed to optimize the contraction order:
- OMEinsumContractionOrder.jl#footnote(text(12pt)[#link("https://github.com/TensorBFS/OMEinsumContractionOrder.jl")],) in Julia
- Cotengra#footnote(text(12pt)[#link("https://github.com/jcmgray/cotengra")],) in Python

#figure(
  image("figs/tn_order.png", width: 280pt),
  caption: [Comparison of different contraction orders.],
)

== Probabilistic Inference via Tensor Network Contraction

=== Probabilistic graphical models

The probabilistic graphical models (PGM) is a probabilistic model for which a graph expresses the conditional dependence structure between random variables.

// A undirected graph can be used to factorize the joint distribution of the random variables, and can be represented as a tensor network.

// A joint distribution can be represented as a tensor network with open edges.

// For example, given the following joint distribution:

#align(center, canvas({
  import draw: *
  let dx = 3
  let dy = 0.5
  for (x, y, name, txt) in ((-7, 0, "a", [Recent trip to #text(red)[A]sia]), (3.5, 0, "b", [Patient is a #text(red)[S]moker]), (-7, -2, "c", [#text(red)[T]uberculosis]), (0, -2, "d", [#text(red)[L]ung cancer]), (7, -2, "e", [#text(red)[B]ronchitis]), (-3.5, -4, "f", [#text(red)[E]ither T or L]), (-7, -6, "g", [#text(red)[X]-Ray is positive]), (0, -6, "h", [#text(red)[D]yspnoea])) {
    rect((x - dx, y - dy), (x + dx, y + dy), stroke: black, name: name, radius: 5pt)
    content((x, y), [#txt])
  }
  for (a, b) in (("a", "c"), ("b", "d"), ("b", "e"), ("c", "f"), ("d", "f"), ("f", "g"), ("f", "h"), ("e", "h")) {
    line(a, b, mark: (end: "straight"))
  }
  content((12, -3), box([*Tensors*\ p(A)\ p(S)\ p(T|A)\ p(L|S)\ p(B|S)\ p(E|T,L)\ p(X|B)\ p(D|E,X)], stroke: blue, inset: 10pt))
})) 

Marginal probability:

#text(19pt)[$P(L) = sum_(A, S, T, B, E, X, D) P(A) P(S) P(T|A) P(L|S) P(B|S) P(E|T,L) P(X|E) P(D|E,B)$]

#pagebreak()
=== TN representation of probabilistic models

A joint distribution can be factorized as:
$
  P(i_0, j_0, k_0, l_0) = P_(i j) (i_0, j_0) P_(j k) (j_0, k_0) P_(k l) (k_0, l_0) P_(l i) (l_0, i_0)
$
It can be represented as a tensor network with open edges:

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

#align(center,
table(
  columns: (auto, auto),
  inset: 10pt,
  align: horizon,
  table.header(
    [*Task*], [*Target*],
  ),
  [
    PR
  ],
  [
    partition function
  ],
  [
    MAR
  ],
  [
    marginal probability distribution
  ],
  [
    MPE
  ],
  [
    most probable explanation
  ],
  [
    MMAP
  ],
  [
    maximum marginal a posteriori
  ],
)
)

In the following, we will show how to perform these tasks via the modern tensor network techniques in our work#footnote(text(12pt)[M. Roa-Villescas, *X.-Z. Gao*, S. Stuijk, H. Corporaal, and J.-G. Liu, Phys. Rev. Research 6, 033261 (2024). Also see #link("https://github.com/TensorBFS/TensorInference.jl")],).

#pagebreak()

=== Solving the inference tasks

// Among these tasks, the PR task can be directly computed by contracting the network with an all-one tensor attached to each indices.
First add a all-one vector to each indices, then:

\

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

  let x = 9

  let center = (x, 4.2)
  content(center, text(20pt, black)[Contraction order optimization $arrow$ PR])
  let center = (x, 2.2)
  content(center, text(20pt, black)[Automatic differentiation + TN#footnote(text(12pt)[H. J. Liao, J. G. Liu, L. Wang, and T. Xiang, Phys. Rev. X 9, 031041 (2019).],) $arrow$ MAR])
  let center = (x, 0.1)
  content(center, text(20pt, black)[Tropical Tensor Network#footnote(text(12pt)[J. G. Liu, X. Gao, M. Cain, M. D. Lukin, and S. T. Wang, SIAM J. Sci. Comput. 45, A1239 (2023).],) + AD $arrow$ MPE])
  let center = (x, -2)
  content(center, text(20pt, black)[Mixed Tensor Network + AD $arrow$ MMAP])
}))

// For the complex network (more than 1000 random variables), the contraction order is optimized automatically by the OMEinsumContractionOrder.jl package via the heuristic local search solvers.


#pagebreak()

=== Numerical results on UAI problems

#figure(
  image("figs/ti_benchmark.png", width: 450pt),
  caption: [The benchmark of the inference tasks on UAI problems against the traditional solvers, an exponential speedup is achieved.],
)

#pagebreak()

=== GPU acceleration

The tropical semiring is not supported by the GPU BLAS libraries, thus we write a custom kernel#footnote(text(12pt)[#link("https://github.com/TensorBFS/CuTropicalGEMM.jl")],) with a C-CUDA backend.

#grid(columns: 3,
image("figs/cutropicalgemm.png", width: 350pt),
h(50pt),
align(horizon, text(20pt, black)[ Reaching $80\%$ of the peak performance on NVIDIA A800 (although fma is not supported) \ \ $4000$ times faster than the CPU version!]),
)

#pagebreak()

=== GPU acceleration

By replacing the matrix multiplication in the contraction with the GPU kernels, the performance is further improved.

#figure(
  image("figs/ti_gpu.png", width: 350pt),
  caption: [Speedup of the MMAP task on GPU against the CPU version.],
)

#pagebreak()

=== Summary

We applied the modern tensor network techniques to solve the probabilistic inference tasks, including contraction order optimization, automatic differentiation, generic tensor network, and GPU acceleration.
It is shown that with these techniques, an exponential speedup is achieved against the traditional solvers.

Advantages:
- powerful contraction order optimization tools
- automatic differentiation techniques
- easily accelerated by GPU

Disadvantages:
- exact inference can not handle ultra-large scale problems
- the contraction order optimization only uses the structural information


== Automatic Discovery of the Optimal Branching Rules

=== Failure of tensor network method

Similar techniques are also used in combinatorial optimization problems#footnote(text(12pt)[*X.-Z. Gao*, X.-F. Li, J.-G. Liu, arXiv:, (2024), also see #link("https://github.com/QuEraComputing/GenericTensorNetworks.jl")],) #footnote(text(12pt)[J.-G. Liu, X. Gao, M. Cain, M. D. Lukin, and S.-T. Wang, SIAM J. Sci. Comput., 45 (2023), pp. A1239–A1270],), such as the maximum independent set (MIS) problem, which can be solved by tropical TN.

MIS problem: given a graph, find the largest subset of vertices such that no two vertices are connected by an edge.

A simple problem: MIS on complete graph.

#align(center,
grid(columns: 3,
  canvas(length: 40pt, {
    import draw: *
    let N = 10
    let locs = ((1.0, 0.0), (0.8090169943749475, 0.5877852522924731), (0.30901699437494745, 0.9510565162951535), (-0.30901699437494734, 0.9510565162951536), (-0.8090169943749473, 0.5877852522924732), (-1.0, 1.2246467991473532e-16), (-0.8090169943749475, -0.587785252292473), (-0.30901699437494756, -0.9510565162951535), (0.30901699437494723, -0.9510565162951536), (0.8090169943749473, -0.5877852522924734))
    let e = ((1, 0), (2, 0), (3, 0), (4, 0), (5, 0), (6, 0), (7, 0), (8, 0), (9, 0), (0, 1), (2, 1), (3, 1), (4, 1), (5, 1), (6, 1), (7, 1), (8, 1), (9, 1), (0, 2), (1, 2), (3, 2), (4, 2), (5, 2), (6, 2), (7, 2), (8, 2), (9, 2), (0, 3), (1, 3), (2, 3), (4, 3), (5, 3), (6, 3), (7, 3), (8, 3), (9, 3), (0, 4), (1, 4), (2, 4), (3, 4), (5, 4), (6, 4), (7, 4), (8, 4), (9, 4), (0, 5), (1, 5), (2, 5), (3, 5), (4, 5), (6, 5), (7, 5), (8, 5), (9, 5), (0, 6), (1, 6), (2, 6), (3, 6), (4, 6), (5, 6), (7, 6), (8, 6), (9, 6), (0, 7), (1, 7), (2, 7), (3, 7), (4, 7), (5, 7), (6, 7), (8, 7), (9, 7), (0, 8), (1, 8), (2, 8), (3, 8), (4, 8), (5, 8), (6, 8), (7, 8), (9, 8), (0, 9), (1, 9), (2, 9), (3, 9), (4, 9), (5, 9), (6, 9), (7, 9), (8, 9))
    show-graph(locs.map(v => (v.at(0) * 2, v.at(1) * 2)), e, radius:0.1)
  }),
  h(60pt),
  align(center + horizon, text(20pt, black)[The solution is trivial. \ TN takes exponential time!]),
)
)

#align(center, text(20pt, black)[The reason is that the TN contraction does not use the information of context.])

#pagebreak()

=== The branching algorithm

Branching algorithm search the solution space in a tree-like structure. 

#figure(
canvas(length: 0.71cm, {
  import draw: *
  let scircle(loc, radius, name) = {
    circle((loc.at(0)-0.1, loc.at(1)-0.1), radius: radius, fill:black)
    circle(loc, radius: radius, stroke:black, name: name, fill:white)
  }
  let s = 1.5
  let dy = 3.0
  let la = (-s, 0)
  let lb = (0, s)
  let lc = (0, 0)
  let ld = (s, 0)
  let le = (s, s)
  scircle((0, 0), (3, 2), "branch")
  for (l, n) in ((la, "a"), (lb, "b"), (lc, "c"), (ld, "d"), (le, "e")){
    circle((l.at(0), l.at(1)-s/2), radius:0.4, name: n, stroke: if n == "a" {red} else {black})
    content((l.at(0), l.at(1)-s/2), text(14pt)[$#n$])
  }
  for (a, b) in (("a", "b"), ("b", "c"), ("c", "d"), ("d", "e"), ("b", "d")){
    line(a, b)
  }
  scircle((-4, -dy), (2, 1.5), "brancha")
  for (l, n) in ((lc, "c"), (ld, "d"), (le, "e")){
    let loc = (l.at(0)-5, l.at(1)-s/2-dy)
    circle(loc, radius:0.4, name: n, stroke: if n == "c" {red} else {black})
    content(loc, text(14pt)[$#n$])
  }
  for (a, b) in (("c", "d"), ("d", "e"), ("c", "d")){
    line(a, b)
  }
  scircle((4, -dy), (1, 1), "branchb")
  circle((4, -dy), radius:0.4, name: "e", stroke: red)
  content((4, -dy), text(14pt)[$e$])
  scircle((-6, -2*dy), (1, 1), "branchaa")
  circle((-6, -2*dy), radius:0.4, name: "e", stroke: red)
  content((-6, -2*dy), text(14pt)[$e$])
  scircle((-2, -2*dy), (0.5, 0.5), "branchab")
  scircle((4, -2*dy), (0.5, 0.5), "branchba")
  scircle((-6, -3*dy), (0.5, 0.5), "branchaaa")
  line("branch", "brancha")
  line("branch", "branchb")
  line("brancha", "branchaa")
  line("brancha", "branchab")
  line("branchb", "branchba")
  line("branchaa", "branchaaa")
  content((-5, -dy/2+0.5), text(12pt)[$G \\ N[a]$])
  content((3.5, -dy/2), text(12pt)[$G \\ N[b]$])
  content((-6.8, -3*dy/2), text(12pt)[$G \\ N[c]$])
  content((-1.5, -3*dy/2-0.4), text(12pt)[$G \\ N[d]$])
  content((-4.8, -5*dy/2-0.4), text(12pt)[$G \\ N[e]$])
  content((5.2, -3*dy/2-0.4), text(12pt)[$G \\ N[e]$])
}),
caption: [The simplest branching algorithm.],
)

Complexity of a branching algorithm is always described as $O(gamma^n)$ where $gamma$ is the branching factor and $n$ is the size of the problem.

#pagebreak()

=== The branching Algorithms for MIS

Various branching algorithms are proposed to improve the complexity, with human-designed rules.

#align(center,
table(
  columns: (auto, auto, auto, auto),
  table.header(hd[Year], hd[Running times], hd[References], hd[Notes]),
  s[1977], s[$O^*(1.2600^n)$], s[Tarjan & Trojanowski, 1977], s[],
  s[1986], s[$O^*(1.2346^n)$], s[Jian, 1986], s[],
  s[1986], s[$O^*(1.2109^n)$], s[Robson, 1986], s[],
  s[1999], s[$O^*(1.0823^m)$], s[Beigel, 1999], s[],
  s[2001], s[$O^*(1.1893^n)$], s[Robson, 2001], s[],
  s[2003], s[$O^*(1.1254^n)$ for 3-MIS], s[Chen, 2003], s[],
  s[2005], s[$O^*(1.1034^n)$ for 3-MIS], s[Xiao, 2005], s[],
  s[2006], s[$O^*(1.2210^n)$], s[Fomin, 2006], s[Measure and conquer, mirror rule],
  s[2006], s[$O^*(1.1225^n)$ for 3-MIS], s[Fomin and Hoie, 2006], s[],
  s[2006], s[$O^*(1.1120^n)$ for 3-MIS], s[Furer, 2006], s[],
  s[2006], s[$O^*(1.1034^n)$ for 3-MIS], s[Razgon, 2006], s[],
  s[2008], s[$O^*(1.0977^n)$ for 3-MIS], s[Bourgeois, 2008], s[],
  s[2009], s[$O^*(1.0919^n)$ for 3-MIS], s[Xiao, 2009], s[],
  s[2009], s[$O^*(1.2132^n)$], s[Kneis, 2009], s[Satellite rule],
  s[2013], s[$O^*(1.0836^n)$ for 3-MIS], s[Xiao, 2013], s[#highlight[SOTA]],
  s[2016], s[$O^*(1.2210^n)$], s[Akiba, 2016], s[#highlight[PACE winner]],
  s[2017], s[$O^*(1.1996^n)$], s[Xiao, 2017], s[SOTA],
)
)

Can we combine the branching algorithm with the tensor network method?

#pagebreak()

=== The optimal branching algorithm

We use the tensor network to extract the local information#footnote(text(12pt)[*X.-Z. Gao*., Y.-J. Wang, P. Zhang, J.-G. Liu, http://arxiv.org/abs/2412.07685, (2024), also see #link("https://github.com/ArrogantGao/OptimalBranching.jl")],), and then automatically search the optimal branching rules.

#align(center, image("figs/ob_rules.png", width: 600pt))


#align(center,
grid(
  columns: 3,
  text(20pt, black)[Naive branching \ 4 branches, each fix 5 variables \ $gamma^n = 4 times gamma^(n-5)$ \ $gamma approx 1.3195$],
  h(60pt),
  text(20pt, black)[Optimal branching \ 3 branches, fix [4, 5, 5] variables \ $gamma^n = gamma^(n-4) + 2 times gamma^(n-5)$ \ $gamma approx 1.2671$],
)
)

#pagebreak()

=== The optimal branching algorithm

The process of finding the optimal branching rules is as the following:

#align(center, canvas({
  import draw: *
  content((0, 0), box([Possible assignments ${bold(s)_1, bold(s)_2, dots, bold(s)_l}$], stroke: black, inset: 10pt), name: "oracle")
  content((0, -3), box([Candidate clauses $cal(C) = {c_1, c_2, dots, c_m}$], stroke: black, inset: 10pt), name: "clauses")
  content((0, -6), box([Optimal branching rule $cal(D) = c_(k_1) or c_(k_2) or dots$], stroke: black, inset: 10pt), name: "branching")
  line("oracle", "clauses", mark: (end: "straight"))
  line("clauses", "branching", mark: (end: "straight"))
}))

The candidate clauses are the combinations of the possible assignments, for example:

$
  "combine"(not a and b and not c and d and not e, a and b and c and not d and not e) = b and not e
$

Since we say it is a combination of the 3rd and 4th assignments, we say it covers ${3, 4}$.
The clauses are generated iteratively.

#pagebreak()

=== The optimal branching algorithm

Then we solve a set covering problem by formulating it as a mixed integer programming problem#footnote(text(12pt)[T. Achterberg, Math. Program. Comput., 1 (2009), pp. 1–41],).

$
min_(gamma, bold(x)) gamma " s.t. " & sum_(i=1)^m gamma^(-Delta rho(c_i)) x_i = 1,\
& union.big_(i = 1, dots, m\ x_i = 1) J_i = {1, 2, dots, n},  #h(50pt) arrow "valid branching rule"\
& x_i in {0, 1} #h(161pt) arrow "a clause is selected or not"
$
where $rho(c_i)$ is the number of variables fixed in the clause $c_i$.

The chosen clauses form the optimal branching rule:
$
  cal(D) = c_(k_1) or c_(k_2) or dots or c_(k_m)
$
with the minimum $gamma$ among all the possible branching rules.

#pagebreak()

=== Example

A bottle neck case has been reported in Xiao's work#footnote(text(12pt)[M. Xiao and H. Nagamochi, Theor. Comput. Sci., 469 (2013), pp. 92–104],), with a branching factor of $1.0836$.

#grid(columns: 3,
  image("figs/ob_bottleneck.png", width: 300pt),
  h(30pt),
  align(horizon, text(20pt, black)[
    - 71 possible assignments, 15782 candidate clauses. However, the optimal branching rule can be solved in few seconds.
    - Optimal branching vector: $[10, 16, 26, 26]$, which $gamma = 1.0817$
  ]),
)

#pagebreak()

=== Numerical results

// The optimal branching algorithm is combined with the traditional branching algorithm to solve the MIS problem.

#align(center,
  image("figs/ob_benchmark_table.png", width: 350pt),
)

#subpar.grid(
  figure(image("figs/ob_benchmark_average.png"), caption: [
    Average case.
  ]), <a>,
  figure(image("figs/ob_benchmark_worst.png"), caption: [
    Worst case.
  ]), <b>,
  columns: (250pt, 250pt),
  label: <full>,
  caption: [Number of branches generated by the branching algorithms on 1000 random graphs of different sizes.],
)

#pagebreak()
=== Summary

A new method to automatically discover the optimal branching rules for the branching algorithm is proposed, by combining the tensor network method and the branching algorithm.
With this method, we achieved an average branching factor of $O(1.0441)$ on random three-regular graphs, which outperforms the SOTA.

Advantages:
- generate the branching rules automatically without human effort
- better scaling than the traditional branching algorithm

Disadvantages:
- vertices selection strategy is needed, the branching rule is not globally optimal
- solving the rule can be computationally expensive

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
// Complexity of a direct sum of the Coulomb interaction in doubly periodic systems is about $O(N^2 epsilon^(-1/3))$.

#pagebreak()

=== Algorithms for Q2D charged systems

Methods have been developed to accelerate the Coulomb interaction in Q2D systems.

The very first method is the Ewald2D#footnote(text(12pt)[David E. Parry, Surf. Sci. 49 (1975), no. 2, 433–440.],) method based on the Ewald splitting of the Coulomb kernel. It is accurate but with $O(N^2)$ complexity.

To reduce the complexity, most methods rely on the following two strategies:

- Fourier spectral methods#footnote(text(12pt)[Ondrej Maxian, Rau´l P. Pel´aez, Leslie Greengard, and Aleksandar Donev, J. Chem. Phys. 154 (2021), no. 20, 204107.],) #footnote(text(12pt)[Franziska Nestler, Michael Pippig, and Daniel Potts, J. Comput. Phys. 285 (2015), 280–315.],) #footnote(text(12pt)[Davoud Saffar Shamshirgar, Joar Bagge, and Anna-Karin Tornberg, J. Chem. Phys. 154 (2021), no. 16, 164109.],): based on Ewald splitting and fast Fourier transform (FFT), with $O(N log N)$ complexity.

- Fast multipole methods#footnote(text(12pt)[Jiuyang Liang, Jiaxing Yuan, Erik Luijten, and Zhenli Xu, J. Chem. Phys. 152 (2020), no. 13, 134109.],) #footnote(text(12pt)[Ruqi Pei, Travis Askham, Leslie Greengard, and Shidong Jiang, J. Comp. Phys. 474 (2023), 111792.],) #footnote(text(12pt)[Wen Yan and Michael Shelley, J. Comput. Phys. 355 (2018), 214–232.],): based on the multipole expansion of the Coulomb kernel, with $O(N)$ complexity.

#pagebreak()

=== Algorithms for Q2D charged systems

For doubly periodic systems, one major challenge is the large prefactor in $O(N)$ or $O(N log N)$ compared to 3D-PBC solvers, especially when the system is strongly confined in the $z$ direction, i.e., $L_z << L_x, L_y$.

For the FFT based methods, the singularity introduced by the Laplacian in the Fourier integral along the free direction leads to huge additional zero-padding#footnote(text(12pt)[Ondrej Maxian, Rau´l P. Pel´aez, Leslie Greengard, and Aleksandar Donev, J. Chem. Phys. 154 (2021), no. 20, 204107.],).

For the FMM based methods, more near field contributions is included#footnote(text(12pt)[Wen Yan and Michael Shelley, J. Comput. Phys. 355 (2018), 214–232.],).

== Fast Spectral Sum-of-Gaussian Method

=== Sum-of-Gaussian approximation of Coulomb kernel

A sum-of-Gaussian (SoG) approximation of the Coulomb kernel is utilized in our method#footnote(text(12pt)[*X.-Z. Gao*, S. Jiang, J. Liang, Z. Xu, Q. Zhou, http://arxiv.org/abs/2412.04595, (2024). Also see #link("https://github.com/HPMolSim/FastSpecSoG.jl")],), where
$
  1 / r approx (2 log b) / sqrt(2 pi sigma^2) sum_(l = - infinity)^(infinity) 1 / (b^l) e^(- r^2 / (sqrt(2) b^l sigma)^2), " with" cal(E)_r < 2 sqrt(2) e^(- pi^2 / (2 log b))
$

Then we can split the potential into three parts:
$
  1 / r approx underbrace(( 1 / r - sum_(l = 0)^M w_l e^(- r^2 / s_l^2)) bb(1)_(r < r_c), "near-field")+ underbrace(sum_(l = 0)^m w_l e^(- r^2 / s_l^2), "mid-range") + underbrace(sum_(l = m + 1)^M w_l e^(- r^2 / s_l^2), "long-range")
$
where $s_l$ and $w_l$ are the nodes and weights of the SOG approximation#footnote(text(12pt)[Gregory Beylkin and Lucas Monz´on, Appl. Comput. Harmon. Anal. 28 (2010), no. 2, 131–149.],).
The nodes $s_l$ are arranged in monotone increasing order with $s_M >> L_x, L_y$, and $r_c$ is the cutoff radius to balance the near-field and far-field contributions.

#pagebreak()

=== Range splitting for far-field potential

How to determine the turning point $m$? 

We select $m$ so that $s_m < eta L_z < s_(m+1)$, where $eta$ is $O(1)$ constant.

\

// #align(center,
  #grid(columns: 3,
    image("figs/farfield.svg", width: 300pt),
    h(50pt),text(20pt)[
      The mid range part decays rapidly in $z$ direction, so that the mid range potential can be accurately periodicized in the z-direction and effectively handled using a pure Fourier spectral solver#footnote(text(12pt)[#link("https://github.com/HPMolSim/ChebParticleMesh.jl")],) with less zero padding.
      
      The long range part is extremely smooth in $[-L_z / 2, L_z / 2]$, and can be solved by a Fourier-Chebyshev method with $O(1)$ number of Chebyshev points.
    ]
  )


#pagebreak()
=== Numerical results

#figure(
  image("figs/sog_benchmark.png", width: 400pt),
  caption: [Total relative error is shown for (a) cubic systems and (b) strongly confined systems with fixed height $L_z = 0.3$, plotted against the number of particles $N$. Time costs are shown in (c) and (d), which shows a linear scaling with $N$.],
)

#pagebreak()

=== Summary

A fast and accurate solver for Q2D charged systems is developed based on the sum-of-Gaussian approximation of the Coulomb kernel.
The Coulomb kernel is splitted into three parts:
- near field terms: solved by real space truncation
- mid-range terms: solved by the Fourier spectral method with less zero padding
- long-range terms: solved by the Fourier-Chebyshev method with $O(1)$ number of Chebyshev points

It has the following advantages:
- spectrally accurate with rigorous error analysis
- not sensitive to the aspect ratio of the system
- smoothness and separability of the Gaussian removes the need of kernel truncation in the free direction
- does not require any upsampling in the gridding step
- all calculations are carried out in the fundamental cell itself

Currently, the major shortcoming of this method is its non-adaptive nature, and has a complexity of $O(N log N)$ rather than $O(N)$.

= Future Research Plans

== Fast Algorithms Based on the DMK Framework

The recently introduced dual-space multilevel kernel-splitting (DMK)#footnote(text(12pt)[S. Jiang, L. Greengard, http://arxiv.org/abs/2308.00292, (2023)],) method offers a powerful framework for solving a variety of long-range interactions and systems with periodic/partially periodic boundary conditions, which can be adaptive and with linear complexity.

#align(center, image("figs/dmk.png", width: 550pt))

== Tensor Network for Quantum Many-Body Systems

Tensor networks with specialized structures have long served as ansatz for quantum many-body systems, including matrix product states (MPS) and projected entangled-pair states (PEPS). 

With the advanced tensor network contraction techniques, I am interested in exploring more flexible tensor network structures as ansatz, which may yield more accurate representations of quantum many-body states.

Additionally, I aim to develop efficient algorithms for evolving these states by integrating them with the Dirac–Frenkel/McLachlan variational principle#footnote(text(12pt)[A. Raab. Chem. Phys. Lett., 319(5):674–678, 2000.],), the automatic differentiation and proper pre-conditioning#footnote(text(12pt)[M. Ganahl, J. Rincón, G. Vidal. Phys. Rev. Lett. 118, 220402 (2017).],).

== Sparse Tensor Networks Contraction via Optimal Branching

The optimal branching algorithm can be applied to contract the sparse tensor networks.
Assume that $T$ is a sparse tensor, the non-zeros values are listed as below:

#align(center,
grid(columns: 3,
  canvas(length: 35pt, {
    import draw: *
    let T = (0, 0)
    let A = (1, 0)
    let B = (0, 1)
    let C = (-1, 0)
  let D = (0, -1)
  let locs = (T, A, B, C, D)
  let mid = (A, B, C, D).map(v => (v.at(0) * 0.5, v.at(1) * 0.5))
  let e = ((0, 1), (0, 2), (0, 3), (0, 4))
  let c = ((0, $T_(i j k l)$), (1, $A_(i *)$), (2, $B_(j *)$), (3, $C_(k *)$), (4, $D_(l *)$))
  let s = 2
  let indices = ($i$, $j$, $k$, $l$)
  show-graph-black(locs.map(v => (v.at(0) * (s + 0.8), v.at(1) * (s + 0.8))), e, radius:0.0)
  show-graph-black(locs.map(v => (v.at(0) * s, v.at(1) * s)), e, radius:0.5)
  show-graph-content(locs.map(v => (v.at(0) * s, v.at(1) * s)), e, c, radius:0.5, fontsize: 16pt)
  for (i, v) in mid.map(v => (v.at(0) * (s), v.at(1) * (s))).enumerate() {
    circle(v, radius:0.3, fill:white, stroke:none)
    content(v, text(15pt, black)[#indices.at(i)])
  }
}),
h(100pt),
align(center + horizon,
canvas({
  import draw: *
  content((0,0), text(20pt, black)[
    #table(
    columns: (auto, auto, auto, auto, auto),
    inset: 10pt,
    align: horizon,
    table.header(
      [i], [j], [k], [l], [value],
      ),
      [1], [1], [0], [1], [0.1],
      [1], [1], [1], [0], [0.2],
      [1], [1], [0], [0], [0.3],
      [1], [0], [0], [0], [0.4],
      [0], [1], [0], [0], [0.5]
    )
  ])
})),
)
)

Such a contraction can be viewed as a branching problem, one can use the optimal branching algorithm to find the optimal way. 

Such sparsity is common in many problems, including probabilistic inference, combinatorial optimization, and quantum circuit simulations#footnote(text(12pt)[I.L. Markov, Y. Shi, SIAM J. Comput. 38, 963–981 (2008).],).

= Acknowledgements

#pagebreak()

#let figsize = 130pt

#align(center,
grid(columns: 5, 
  grid(rows : 3, image("photos/zechenggan.jpeg", width : figsize), "" ,  text[Prof. Zecheng Gan \ HKUST(GZ)]),
  h(50pt), 
  grid(rows : 3, image("photos/zhenlixu.jpg", width : figsize), "" ,  text[Prof. Zhenli Xu \ SJTU]),
  h(50pt), 
  grid(rows : 3, image("photos/shidongjiang.jpeg", width : figsize), "" ,  text[Prof. Shidong Jiang \ Flatiron Institute]),
  )
)

// #pagebreak()

// \

#align(center,
grid(columns: 3, 
  grid(rows : 3, image("photos/jinguoliu.jpeg", width : figsize), "" ,  text[Prof. Jinguo Liu, \ HKUST(GZ)]),
  h(60pt), 
  grid(rows : 3, image("photos/panzhang.jpeg", width : figsize), "" ,  text[Prof. Pan Zhang, \ ITP, CAS]),
  )
)

Also thank Jiuyang Liang (SJTU & FI), Qi Zhou (SJTU), Martin Roa-Villescas (TUE) and Yijia Wang (ITP, CAS) for collaboration.

= Appendix

== Details of Tensor Network Based Algorithms

=== Tropical Tensor Network

In tropical semiring, the multiplication and addition are defined as:
$
  a times.circle b &= a + b \ a plus.circle b &= max(a, b)
$

By replacing the matrix multiplication with the tropical matrix multiplication, we get the tropical TN, where the contraction results:
$
  T_(i_y ...) = max_(i in.not {i_y ...}) (A_(i_a ...) + B_(i_b ...) + C_(i_c ...) + ...)
$
which is the maximum of the sum of the elements among all the possible assignments.

Very useful in combinatorial optimization problems and ground state search.

#pagebreak()

=== Marginal probability by differential programming

A direct forward way: contract $N$ networks with $1$ open edge.

Our way: use the backward mode automatic differentiation, since:
$
  Z = sum_(i_0) Z_(i = i_0) times bb(1)_(i_0)
$
then
$
  frac(partial Z, partial bb(1)_(i_0)) = Z(i = i_0)
$
which gives the marginal probability together with the partition: $P_(i = i_0) = Z(i = i_0) / Z$.

With the backward mode automatic differentiation, the marginal probability of all variables can be computed in a single contraction and a backward pass.

#pagebreak()

=== Most probable explanation

The MPE task is to find a most likely assignment to all variables so that the probability is maximized:
$
  arg max_(i_0, j_0, ...) P(i = i_0, j = j_0, ...)
$
Solving the maximum $P$ is equivalent to contracting the logarithm of the network under a tropical semiring.
The contraction order is exactly the same as the one in MAR.

Utilizing a similar backpropagation method#footnote(text(12pt)[J. G. Liu, X. Gao, M. Cain, M. D. Lukin, and S. T. Wang, SIAM J. Sci. Comput. 45, A1239 (2023).],) as employed in MAR, one can derive the gradient with respect to the all-one tensor. 
In each gradient computation, there exists precisely one non-zero element, representing the most likely assignment.

#pagebreak()

=== Tensor Network for Maximum Independent Set Problem

A tropical TN can be used to solve the MIS problem, a simple example is shown below:

#align(center, canvas(length:40pt, {
  import draw: *
  // petersen graph
  let vertices1 = ((0, 0), (5, 0), (5, 5), (0, 5))
  let edges = ((0, 1), (1, 2), (2, 3), (3, 0))
  show-graph((vertices1).map(v=>(0.5 * (v.at(0)), 0.5 * (v.at(1)) + 1.24)), edges, radius:0.1)
  // content((2.5, -1.5), text(16pt)[(a)])

  content((4, 2.5), text(30pt)[$arrow$])

  set-origin((6, 1.25))
  show-graph((vertices1).map(v=>(0.5 * (v.at(0)), 0.5 * (v.at(1)))), edges, radius:0.0)

  let vertices3 = ((0, 0), (5, 0), (5, 5), (0, 5), (-1.5, -1.5), (6.5, -1.5), (6.5, 6.5), (-1.5, 6.5))
  let edges3 = ((0, 4), (1, 5), (2, 6), (3, 7))
  show-graph((vertices3).map(v=>(0.5 * (v.at(0)), 0.5 * (v.at(1)))), edges3, radius:0.0)
  // 
  for (i, v) in (vertices1.map(v=>(0.5 * (v.at(0)), 0.5 * (v.at(1))))).enumerate() {
    circle(v, radius:0.3, fill:white, stroke:none)
    content(v, text(20pt, black)[#(4 - i)])
  }

  let vertices2 = ((2.5, 0), (5, 2.5), (2.5, 5), (0, 2.5), (-1.5, -1.5), (6.5, -1.5), (6.5, 6.5), (-1.5, 6.5))
  let edges2 = ()
  show-graph-content((vertices2).map(v=>(0.5 * (v.at(0)), 0.5 * (v.at(1)))), edges2, ((0, "B"), (1, "B"), (2, "B"), (3, "B"), (4, "W"), (5, "W"), (6, "W"), (7, "W")), radius:0.4, fontsize: 20pt)

  // content((1.25, -2.75), text(16pt)[(b)])
}))

with 
$
  B = mat(
    bb(1), bb(1);
    bb(1), bb(0);
  ), " and " 
  W = mat(
    bb(0);
    bb(1);
  )
$
where $bb(1)$ and $bb(0)$ are identity element of tropical multiplication and addition, corresponding to $0$ (add) and $- infinity$ (max) of normal algebra, respectively.

#pagebreak()

=== Solving the set covering problem via Linear Integer Programming

In the original set covering problem, the function to be optimized is not linear. We solve it iteratively.

Fixed the branching complexity $gamma$, find a solution to $x$ that satisfies
$
min_(x) sum_(i=1)^(|cal(C)|) gamma^(-Delta rho(c_i)) x_i,  "s.t." union.big_(i = 1, dots, |cal(D)|,\ x_i = 1) J_i = {1, 2, dots, |cal(S)_R|}
$
It corresponds to the following WMSC problem:
$
cases(
"Alphabet:" {1, 2, dots, |cal(S)_R|},
"Sets:" {J_1, J_2, dots, J_(|cal(C)|)},
"Weights:" i arrow.r.bar gamma^(-Delta rho(c_i))
)
$
After each iteration, the branching complexity $gamma$ is updated, coverge in a few iterations.

#pagebreak()

=== The branching algorithm used in benchmark

#let s(name) = table.cell(align(center + horizon, text(16pt)[#name]))
#align(center,
table(
    columns: (150pt, auto, auto, auto),
    table.header(
    bdiagbox(s[Reduction], s[Branching], width:150pt), s[Optimal Branching \ (this work)], s[Xiao 2013], s[Akiba 2015],
    ),
    s[d1/d2 reduction], s[ob], s[-], s[akiba2015],
    s[d1/d2 reduction\ Xiao's rules], s[ob+xiao], s[xiao2013], s[-],
    s[d1/d2 reduction\ Xiao's rules\ packing rule], s[-], s[-], s[akiba2015+xiao&packing],
))

// === Ewald summation

// A traditional method to deal with the long range Coulomb interaction is Ewald summation, where the Coulomb kernel is split into short-range and long-range parts:
// $
//   1 / r = (#text[erfc($alpha r$)]) / r + (#text[erf($alpha r$)]) / r
// $
// where $#text[erfc($x$)]$ is the complementary error function and $#text[erf($x$)]$ is the error function.

// #align(center,
//   image("figs/erfc.png", width: 350pt),
// )


// #pagebreak()

// === Ewald summation

// The short-range part then is truncated in real space, and the long-range part is truncated in Fourier space.

// Due to the periodicity in $x$ and $y$, the long-range interaction energy can be given by a 
// $
//   E_("long") = 1/(L_x L_y) sum_(i, j) q_i q_j sum_(k_x, k_y) integral_(- infinity)^(infinity) e^(-i arrow(k) dot arrow(r)_(i j))/(k^2) e^(-k^2/(4 alpha^2)) d k_z - alpha / sqrt(pi) sum_i q_i^2
// $
// where the infinite integral can be calculated analytically:
// $
//   h / pi integral_(- infinity)^(infinity) e^(-i k_z z_(i j))/(k^2) e^(-k_z^2/(4 alpha^2)) d k_z = underbrace(e^(h abs(z_(i j)))"erfc"(h/(2 alpha) + alpha abs(z_(i j))), xi^+(z_(i j))) + underbrace(e^(-h abs(z_(i j)))"erfc"(h/(2 alpha) - alpha abs(z_(i j))), xi^-(z_(i j)))
// $
// where $h = sqrt(k_x^2 + k_y^2)$. This function decay exponentially fast in the reciprocal space.

// The resulting algorithm is called Ewald2D#footnote(text(12pt)[D. Parry, Surf. Sci. 49 (2) (1975) 433–440.],), with a complexity of $O(N^2 log epsilon)$.
// // , which prohibits the simulation of large systems.





// == Sum-of-Exponential Ewald2D Method

// === The SOE approximation

// Inspired by the fast Gaussian transform via sum-of-exponential (SOE) approximation#footnote(text(12pt)[S. Jiang, L. Greengard, Commun. Comput. Phys. 31 (1) (2022) 1–26],), we propose a sum-of-exponential Ewald2D method#footnote(text(12pt)[Z. Gan, X. Gao, J. Liang, and Z. Xu. arXiv preprint arXiv:2405.06333, 2024; also see #link("https://github.com/HPMolSim/SoEwald2D.jl")],) utilizing VPMR#footnote(text(12pt)[Z. Gao, J. Liang, Z. Xu, J. Sci. Comput. 93 (2) (2022) 40],) method.

// The Guassian kernel is approximated by a sum of exponential functions as follows:
// $
//   e^(- alpha^2 t^2) approx sum_(l = 1)^M  w_l e^(- s_l alpha abs(t))
// $
// where $w_l$ and $s_l$ are complex parameters of the SOE approximation, and the following error bound is satisfied:
// $
//   abs( e^(- alpha^2 t^2) - sum_(l = 1)^M  w_l e^(- s_l alpha abs(t)) ) < epsilon.
// $
// In VPMR, $M = 4$, $8$ and $16$ gives the error bound $epsilon = 10^(-4)$, $10^(-8)$ and $10^(-14)$, respectively.

// #pagebreak()

// === Approximation of the Ewald2D kernel

// In Ewald2D summation, we observed the following identity:
// $
//   xi^+ (z) = e^( h abs(z))"erfc"(h/(2 alpha) + alpha abs(z)) & = (2 alpha) / sqrt(pi) e^(-h^2/(4 alpha^2)) e^(h abs(z)) integral_(abs(z))^(infinity) e^(-alpha^2 t^2 - h t) d t\
// $
// and the SOE approximation is applied to get:
// $
//    xi^+ (z) approx (2 alpha) / sqrt(pi) e^(-h^2/(4 alpha^2)) e^(h abs(z)) integral_(abs(z))^(infinity) sum_(l = 1)^M  w_l e^(- s_l alpha abs(t)  - h t) d t = (2 alpha) / sqrt(pi) e^(-h^2/(4 alpha^2)) sum_(l = 1)^M  w_l e^(- s_l alpha abs(z)) / (alpha s_l + h)
// $
// For $xi^-(z)$, the same approximation is applied.

// Now the double summation over $i$ and $j$ in energy is approximated by
// $
//   sum_(i, j) q_i q_j xi (z_(i j)) approx sum_(l = 1)^M C_l sum_(i, j) q_i q_j e^(- s_l alpha abs(z_(i j)))
// $
// where $C_l$ is a constant irrelevant to $i$ and $j$.

// #pagebreak()

// === Summing up the exponentials via sorting

// The double summation can be simplified as 
// $
//   S = sum_(i, j) q_i q_j e^(- abs(z_i - z_j))
// $
// which still takes $O(N^2)$ operations due to the absolute value.

// For further acceleration, we reorder the indices via sorting (at most $O(N log N)$ operations):
// $
//   z_1 < z_2 < ... < z_N
// $
// so that the absolute value can be removed: 
// $
//   S = sum_(i = 1)^N q_i e^(-z_i) underbrace(sum_(j = 1)^i q_j e^(z_j), A_i) = sum_(i = 1)^N q_i e^(-z_i) A_i
// $
// where the array $A$ can be computed iteratively in $O(N)$ operations, then calculating $S$ takes another $O(N)$ operations.

// #pagebreak()

// === Complexity of the SOEwald2D method

// Utilizing the SOE approximation, the double summation can be calculated in $O(N)$ operations.
// Calculation of the energy can be simplified as
// $
//   sum_(k_x, k_y) e^(-h^2 / (4 alpha^2)) sum_(i, j) S(k_x, k_y, r_i, r_j) = underbrace(sum_(k_x, k_y), O(N^0.4)) e^(-h^2 / (4 alpha^2)) underbrace(f(k_x, k_y), O(N))
// $
// However, a summation in the reciprocal space is still required, and the resulting total complexity is $O(N^(1.4))$ since number of the Fourier mode to be summed is of $O(N^(0.4))$.

// #pagebreak()

// === Random batch sampling

// To further reduce the complexity, we use the the random batch sampling technique#footnote(text(12pt)[S. Jin, L. Li, Z. Xu, Y. Zhao, SIAM J. Sci. Comput. 43 (4) (2021) B937–B960.],), where a importance sampling is applied:
// $
//   sum_(k_x, k_y) e^(-h^2 / (4 alpha^2)) f(k_x, k_y) approx H / P sum_((k_x, k_y) in cal(K)_P) f(k_x, k_y)
// $
// where $P$ Fourier modes are selected using the Gaussian as distribution, $H$ is the summation of Guassians.

// It has been shown that for a system with fixed density, $P ~ O(1)$ is sufficient to achieve the correct equilibrium state, and gives an accurate ensemble average.

// The resulting algorithm is called RBSE2D, with a complexity of $O(N)$.


// #pagebreak()

// === Numerical Results

// The accuracy of the SOEwald2D method is determined by the error bound of the SOE approximation, and coverage exponentially fast as the Ewald2D method.

// #figure(
//   image("figs/soewald2d_error.png", width: 600pt),
//   caption: [Error of the SOEwald2D method (a) with different number of term $M$ used in the SOE approximation (b) with different system sizes ($L_z$ is fixed).],
// )

// #pagebreak()

// === Numerical Results

// We use the RBSE2D method to simulate a typical Q2D electrolyte system as shown below:

// #figure(
//   image("figs/soewald2d_md.png", width: 500pt),
//   caption: [The charge density of the Q2D electrolyte system in $z$ direction with different $P$, inset shows the relative error of the average energy.],
// )

// #pagebreak()

// === Numerical Results

// Time complexity of the method is verified as below:

// #figure(
//   image("figs/soewald2d_complexity.png", width: 500pt),
//   caption: [Complexity of the Ewald2D method, SOEwald2D method and RBSE2D method.],
// )

// #pagebreak()

// === Summary

// The SoEwald2D utilized the SOE approximation and spirit of the FGT to accelerate the Ewald2D method.

// It has the following advantages:
// - Reducing the complexity from $O(N^2)$ to $O(N^1.4)$ with a controlled error bound.
// - With the random batch sampling technique, the complexity can be further reduced to $O(N)$.
// - Not sensitive to the aspect ratio of the system.

// However, it also has some disadvantages:
// - A sorting step is required, which is not efficient for parallel computing.
// - Importance sampling is necessary to achieve linear complexity, which may not be efficient for cases requiring exact results.

== Details of Fast Summation Algorithms

=== Mid range part

The mid range energy is given by:
$
  U_("mid") = sum_(i = 1)^N q_i cal(F)^(-1) [hat(W) e^(i k r_i) sum_(l = 0)^m w_l s_l^3 e^(- (k^2 s_l^2) / 4) abs(hat(W))^(-2) sum_(j = 1)^N q_j hat(W) e^( - i k r_j)]
$

A simple example of the Fourier spectral method:

#align(center,
grid(columns: 9,
  align(center + horizon, image("figs/particlemesh_1.svg", width: 100pt)),
  align(center + horizon, text(30pt)[$arrow$]),
  align(center + horizon, image("figs/particlemesh_2.svg", width: 100pt)),
  align(center + horizon, text(30pt)[$arrow$]),
  align(center + horizon, image("figs/particlemesh_3.svg", width: 100pt)),
  align(center + horizon, text(30pt)[$arrow$]),
  align(center + horizon, image("figs/particlemesh_4.svg", width: 100pt)),
  align(center + horizon, text(30pt)[$arrow$]),
  align(center + horizon, image("figs/particlemesh_5.svg", width: 100pt)),
)
)

#pagebreak()

=== Fourier spectral solver for mid-range terms

After manually adding a vacuum layer in the $z$ direction, the mid-range potential can be solved by a standard Fourier spectral solver#footnote(text(12pt)[#link("https://github.com/HPMolSim/ChebParticleMesh.jl")],) in the following procedure:
- interpolate the particles using the window function to grid points
- perform 3D FFT
- scale the Fourier coefficients with the Green's function
- perform inverse FFT
- integrate from the 3D grid back to particles

The process is similar to that of the NUFFT#footnote(text(12pt)[Franziska Nestler, Michael Pippig, and Daniel Potts, J. Comput. Phys. 285 (2015), 280–315.],).
However, our method does not require any upsampling in the gridding step since the Gaussian decays quickly and it compensates the loss of accuracy in calculating the Fourier transform of the data.
// This is because we have taken advantage of the fact that the Fourier transform of the Gaussian decays quickly and it compensates the loss of accuracy in calculating the Fourier transform of the data.

#pagebreak()

=== Fourier-Chebyshev solver for long-range terms

For the long-range part, we map the charge density onto the Chebyshev proxy points in $z$, 
$
  phi_l (r) & = sum_(j) sum_(m) w_l q_j e^(- (r - r_j + m L)^2 / s_l^2) \
  & = pi / (L_x L_y) sum_(k_x, k_y) w_l s_l^2 e^(- (h^2 s_l^2) / 4) e^(i (k_x x + k_y y)) sum_(j) q_j e^(- i (k_x x_j + k_y y_j)) e^(- (z - z_j)^2 / s_l^2)
$
where $m$ is the periodic index, and $h = (k_x, k_y)$ is the wave vector in the $x y$ plane.

#pagebreak()

=== Fourier-Chebyshev solver for long-range terms

// Use Chebyshev function for fast summation in $z$ direction:
// $
//   sum_(i, j) f(x_i - x_j) = sum_(i) g(x_i) 
// $
// where 
// $
//   g(x) = sum_j f(x - x_j) approx 
// $

Consider the following summation in 1D:
$
  S = sum_(i, j) f(x_i - x_j)
$
where $x in (-1, 1)$. If we now the following function:
$
  g(x) = sum_j f(x - x_j)
$
then we can directly evaluate all $g(x_i)$ and sum them up.
Chebyshev interpolation is a good choice for this purpose, especially when $f$ is extremely smooth.

In our case, one first evaluate the values of the long-range Gaussians at the Chebyshev points, and then use the resulting Chebyshev series to evaluate the potential on the particle positions.

#pagebreak()

=== Fourier-Chebyshev solver for long-range terms

This changes the problem into a set of 2D problems.

#align(center,
  image("figs/sog_long.png", width: 500pt),
)

On each layer the problem can be solved by a standard 2D Fourier spectral solver.
Then the potential is integrated back to particles using the inverse Chebyshev transform.