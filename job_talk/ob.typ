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
// #let problem = thmplain("problem", "Problem").with(numbering: none)

#let hd(name) = table.cell(text(12pt)[#name], fill: green.lighten(50%))
#let s(name) = table.cell(text(12pt)[#name])

#set page(height: auto)
#set par(justify: true)

#let globalvars = state("t", 0)
#let timecounter(minutes) = [
  #globalvars.update(t => t + minutes)
  #place(dx: 95%, dy: -10%, align(right, text(16pt, red)[#context globalvars.get()min]))
]
#let clip(image, top: 0pt, bottom: 0pt, left: 0pt, right: 0pt) = {
  box(clip: true, image, inset: (top: -top, right: -right, left: -left, bottom: -bottom))
}

#import "@preview/touying-flow:1.1.0":*
#show: flow-theme.with(
  aspect-ratio: "4-3",
  footer: self => self.info.title,
  // footer-alt: self => self.info.subtitle,
  navigation: "none",
  primary:rgb("#004098"),//rgb(0,108,57),//rgb("#006c39"),
  secondary:rgb("#004098"),//rgb(161,63,61),//rgb("#a13f3d"),
  text-font: ("Libertinus Serif"),
  text-size: 20pt,

  config-info(
    title: [Automated Discovery of the Optimal Branching Rules],
    subtitle: [arXiv:2412.07685],
    author: text(23pt)[Xuanzhao Gao],
    institution: text(20pt)[Hong Kong University of Science and Technology],
    date: text(23pt)[2025-1-10],
  ),
  // config-common(show-notes-on-second-screen: right),
)

#show link: underline

#title-slide()

= Discovering the Optimal Branching Rules

== Background

=== The maximum independent set problem
#timecounter(1)

#align(center, box([One of the first batch of 21 NP-hard problems proved by @Karp1972.], stroke: black, inset: 10pt))

An independent set is a set of vertices in a graph, no two of which are adjacent.

#align(center, canvas({
  import draw: *
  for (loc, name, color) in (((0, 0), "a", black), ((0, 3), "b", red), ((3, 0), "d", red), ((3, 3), "c", black), ((5, 1.5), "e", red)) {
    circle(loc, radius:0.5, name: name, stroke: color)
    content(loc, [$#name$])
  }
  for (a, b) in (
    ("a", "b"),
    ("a", "d"),
    ("a", "c"),
    ("c", "d"),
    ("b", "c"),
    ("c", "e"),
    ) {
    line(a, b)
  }
  content((1.5, -1.5), [$G = (V, E)$, MIS = ${b, d, e}$, size $alpha(G) = 3$])
}))

MIS problem has an exponential large solution space and no polynomial-time algorithm is known to solve it exactly.

The branching algorithms @Fomin2013 and the tensor network methods are two popular methods to solve this problem.

#pagebreak()

=== Branching algorithm

Branching algorithm search the solution space in a tree-like structure. Complexity of a branching algorithm is always described as $O(gamma^n)$ where $gamma$ is the branching factor and $n$ is the size of the problem.

#grid(columns: 3,
align(center + horizon,
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
})
),
h(20pt),
align(center,
table(
  columns: (auto, auto, auto, auto),
  table.header(hd[Year], hd[Running times], hd[References], hd[Notes]),
  // s[1977], s[$O^*(1.2600^n)$], s[@Tarjan1977], s[],
  s[1986], s[$O^*(1.2346^n)$], s[@Jian1986], s[],
  s[1986], s[$O^*(1.2109^n)$], s[@Robson1986], s[],
  s[1999], s[$O^*(1.0823^m)$], s[@Beigel1999], s[],
  s[2001], s[$O^*(1.1893^n)$], s[@Robson2001], s[],
  s[2003], s[$O^*(1.1254^n)$ for 3-MIS], s[@Chen2003], s[],
  s[2005], s[$O^*(1.1034^n)$ for 3-MIS], s[@Xiao2005], s[],
  s[2006], s[$O^*(1.2210^n)$], s[@Fomin2006], s[],
  s[2006], s[$O^*(1.1225^n)$ for 3-MIS], s[@Fomin2006b], s[],
  s[2006], s[$O^*(1.1120^n)$ for 3-MIS], s[@Furer2006], s[],
  s[2006], s[$O^*(1.1034^n)$ for 3-MIS], s[@Razgon2006], s[],
  s[2008], s[$O^*(1.0977^n)$ for 3-MIS], s[@Bourgeois2008], s[],
  s[2009], s[$O^*(1.0919^n)$ for 3-MIS], s[@Xiao2009], s[],
  s[2009], s[$O^*(1.2132^n)$], s[@Kneis2009], s[Satellite rule],
  s[2013], s[$O^*(1.0836^n)$ for 3-MIS], s[@Xiao2013], s[#highlight[SOTA]],
  s[2016], s[$O^*(1.2210^n)$], s[@Akiba2016], s[#highlight[PACE winner]],
  s[2017], s[$O^*(1.1996^n)$], s[@Xiao2017], s[SOTA],
),
)
)

They rely on predesigned rules and search for special structures on the graph.

#pagebreak()

=== Tensor networks
#timecounter(1)

Tensor network (TN) is a powerful tool to represent and manipulate high-dimensional data, and have been used in various fields such as quantum physics @nielsen2010quantum and machine learning @stoudenmire2016supervised.
// Tensor networks have proven to be an powerful tool for solving large-scale problems.
// In 2019, Google announced quantum supremacy with their Sycamore chip @arute2019quantum; however, by 2022, their results were classically simulated using tensor network-based methods @pan2022simulation.

#align(center,
grid(columns: 3,
  figure(image("figs/big_batch.png", height: 280pt), caption: [
    #text(15pt)[Classical simulation of Google's Sycamore chip via tensor network @pan2022simulation.]
  ]),
  h(20pt),
  figure(image("figs/tn_inference.png", height: 280pt), caption: [
    #text(15pt)[Exact solve large-scale PACE inference problem via tensor network @Roa2024.]
  ]),
)
)

#pagebreak()

=== Tensor networks for the MIS problem
#timecounter(1)

The MIS problem can be efficiently solved and analyzed by tensor networks @ebadi2022quantum  @liu2023computing.

#align(center,
grid(columns: 5,
  align(center + horizon, image("figs/kingsubgraph.png", width:150pt)),
  h(30pt),
  align(center + horizon, text(40pt)[$arrow$]),
  h(30pt),
  align(center + horizon, image("figs/tn_mis_solutionspace.png", width:400pt)),
)
)

However, the tensor network method does not work well for non-geometric graphs.
Its complexity on 3-regular graphs is about $O(1.1224^n)$, far from the SOTA ($O^*(1.0836^n)$).

#pagebreak()

=== What is the difference?

#timecounter(1)

#align(center,
  canvas({
    import draw: *
    content((-11, -1), image("figs/ob_reduce_subgraph.png", width:150pt))

    
    content((8, 1), text(20pt, black)[Branching algorithm])
    content((8, -0.5), box(width: 280pt, text(20pt)[1. Search for structures in the sub-graph.]))
    content((8, -3.5), box(width: 280pt, text(20pt)[2. Find $N[v] subset N[w]$, satisfying the domination rule, fix $w = 0$ $ gamma = 1.0 $]))

    content((-3, 1), text(20pt, black)[Tensor network approach])
    content((-3, -1), box(width: 280pt, text(20pt)[1. Contract the local tensor network for the sub-graph and pick the non-zero elements.]))
    content((-3, -4), image("figs/ob_reduce_table.png", width:220pt))
    content((-3, -8), box(width: 280pt, text(20pt)[2. For four possible boundaries, fix all the variables and continue the contraction$ gamma^n = 4 times gamma^(n-5) arrow gamma approx 1.3195 $]))
    
  })
)

#align(center, box([Key point: No need to use all results, find the correct pattern!], stroke: black, inset: 10pt))

// #pagebreak()
== The optimal branching algorithm

===

#timecounter(4)

We use the tensor network to extract the local information, and then automatically search the optimal branching rules #footnote(text(12pt)[#link("https://github.com/ArrogantGao/OptimalBranching.jl")],) @Gao2024.

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

===

#timecounter(2)

The process of finding the optimal branching rules is as the following:

#align(center, canvas({
  import draw: *
  content((0, 0), box([Possible assignments ${bold(s)_1, bold(s)_2, dots, bold(s)_l}$], stroke: black, inset: 10pt), name: "oracle")
  content((0, -3), box([Candidate clauses $cal(C) = {c_1, c_2, dots, c_m}$], stroke: black, inset: 10pt), name: "clauses")
  content((0, -6), box([Optimal branching rule $cal(D) = c_(k_1) or c_(k_2) or dots$], stroke: black, inset: 10pt), name: "branching")
  line("oracle", "clauses", mark: (end: "straight"))
  line("clauses", "branching", mark: (end: "straight"))
}))

The candidate clauses are the combinations of the possible assignments, for example, a combination of the 3rd and 4th assignments is given by:

$
  "combine"(not a and b and not c and d and not e, a and b and c and not d and not e) = b and not e
$

we say $b and not e$ covers ${3, 4}$.
The clauses are generated iteratively.

#pagebreak()

===

#timecounter(2)

Then we solve a set covering problem by formulating it as a mixed integer programming problem.

$
min_(gamma, bold(x)) gamma " s.t. " & sum_(i=1)^m gamma^(-Delta rho(c_i)) x_i = 1,\
& union.big_(i = 1, dots, m\ x_i = 1) J_i = {1, 2, dots, n},  #h(50pt) arrow "valid branching rule"\
& x_i in {0, 1} #h(161pt) arrow "a clause is selected or not"
$
where $rho(c_i)$ is the size reduced by the clause $c_i$ of the problem.

The chosen clauses form the optimal branching rule:
$
  cal(D) = c_(k_1) or c_(k_2) or dots or c_(k_m)
$
with the minimum $gamma$ among all the possible branching rules.

== Numerical results

===

#timecounter(3)

A bottle neck case has been reported in Xiao's work @Xiao2013, with a branching factor of $1.0836$.

#grid(columns: 3,
  image("figs/ob_bottleneck.png", width: 300pt),
  h(30pt),
  align(horizon, text(20pt, black)[
    - 71 possible assignments, 15782 candidate clauses. 
    - The optimal branching rule can be solved in few seconds.
    - Size reduced by branches: $[10, 16, 26, 26]$, with $gamma = 1.0817$
  ]),
)

#pagebreak()
// The optimal branching algorithm is combined with the traditional branching algorithm to solve the MIS problem.

// #align(center,
//   image("figs/ob_benchmark_table.png", width: 450pt),
// )

// #subpar.grid(
//   figure(image("figs/ob_benchmark_average.png"), caption: [
//     Average case.
//   ]), <a>,
//   figure(image("figs/ob_benchmark_worst.png"), caption: [
//     Worst case.
//   ]), <b>,
//   columns: (280pt, 280pt),
//   label: <full>,
//   caption: [#text(15pt)[Number of branches generated by the branching algorithms on 1000 random graphs of different sizes.]],
// )

#grid(columns: 3,
  text(20pt)[
    The branch-and-reduce algorithm is improved by including our optimal branching algorithm to generate the branching rules on-the-fly.
    The resulting method is compared with the existing methods on different type of graphs.

    The results show that the optimal branching algorithm can significantly reduce the number of branches and reaches a better scaling.

    On 3-regular graphs, a complexity of $O(1.0441^n)$ is achieved, which outperforms the SOTA (xiao2013) with a factor of $1.0487$.
  ],
  h(30pt),
  grid(rows: 3,
    image("figs/ob_benchmark_table.png", width: 320pt),
    v(10pt),
    figure(image("figs/ob_benchmark_average.png", width: 340pt), caption : [
      #text(15pt)[Average number of branches generated by different branching algorithms on 1000 random graphs.]
    ])
  )
)

== Potential applications

=== Sparse Tensor Networks Contraction via Optimal Branching

#timecounter(2)

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

Such sparsity is common in many problems, including probabilistic inference, combinatorial optimization, and quantum circuit simulations @Markov2008.

== Summary

===

#timecounter(1)

A new method to automatically discover the optimal branching rules is proposed, by combining the tensor network method and the branching algorithm.

With this method, we achieved an average complexity of $O(1.0441^n)$ on random 3-regular graphs, which outperforms the SOTA.

Advantages:
- generate the branching rules automatically without human effort
- fully utlize the information of the sub-graph, which can be selected flexibly
// - no need to pick rules
- can be applied to different problems, not only the MIS problem

Disadvantages:
- solving the rule can be computationally expensive

= Other works

== Probabilistic inference in the era of tensor networks and differential programming

#grid(columns: 3,
  text(20pt)[
    The modern tensor network techniques are applied to solve the probabilistic inference tasks#footnote(text(12pt)[#link("https://github.com/TensorBFS/TensorInference.jl")],):
    - Contraction order optimization @gray2021hyper
    - Automatic differentiation @liao2019differentiable
    - Tropical tensor network @liu2023computing

    An exponential speedup is achieved for the typical inference tasks against the widely used solvers on the UAI problem set#footnote(text(12pt)[https://auai.org/uai2014/competition.shtml],).
  ],
  h(20pt),
  figure(
    image("figs/ti_benchmark.png", height: 330pt),
    caption: [#text(15pt)[Speedup of the inference tasks against the previous solvers.]],
  )
)

#pagebreak()

A custom GPU kernel#footnote(text(12pt)[Code: #link("https://github.com/TensorBFS/CuTropicalGEMM.jl"). Blog: #link("https://arrogantgao.github.io/blogs/CuTropicalGEMM/")],) for tropical matrix multiplication is implemented to accelerate the inference. It can achieve a 4000x speedup for tropical mat-mul compared to CPU.

#grid(columns: 3,
  figure(
    image("figs/cutropicalgemm.png", height: 280pt),
    caption: [#text(15pt)[Benchmark of the custom GPU kernel for tropical matrix multiplication.]],
  ),
  h(20pt),
  figure(
    image("figs/ti_gpu.png", height: 280pt),
    caption: [#text(15pt)[Speedup of the MMAP inference tasks by GPU againist the CPU version.]],
  )
)


== A fast spectral sum-of-Gaussians method for electrostatic summation in doubly periodic 3D systems

A sum-of-Gaussian (SOG) approximation @beylkin2010approximation of the Coulomb kernel is utilized in our method @gao2024fast, where
$
  1 / r approx (2 log b) / sqrt(2 pi sigma^2) sum_(l = - infinity)^(infinity) 1 / (b^l) e^(- r^2 / (sqrt(2) b^l sigma)^2), " with" cal(E)_r < 2 sqrt(2) e^(- pi^2 / (2 log b))
$

Then we can split the potential into three parts:
$
  1 / r approx underbrace(( 1 / r - sum_(l = 0)^M w_l e^(- r^2 / s_l^2)) bb(1)_(r < r_c), "near-field")+ underbrace(sum_(l = 0)^m w_l e^(- r^2 / s_l^2), "mid-range") + underbrace(sum_(l = m + 1)^M w_l e^(- r^2 / s_l^2), "long-range")
$
where $s_l$ and $w_l$ are the nodes and weights of the SOG approximation.
// The nodes $s_l$ are arranged in monotone increasing order with $s_M >> L_x, L_y$, and $r_c$ is the cutoff radius to balance the near-field and far-field contributions.

#pagebreak()

#grid(columns: 3,
  text(20pt)[
    Different strategies are used to handle the different parts of the potential:
    - The near field is then truncated in the real space.
    - The mid-range is computed by a Fourier spectral solver#footnote(text(12pt)[#link("https://github.com/HPMolSim/ChebParticleMesh.jl")],) with very little zero padding and no upsampling.
    - The long-range part is computed by a Fourier-Chebyshev solver with $O(1)$ number of Chebyshev points in $z$.

    The resulting method#footnote(text(12pt)[#link("https://github.com/HPMolSim/FastSpecSoG.jl")],) is spectrally accurate with rigorous error analysis, has a complexity of $O(N log N)$ and is not sensitive to the aspect ratio of the system.
  ],
  h(20pt),
  figure(
    image("figs/sog_benchmark.png", width: 350pt),
    caption: [#text(15pt)[Error and time cost for the SOG method in the (a,c) cubic and (b,d) strongly confined systems.]],
  )
)

#pagebreak()

= Future Research Plans

#pagebreak()

=== Tensor Network Algorithms

// With the advanced tensor network techniques, I am interested in exploring more flexible tensor network structures as ansatz, which may yield more accurate representations of quantum many-body states.

// I am also interested in developing efficient algorithms for evolving these states by integrating them with the Dirac–Frenkel/McLachlan variational principle @raab2000dirac, the automatic differentiation and proper pre-conditioning @ganahl2017continuous.

- branching based sparse tensor network contraction
- more flexible quantum many-body ansatz
// - proper pre-conditioning @ganahl2017continuous


=== Fast Summation Algorithms

- efficient methods based on the DMK framework @jiang2023dmk

// I am keen on developing efficient summation methods for long-range interacting systems, such as those described by Poisson’s equation and the Helmholtz equation. For example, I am interested in extending the recently introduced dual-space multilevel kernel-splitting (DMK) @jiang2023dmk framework to a variety of long-range interactions and systems with periodic/partially periodic boundary conditions.

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

#align(center,
grid(columns: 3, 
  grid(rows : 3, image("photos/jinguoliu.jpeg", width : figsize), "" ,  text[Prof. Jinguo Liu, \ HKUST(GZ)]),
  h(60pt), 
  grid(rows : 3, image("photos/panzhang.jpeg", width : figsize), "" ,  text[Prof. Pan Zhang, \ ITP, CAS]),
  )
)

Also thank Jiuyang Liang (SJTU & FI), Qi Zhou (SJTU), Martin Roa-Villescas (TUE) and Yijia Wang (ITP, CAS) for collaboration.

= Appendix

#pagebreak()

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

=== Tensor Network for Maximum Independent Set Problem

A tropical TN can be used to solve the MIS problem, a simple example is shown below:

#align(center, 
  grid(columns: 5,
  align(center + horizon, align(center, canvas({
  import draw: *
  for (loc, name, color) in (((0, 0), "a", black), ((0, 3), "b", black), ((3, 0), "d", black), ((3, 3), "c", black), ((5, 1.5), "e", black)) {
    circle(loc, radius:0.5, name: name, stroke: color)
    content(loc, [$#name$])
  }
  for (a, b) in (
    ("a", "b"),
    ("a", "d"),
    ("a", "c"),
    ("c", "d"),
    ("b", "c"),
    ("c", "e"),
    ) {
    line(a, b)
  }
}))),
h(30pt),
align(center + horizon, text(50pt)[$arrow$]),
h(30pt),
image("figs/tn_mis.png", width: 300pt),
)
)

// #align(center, canvas(length:40pt, {
//   import draw: *
//   // petersen graph
//   let vertices1 = ((0, 0), (5, 0), (5, 5), (0, 5))
//   let edges = ((0, 1), (1, 2), (2, 3), (3, 0))
//   show-graph((vertices1).map(v=>(0.5 * (v.at(0)), 0.5 * (v.at(1)) + 1.24)), edges, radius:0.1)
//   // content((2.5, -1.5), text(16pt)[(a)])

//   content((4, 2.5), text(30pt)[$arrow$])

//   set-origin((6, 1.25))
//   show-graph((vertices1).map(v=>(0.5 * (v.at(0)), 0.5 * (v.at(1)))), edges, radius:0.0)

//   let vertices3 = ((0, 0), (5, 0), (5, 5), (0, 5), (-1.5, -1.5), (6.5, -1.5), (6.5, 6.5), (-1.5, 6.5))
//   let edges3 = ((0, 4), (1, 5), (2, 6), (3, 7))
//   show-graph((vertices3).map(v=>(0.5 * (v.at(0)), 0.5 * (v.at(1)))), edges3, radius:0.0)
//   // 
//   for (i, v) in (vertices1.map(v=>(0.5 * (v.at(0)), 0.5 * (v.at(1))))).enumerate() {
//     circle(v, radius:0.3, fill:white, stroke:none)
//     content(v, text(20pt, black)[#(4 - i)])
//   }

//   let vertices2 = ((2.5, 0), (5, 2.5), (2.5, 5), (0, 2.5), (-1.5, -1.5), (6.5, -1.5), (6.5, 6.5), (-1.5, 6.5))
//   let edges2 = ()
//   show-graph-content((vertices2).map(v=>(0.5 * (v.at(0)), 0.5 * (v.at(1)))), edges2, ((0, "B"), (1, "B"), (2, "B"), (3, "B"), (4, "W"), (5, "W"), (6, "W"), (7, "W")), radius:0.4, fontsize: 20pt)

//   // content((1.25, -2.75), text(16pt)[(b)])
// }))

with 
$
  B = mat(
    0, 0;
    0, -infinity;
  ), " and " 
  W = mat(
    0;
    1;
  )
$
under a tropical semiring.

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

#pagebreak()

#bibliography("ref.bib", style: "apa")