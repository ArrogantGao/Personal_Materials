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

#let redtext(name) = text(fill:red)[#name]
#let bluetext(name) = text(fill:blue)[#name]

#let hd(name) = table.cell(text(12pt)[#name], fill: green.lighten(50%))
#let s(name) = table.cell(text(12pt)[#name])

#let hdl(name) = table.cell(text(15pt)[#name], fill: green.lighten(50%))
#let sl(name) = table.cell(text(15pt)[#name])

#set page(height: auto)
#set par(justify: true)

#let globalvars = state("t", 0)
#let timecounter(minutes) = [
  #globalvars.update(t => t + minutes)
  #place(top + right,text(16pt, red)[#context globalvars.get()min])
]
#let clip(image, top: 0pt, bottom: 0pt, left: 0pt, right: 0pt) = {
  box(clip: true, image, inset: (top: -top, right: -right, left: -left, bottom: -bottom))
}

#let leftright(a, b) = align(center, grid(columns: 2, gutter: 30pt, align(left, box(width: 350pt)[#a]), align(left, box(width: 350pt)[#b])))

#import "@preview/touying:0.5.5": *
#import themes.metropolis: *
#show: metropolis-theme.with(
  aspect-ratio: "4-3",
  footer: self => self.info.title,
  navigation: "none",
  text-font: ("Libertinus Serif"),
  text-size: 20pt,
  align: horizon,

  config-info(
    title: [Fast Summation Algorithms and Tensor Networks Methods for  Scientific Problems],
    // subtitle: [],
    author: text(23pt)[Xuanzhao Gao],
    institution: text(20pt)[Hong Kong University of Science and Technology],
    date: text(23pt)[2025-1-10],
  ),

  config-colors(
    primary: rgb("#eb811b"),
    primary-light: rgb("#d6c6b7"),
    secondary: rgb("#343d73"),
    neutral-lightest: rgb(255, 255, 255),
    neutral-dark: rgb("#24233b"),
    neutral-darkest: rgb("#23373b"),
  )
)

#show link: underline

#title-slide()

= Outline <touying:hidden>

#outline(title: none, indent: 1em, depth: 1)

= A Fast Spectral Sum-of-Gaussians Method for Coulomb Interaction

#text(20pt)[Joint work with Shidong Jiang, Jiuyang Liang, Zhenli Xu, and Qi Zhou]

#text(20pt)[arXiv:2412.04595]

== Background

#pagebreak()

=== Quasi-2D charged systems

#timecounter(1)

Quasi-2D systems @mazars2011long are at the macroscopic scale in $x y$, but microscopic in $z$, so that are always modeled as doubly periodic in numerical simulations.
// Q2D systems are widely exist in nature and engineering, for example, cell membranes and electrolyte near surfaces.

#figure(
  image("figs/Q2D.png", width: 400pt),
  caption: [Illustration of a quasi-2D charged system.],
)

Coulomb interaction plays a key role in nature, leading to effect such as ion transportation and self-assembly @luijten2014.

However, the Coulomb interaction decays as $r^(-1)$ in 3D, so that it is long ranged and singular at $r=0$, which make such simulation computationally expensive.

// #pagebreak()
// === Coulomb interaction

// #timecounter(1)

// Coulomb interaction plays a key role in nature, leading to effect such as ion transportation and self-assembly.

// #figure(
//   image("figs/self-assembly.png", width: 450pt),
//   caption: [Self-assembly of nanoparticles via Coulomb interaction @luijten2014.
// ],)

// However, the Coulomb interaction decays as $r^(-1)$ in 3D, so that is long ranged and singular at $r=0$, which make such simulation computationally expensive.
// Complexity of a direct sum of the Coulomb interaction in doubly periodic systems is about $O(N^2 epsilon^(-1/3))$.

#pagebreak()

=== Algorithms for Q2D charged systems
#timecounter(1)

Methods have been developed to accelerate the Coulomb interaction in Q2D systems.

The very first method is the Ewald2D @parry1975electrostatic method based on the Ewald splitting of the Coulomb kernel. It is accurate but with $O(N^2)$ complexity.

To reduce the complexity, most methods rely on the following three strategies:

- *Fourier spectral method* (Lindbo & Tornberg, 2011; 2012; Nestler et al., 2015; Shamshirgar & Tornberg, 2017; Shamshirgar et al., 2021; Maxian et al., 2021): based on Ewald splitting and fast Fourier transform (FFT), with $O(N log N)$ complexity.

// @lindbo2011spectral @lindbo2012fast, @nestler2015fast, @shamshirgar2017spectral, @shamshirgar2021fast, @maxian2021fast

- *Fast multipole methods* (Greengard & Rokhlin, 1987; Berman & Greengard, 1994; Yan & Shelley, 2018; Liang et al., 2020): accelerated by hierarchical low-rank compression, adaptive and with $O(N)$ complexity.

// @greengard1987fast, @yan2018flexibly, @liang2020harmonic, @liang2022hsma, @jiang2023jcp, @berman1994renormalization

- *Random batch Ewald* (Jin et al., 2021; Liang et al., 2022; Gan et al., 2024a; Gan et al., 2024b): based on Ewald splitting and random batch sampling, stochastic and with $O(N)$ complexity, efficient parallelization.

// @Jin2020SISC, @gan2024fast, @gan2024random, @liang2022superscalability

#pagebreak()

=== Algorithms for Q2D charged systems
#timecounter(1)

For doubly periodic systems, one major challenge is the large prefactor in $O(N)$ or $O(N log N)$ compared to 3D-PBC solvers @mazars2011long, especially when the system is strongly confined in the $z$ direction, i.e., $L_z << L_x, L_y$.

- For the FFT based methods, *huge zero-padding* is required.

- For the FMM based methods, *more near field contributions* is needed.

Some recently developed methods offer potential solutions to this challenge, including:
- Anisotropic truncation kernel method @greengard2018anisotropic
- Periodic FMM @jiang2023jcp
- Dual-space multilevel kernel-splitting method @greengard2023dual

However, these methods have not yet been extended to handle quasi-2D systems.

== The algorithm

=== The sum-of-Gaussians approximation

#timecounter(1)

In our work, we use the bilateral series approximation @beylkin2010approximation of the Coulomb kernel, where
$
  1 / r approx (2 log b) / sqrt(2 pi sigma^2) sum_(l = - infinity)^(infinity) 1 / (b^l) e^(- r^2 / (sqrt(2) b^l sigma)^2), " with" cal(E)_r < 2 sqrt(2) e^(- pi^2 / (2 log b)), r >0
$

Based on the u-series decomposition @DEShaw2020JCP, we further split the potential into three parts:
$
  1 / r approx underbrace(( 1 / r - sum_(l = 0)^M w_l e^(- r^2 / s_l^2)) bb(1)_(r < r_c), "near-field")+ underbrace(sum_(l = 0)^m w_l e^(- r^2 / s_l^2), "mid-range") + underbrace(sum_(l = m + 1)^M w_l e^(- r^2 / s_l^2), "long-range")
$
The weight of the narrowest Gaussian is modified to be
$
  w_0 = omega (2 log b) / (sqrt(2 pi sigma^2))
$
to enforce the $C^0$ and $C^1$ continuity of the near-field potential at $r = r_c$, which is important for MD simulations @shamshirgar2019regularizing.

// #pagebreak()

// === Near-field potential

// #timecounter(1)

// The near-field potential is computed by a real space truncation.

#pagebreak()

=== Splitting the far-field potential

#timecounter(1)

// How to determine the turning point $m$? 

Selecting $m$ so that $s_m < eta L_z < s_(m+1)$, where $eta$ is $O(1)$ constant.

#figure(
  image("figs/farfield.svg", width: 400pt),
  caption: [Mid-range part and long-range part of the potential, $L_z = 10$, $eta approx 0.6$.],
)

#pagebreak()

=== Mid-range potential

#timecounter(1)

The mid-range potential is computed by a standard Fourier spectral solver #footnote(text(12pt)[#link("https://github.com/HPMolSim/ChebParticleMesh.jl")],) (type-1 and type-2 NUFFT in 3D @barnett2019parallel) with *little zero padding* ($lambda_z < 2$ for double precision).
$
  Phi_("mid")^l (arrow(r)) = sum_(arrow(n)) sum_(j = 1)^N q_j w_l e^(- (arrow(r) - arrow(r)_j + arrow(n) circle arrow(L))^2 / s_l^2), quad s_l < eta L_z
$

#align(center, canvas({
  import draw: *
  let box_loc = 3.5
  let box_1 = (-2 * box_loc, 0)
  let box_2 = (-box_loc, 0)
  let box_3 = (0, 0)
  let box_4 = (box_loc, 0)
  let box_5 = (2 * box_loc, 0)

  content(box_1, box(text(15pt)[Gridding], stroke: black, inset: 10pt), name: "box1")
  content(box_2, box(text(15pt)[FFT], stroke: black, inset: 10pt), name: "box2")
  content(box_3, box(text(15pt)[Scaling], stroke: black, inset: 10pt), name: "box3")
  content(box_4, box(text(15pt)[IFFT], stroke: black, inset: 10pt), name: "box4")
  content(box_5, box(text(15pt)[Gathering], stroke: black, inset: 10pt), name: "box5")

  line("box1", "box2", mark: (end: "straight"))
  line("box2", "box3", mark: (end: "straight"))
  line("box3", "box4", mark: (end: "straight"))
  line("box4", "box5", mark: (end: "straight"))
}))


- *No need of the kernel truncation* in the free direction due to the smoothness and separability of the Gaussian.
- *No upsampling* is needed since the Fourier transform of the Gaussian decays quickly and it compensates the loss of accuracy in calculating the Fourier transform of the data.

#pagebreak()

=== Long-range potential

#timecounter(1)

The long-range potential is computed by a Fourier-Chebyshev solver.
$
  Phi_("long")^l (arrow(r)) = sum_(arrow(n)) sum_(j = 1)^N q_j w_l e^(- (arrow(r) - arrow(r)_j + arrow(n) circle arrow(L))^2 / s_l^2), quad s_l > eta L_z
$

The extremely smooth long-range Gaussians are interpolated on the Chebyshev proxy points in $z$, similar to that of the periodic FMM @jiang2023jcp, and only *$O(1)$ number of Chebyshev points are required*.

Then 2D NUFFT can be used to evaluate the potential on a tensor-product grid.

#align(center, canvas({
  import draw: *
  let box_loc = 6
  let box_1 = (-box_loc, 0)
  let box_2 = (0, 0)
  let box_3 = (box_loc, 0)

  content(box_1, box(text(15pt)[Chebyshev interpolation], stroke: black, inset: 10pt), name: "box1")
  content(box_2, box(text(15pt)[2D NUFFT], stroke: black, inset: 10pt), name: "box2")
  content(box_3, box(text(15pt)[Evaluating polynomial], stroke: black, inset: 10pt), name: "box3")

  line("box1", "box2", mark: (end: "straight"))
  line("box2", "box3", mark: (end: "straight"))

}))

In cubic systems, $L_x ~ L_y ~ L_z$, $O(1)$ Fourier modes in $x y$ and $O(1)$ Chebyshev points in $z$, no need for NUFFT.

#pagebreak()

=== Complexity

#timecounter(1)

Using DFT for long-range potential, the complexity is
$
  O"("underbrace(4 pi r_c^3 rho_r N, "near-field") + underbrace(cal(P)_x cal(P)_y cal(P)_z N + (lambda_z (1 + delta / L_z)) / (r_c^3 rho_r) N log N, "mid-range by 3D-NUFFT") + underbrace((P L_x L_y) / (eta^2 L_z^2) N, "long-range by DFT") ")"
$
where $cal(P)_x, cal(P)_y, cal(P)_z$ are the window supports, $lambda_z$ is the padding ratio, $delta$ is the extended length of the box in the free direction to accommodate the support of the window function, $P$ is the number of Chebyshev points.
By taking $r_c ~ O(1)$ and assume $L_z ~ O(sqrt(L_x L_y))$, the complexity is $O(N log N)$.

Using 2D-NUFFT for long-range potential, the complexity is
$
  O"("underbrace(4 pi r_c^3 rho_r N, "near-field") + underbrace(cal(P)_x cal(P)_y cal(P)_z N + (lambda_z (1 + delta / L_z)) / (r_c^3 rho_r) N log N, "mid-range by 3D-NUFFT") + underbrace(cal(P)_x cal(P)_y P N +  (P L_x L_y) / (s_(m + 1)^2) N log N, "long-range by 2D-NUFFT") ")"
$
which is needed when $L_z << L_x, L_y$, the total complexity is also $O(N log N)$.

== Numerical results

#timecounter(2)

The method #footnote(text(12pt)[#link("https://github.com/HPMolSim/FastSpecSoG.jl")],) is benchmarked on the following systems:
- Cubic systems with fixed aspect ratio equals to $1$.
- Strongly confined systems with fixed $L_z$, aspect ratio up to $10^(3.5)$

#figure(
  image("figs/sog_benchmark.png", width: 370pt),
  caption: [#text(15pt)[Error and time cost for the SOG method in the (a,c) cubic and (b,d) strongly confined systems.]],
)

#pagebreak()

== Summary

#timecounter(1)

A fast and accurate solver for Q2D charged systems is developed based on the sum-of-Gaussian approximation of the Coulomb kernel and the kernel splitting technique.
The method can be regarded as a 2-level DMK method @greengard2023dual.
// The Coulomb kernel is splitted into three parts:
// - near field terms: solved by real space truncation
// - mid-range terms: solved by the Fourier spectral method with little zero padding and no upsampling
// - long-range terms: solved by the Fourier-Chebyshev method with $O(1)$ number of Chebyshev points

It has the following advantages:
- spectrally accurate with rigorous error analysis @liang2025errorestimateuseriesmethod
- need little/no zero-padding for systems that are confined in a rectangular box of high aspect ratio
- no need for upsampling in the gridding step
- all calculations are carried out in the fundamental cell itself
- easy to be implemented and parallelized for large-scale MD simulations

Currently, the major shortcoming of this method is its non-adaptive nature, and has a complexity of $O(N log N)$ rather than $O(N)$.

= Automated Discovery of the Optimal Branching Rules 

#text(20pt)[Joint work with Yi-Jia Wang, Pan Zhang, and Jin-Guo Liu]

#text(20pt)[arXiv:2412.07685]

== Background

=== The maximum independent set problem
#timecounter(1)

#align(center, box([One of the first batch of 21 NP-hard problems proved by @Karp1972.], stroke: black, inset: 10pt))

An independent set is a set of vertices in a graph, no two of which are adjacent. 

#align(center,
  grid(columns: 3,
    align(center, 
    canvas(length: 25pt, {
      import draw: *
      let s = 3
      let dy = 3.0
      let la = (-s, 0)
      let lb = (0, s)
      let lc = (0, 0)
      let ld = (s, 0)
      let le = (s, s)
      for (loc, name, color) in ((la, "a", red), (lb, "b", black), (lc, "c", red), (ld, "d", black), (le, "e", red)) {
        circle(loc, radius:0.5, name: name, fill: color)
        content(loc, text(20pt, fill:white)[$#name$])
      }
      for (a, b) in (("a", "b"), ("b", "c"), ("c", "d"), ("d", "e"), ("b", "d")) {
        line(a, b)
      }
      content((1.5, -1.5), text(15pt)[$G = (V, E)$, MIS = ${a, c, e}$, size $alpha(G) = 3$])
    })),
    h(50pt),
    align(center, image("figs/ob_rydberg.png", width: 250pt))
  )
)
MIS problem has an exponential large solution space, since it is NP-hard, no polynomial-time algorithm is known to solve it exactly.

It caught much attention since the Rydberg atom systems realize spin models that naturally MIS problem @Nguyen2023QuantumOptimization.

// The branching algorithms @Fomin2013 is the most widely used method to solve this problem.

#pagebreak()

=== Branching algorithm

#timecounter(2)

The branching algorithm @Fomin2013 explores the solution space using a tree-like structure, relying on predesigned rules.

Complexity of a branching algorithm is always described as *$O(gamma^n)$* where $gamma$ is the branching factor and $n$ is the size of the problem.

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
    circle((l.at(0), l.at(1)-s/2), radius:0.4, name: n, fill: black)
    content((l.at(0), l.at(1)-s/2), text(14pt, fill:white)[$#n$])
  }
  for (a, b) in (("a", "b"), ("b", "c"), ("c", "d"), ("d", "e"), ("b", "d")){
    line(a, b)
  }
  scircle((-4, -dy), (2, 1.5), "brancha")
  for (l, n) in ((lc, "c"), (ld, "d"), (le, "e")){
    let loc = (l.at(0)-5, l.at(1)-s/2-dy)
    circle(loc, radius:0.4, name: n, fill: black)
    content(loc, text(14pt, fill:white)[$#n$])
  }
  for (a, b) in (("c", "d"), ("d", "e"), ("c", "d")){
    line(a, b)
  }
  scircle((4, -dy), (1, 1), "branchb")
  circle((4, -dy), radius:0.4, name: "e", fill: black)
  content((4, -dy), text(14pt, fill:white)[$e$])
  scircle((-6, -2*dy), (1, 1), "branchaa")
  circle((-6, -2*dy), radius:0.4, name: "e", fill: black)
  content((-6, -2*dy), text(14pt, fill:white)[$e$])
  scircle((-2, -2*dy), (0.5, 0.5), "branchab")
  content((-2, -2*dy), text(14pt)[$2$])
  scircle((4, -2*dy), (0.5, 0.5), "branchba")
  content((4, -2*dy), text(14pt)[$2$])
  scircle((-6, -3*dy), (0.5, 0.5), "branchaaa")
  content((-6, -3*dy), text(14pt)[$3$])
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
  s[1977], s[$O^*(1.2600^n)$], s[@Tarjan1977], s[],
  s[1986], s[$O^*(1.2346^n)$], s[@Jian1986], s[],
  s[1986], s[$O^*(1.2109^n)$], s[@Robson1986], s[],
  s[1999], s[$O^*(1.0823^m)$], s[@Beigel1999], s[num of edges],
  s[2001], s[$O^*(1.1893^n)$], s[@Robson2001], s[],
  s[2003], s[$O^*(1.1254^n)$ for 3-MIS], s[@Chen2003], s[],
  s[2005], s[$O^*(1.1034^n)$ for 3-MIS], s[@Xiao2005], s[],
  s[2006], s[$O^*(1.2210^n)$], s[@Fomin2006], s[],
  s[2006], s[$O^*(1.1225^n)$ for 3-MIS], s[@Fomin2006b], s[],
  s[2006], s[$O^*(1.1120^n)$ for 3-MIS], s[@Furer2006], s[],
  s[2006], s[$O^*(1.1034^n)$ for 3-MIS], s[@Razgon2006], s[],
  s[2008], s[$O^*(1.0977^n)$ for 3-MIS], s[@Bourgeois2008], s[],
  s[2009], s[$O^*(1.0919^n)$ for 3-MIS], s[@Xiao2009], s[],
  s[2009], s[$O^*(1.2132^n)$], s[@Kneis2009], s[],
  s[2013], s[$O^*(1.0836^n)$ for 3-MIS], s[@Xiao2013], s[#highlight[SOTA]],
  s[2016], s[$O^*(1.2210^n)$], s[@Akiba2016], s[#highlight[PACE winner]],
  s[2017], s[$O^*(1.1996^n)$], s[@Xiao2017], s[SOTA],
),
)
)

#pagebreak()

// === Tensor networks
// #timecounter(1)

// Tensor network (TN) is a powerful tool to represent and manipulate high-dimensional data, and have been used in various fields such as quantum physics @nielsen2010quantum and machine learning @stoudenmire2016supervised.
// // Tensor networks have proven to be an powerful tool for solving large-scale problems.
// // In 2019, Google announced quantum supremacy with their Sycamore chip @arute2019quantum; however, by 2022, their results were classically simulated using tensor network-based methods @pan2022simulation.

// #align(center,
// grid(columns: 3,
//   figure(image("figs/big_batch.png", height: 280pt), caption: [
//     #text(15pt)[Classical simulation of Google's Sycamore chip via tensor network @pan2022simulation.]
//   ]),
//   h(20pt),
//   figure(image("figs/tn_inference.png", height: 280pt), caption: [
//     #text(15pt)[Exact solve large-scale PACE inference problem via tensor network @Roa2024.]
//   ]),
// )
// )

// #pagebreak()

== The algorithm

=== Tensor networks for the MIS problem
#timecounter(2)

// The MIS problem can be efficiently solved and analyzed by tensor networks @ebadi2022quantum  @liu2023computing.

// #align(center,
// grid(columns: 5,
//   align(center + horizon, image("figs/kingsubgraph.png", width:150pt)),
//   h(30pt),
//   align(center + horizon, text(40pt)[$arrow$]),
//   h(30pt),
//   align(center + horizon, image("figs/tn_mis_solutionspace.png", width:400pt)),
// )
// )

Tensor networks can be used to extract the local information of the sub-graph @gao2024programmingguidesolvingconstraint @liu2023computing.

#align(center, 
  grid(
    columns: 5,
    align(center + horizon, 
      // canvas({
      //   import draw: *

      //   circle((0, 0), radius: 2.7, fill: rgb("#63e64269"))
      //   circle((0, -0.2), radius: 1.9, fill: white, stroke: white)

      //   let la = (0, 1)
      //   let lb = (-1, -1)
      //   let lc = (1, -1)
      //   let ld = (0, 0)
      //   let le = (1, 0)
      //   let s = 1.7
      //   let lao = (0, 2.0)
      //   let lbo = (-s, -s)
      //   let lco = (s, -s)

      //   line(la, lao, stroke: (dash : "dashed"))
      //   line(lb, lbo, stroke: (dash : "dashed"))
      //   line(lc, lco, stroke: (dash : "dashed"))

      //   for (l, n) in ((la, "a"), (lb, "b"), (lc, "c"), (ld, "d"), (le, "e")){
      //     circle((l.at(0), l.at(1)), radius:0.3, name: n, fill: black)
      //     content((l.at(0), l.at(1)), text(14pt, fill:white)[$#n$])
      //   }
      //   for (a, b) in (("a", "e"), ("b", "c"), ("c", "d"), ("d", "a"), ("b", "d"), ("d", "e"), ("c", "e")){
      //     line(a, b)
      //   }

      //   content((-1.2, 2), text(16pt, fill:black)[$G$])
      //   content((-1, 1), text(16pt, fill:black)[$R$])
      // })
      image("figs/ob_graph.png", width: 150pt)
    ),
    h(30pt),
    align(center + horizon, 
      canvas({
        import draw: *
        let a = (-2, 0)
        let b = (2, 0)
        line(a, b, mark: (end: "straight"))
        content((0, -1.2), text(15pt)[Contract the local \ tensor network])
      })
    ),
    h(30pt),
    table(
      columns: (auto, auto),
      table.header(hdl[Boundary configuration: $s_(a b c)$], hdl[Possible assignments: $S_(a b c d e)$]),
      sl[000], sl[00001, 00010],
      sl[001], sl[00101],
      sl[010], sl[01010],
      sl[111], sl[11100],
    )
  )
)

However, a pure tensor network method does not work well for non-geometric graphs.
Its complexity on 3-regular graphs is about $O(1.1224^n)$, far from the SOTA ($O^*(1.0836^n)$).

#pagebreak()

// === What is the difference?

// #timecounter(2)

//   #leftright(
//     text(20pt)[Tensor network approach \ 1. Contract the local tensor network for the sub-graph and pick the non-zero elements. \ 2. For all possible boundaries, fix all the variables and continue the contraction$ gamma^n = 3 times gamma^(n-5) arrow gamma approx 1.2457 $],
//     text(20pt)[Branching algorithm \ 1. Search for structures in the sub-graph. \ 2. Find $d$ and $e$ are connected and $N[e] subset N[d]$, satisfying the domination rule, fix $d = 0$, i.e., not in the MIS $ gamma^n = gamma^(n-1) arrow gamma = 1.0 $]
//   )

// #align(center, 
//   grid(
//     columns: 3,
//     align(center + horizon, 
//       canvas({
//         import draw: *

//         circle((0, 0), radius: 2.7, fill: rgb("#63e64269"))
//         circle((0, -0.2), radius: 1.9, fill: white, stroke: white)

//         let la = (0, 1)
//         let lb = (-1, -1)
//         let lc = (1, -1)
//         let ld = (0, 0)
//         let le = (1, 0)
//         let s = 1.7
//         let lao = (0, 2.0)
//         let lbo = (-s, -s)
//         let lco = (s, -s)

//         line(la, lao, stroke: (dash : "dashed"))
//         line(lb, lbo, stroke: (dash : "dashed"))
//         line(lc, lco, stroke: (dash : "dashed"))

//         for (l, n) in ((la, "a"), (lb, "b"), (lc, "c"), (ld, "d"), (le, "e")){
//           circle((l.at(0), l.at(1)), radius:0.3, name: n, fill: black)
//           content((l.at(0), l.at(1)), text(14pt, fill:white)[$#n$])
//         }
//         for (a, b) in (("a", "e"), ("b", "c"), ("c", "d"), ("d", "a"), ("b", "d"), ("d", "e"), ("c", "e")){
//           line(a, b)
//         }

//         content((-1.2, 2), text(16pt, fill:black)[$G$])
//         content((-1, 1), text(16pt, fill:black)[$R$])
//       })
//     ),
//     h(30pt),
//     table(
//       columns: (auto, auto),
//       table.header(hdl[Boundary configuration: $s_(a b c)$], hdl[Possible assignments: $S_(a b c d e)$]),
//       sl[000], sl[00010, 000#redtext(0)1],
//       // sl[100], sl[100#redtext(0)0],
//       sl[010], sl[010#redtext(0)1],
//       // sl[001], sl[001#redtext(0)0],
//       // sl[110], sl[110#redtext(0)0],
//       sl[101], sl[101#redtext(0)0]
//     )
//   )
// )

// #align(center, box([Key point: No need to use all results, find the #emph[correct pattern]!], stroke: black, inset: 10pt))

// #pagebreak()

=== The optimal branching algorithm

#timecounter(1)

We use the tensor network to extract the local information, and then automatically search the optimal branching rules #footnote(text(12pt)[#link("https://github.com/OptimalBranching/OptimalBranching.jl")],) @Gao2024.

#align(center, image("figs/ob_rules.png", width: 600pt))


#align(center,
grid(
  columns: 3,
  text(20pt, black)[Naive branching \ 4 branches, each fix 5 variables \ $gamma^n = 4 times gamma^(n-5)$ \ $gamma approx 1.3195$],
  h(60pt),
  text(20pt, black)[Optimal branching \ 3 branches, fix [4, 5, 5] variables \ $gamma^n = gamma^(n-4) + 2 times gamma^(n-5)$ \ $gamma approx 1.2671$],
)
)

#align(center, box([Question: How to solve the rules from the assignments?], stroke: black, inset: 10pt))

#pagebreak()

=== Finding the optimal branching rule

#timecounter(2)

Bruteforce search? $arrow$ Given $l$ assignments, possible rules: $O(2^(2^l))$.

The process of finding the optimal branching rules is as the following:

#grid(columns: 3,
  align(center, canvas({
    import draw: *
    content((0, 0), box(text(15pt)[Possible assignments $ {bold(s)_1, bold(s)_2, dots, bold(s)_l} $], stroke: black, inset: 10pt), name: "oracle")
    content((0, -4.5), box(text(15pt)[Candidate clauses $ cal(C) = {c_1, c_2, dots, c_m} $], stroke: black, inset: 10pt), name: "clauses")
    content((0, -9), box(text(15pt)[Optimal branching rule $ cal(D) = c_(k_1) or c_(k_2) or dots $], stroke: black, inset: 10pt), name: "branching")
    line("oracle", "clauses", mark: (end: "straight"))
    line("clauses", "branching", mark: (end: "straight"))
  })),
  h(40pt),
  text(20pt)[
    The first step is very direct forward, we generate all possible combination of the assignments (the candidate clauses).

    The second step is formulated as a *set covering problem*, which can be solved by mixed integer programming solvers @Achterberg2009.
    $
      min_(gamma, bold(x)) gamma " s.t. " & sum_(i=1)^m gamma^(-Delta rho(c_i)) x_i = 1,\
      & union.big_(i = 1, dots, m\ x_i = 1) J_i = {1, 2, dots, l},  #h(10pt) arrow "valid branching rule"\
      & x_i in {0, 1} #h(10pt) arrow "a clause is selected or not"
    $
    where $Delta rho(c_i)$ is the size reduced by the clause $c_i$ of the problem.
  ]
)

// The candidate clauses are the combinations of the possible assignments, for example, a combination of the 3rd and 4th assignments is given by:

// $
//   "combine"(not a and b and not c and d and not e, a and b and c and not d and not e) = b and not e
// $

// we say $b and not e$ covers ${3, 4}$.
// The clauses are generated iteratively.

// #pagebreak()

// === Set covering

// #timecounter(1)

// Then we formulate the problem as a set covering problem, which can be solved by MIP solvers @Achterberg2009.

// $
// min_(gamma, bold(x)) gamma " s.t. " & sum_(i=1)^m gamma^(-Delta rho(c_i)) x_i = 1,\
// & union.big_(i = 1, dots, m\ x_i = 1) J_i = {1, 2, dots, l},  #h(50pt) arrow "valid branching rule"\
// & x_i in {0, 1} #h(161pt) arrow "a clause is selected or not"
// $
// where $Delta rho(c_i)$ is the size reduced by the clause $c_i$ of the problem.

// The chosen clauses form the optimal branching rule:
// $
//   cal(D) = c_(k_1) or c_(k_2) or dots or c_(k_m)
// $
// with the minimum $gamma$ among all the possible branching rules.

== Numerical results

=== A bottleneck case

#timecounter(1)

A bottle neck case has been reported in Xiao's work @Xiao2013, with a branching factor of $1.0836$.

#grid(columns: 3,
  image("figs/ob_bottleneck.png", width: 300pt),
  h(30pt),
  align(horizon, text(20pt, black)[
    - 71 possible assignments, 15782 candidate clauses. 
    // - The optimal branching rule can be solved in few seconds.
    - 4 branches, size reduced by branches: $[10, 16, 26, 26]$, with 
    *$ gamma = 1.0817 < 1.0836 $*
    which indicates our method can find better branching rules than the predesigned rules.
  ]),
)

#pagebreak()

=== Branchmark on random graphs
#timecounter(1)

#leftright(
  grid(rows: 5,
  align(center,
    canvas(length: 25pt, {
      import draw: *
      content((0, 0), box(text(15pt)[Reduction by rules], stroke: black, inset: 10pt), name: "reduce")
      content((0, -2), box(text(15pt)[Selecting subgraph], stroke: black, inset: 10pt), name: "select")
      content((0, -4), box(text(15pt)[Generating branching rules], stroke: black, inset: 10pt), name: "rules")
      content((0, -6), box(text(15pt)[Branching], stroke: black, inset: 10pt), name: "branching")
      // content((0, -8), ortho({on-xz({rect((-1.5,-1.5), (1.5,1.5))})}))

      let ps_1 = (0, -7.5)
      let ps_2 = (0, -9.5)
      let ps_3 = (-2, -8.5)
      let ps_4 = (2, -8.5)

      line(ps_1, ps_3)
      line(ps_3, ps_2)
      line(ps_2, ps_4)
      line(ps_4, ps_1)

      content((0, -8.5), text(15pt)[Is solved?], stroke: black, inset: 10pt)

      let point_1 = (5, -8.5)
      let point_2 = (5, 0)
      let point_start = (0, 1.2)
      let point_end = (0, -10.0)
      line("reduce", "select", mark: (end: "straight"))
      line("select", "rules", mark: (end: "straight"))
      line("rules", "branching", mark: (end: "straight"))
      line(ps_4, point_1)
      line(point_1, point_2)
      line(point_2, "reduce.east", mark: (end: "straight"))
      line(point_start, "reduce", mark: (end: "straight"))
      line("branching", ps_1, mark: (end: "straight"))
      line(ps_2, point_end, mark: (end: "straight"))
    })
  ),
  v(10pt),
  text(15pt)[The resulting methods are denoted as *ob* and *ob+xiao* and the average branching factor is shown in the table.],
  v(10pt),
  align(center, table(
    columns: (auto, auto, auto, auto, auto, auto),
    table.header(hd[], hd[*ob*], hd[*ob+\ xiao*], hd[xiao2013], hd[akiba2015], hd[akiba2015+\ xiao&packing]),
    s[3RR], s[1.0457], s[*1.0441*], s[*1.0487*], s[-], s[-],
    s[ER], s[1.0011], s[1.0002], s[-], s[1.0044], s[1.0001],
    s[KSG], s[1.0116], s[1.0022], s[-], s[1.0313], s[1.0019],
    s[Grid], s[1.0012], s[1.0009], s[-], s[1.0294], s[1.0007],
  ))),
  figure(image("figs/ob_benchmark_average.png", width: 380pt), caption : [#text(15pt)[Average number of branches generated by different branching algorithms on 1000 random graphs.]])
)

#pagebreak()

== Potential applications

=== Sparse Tensor Networks Contraction via Optimal Branching
#timecounter(2)

The optimal branching algorithm can be applied to contract the sparse tensor networks.
// Assume that $T$ is a sparse tensor, the non-zeros values are listed as below:

#align(center,
grid(columns: 9,
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
    // show-graph-content(locs.map(v => (v.at(0) * s, v.at(1) * s)), e, c, radius:0.5, fontsize: 16pt)
    content((T.at(0) * s, T.at(1) * s), text(16pt, fill:white)[$T$])
    for (i, v) in mid.map(v => (v.at(0) * (s), v.at(1) * (s))).enumerate() {
      circle(v, radius:0.3, fill:white, stroke:none)
      content(v, text(15pt, black)[#indices.at(i)])
    }
  }),
  h(20pt),
  align(center + horizon,text(30pt)[$arrow$]),
  h(20pt),
  align(center + horizon,
  canvas({
    import draw: *
    content((0, 3.5), text(20pt, black)[$T$])
    content((0,0), text(20pt, black)[
      #table(
      columns: (auto, auto),
      inset: 8pt,
      align: horizon,
      table.header(hdl[$i j k l$], hdl[value]),
      sl[#redtext(11)01], sl[0.1],
      sl[#redtext(11)10], sl[0.2],
      sl[#redtext(11)#bluetext(0)#bluetext(0)], sl[0.3],
      sl[10#bluetext(0)#bluetext(0)], sl[0.4],
      sl[01#bluetext(0)#bluetext(0)], sl[0.5]
      )
    ])
  })),
  h(20pt),
  align(center + horizon,text(30pt)[$arrow$]),
  h(20pt),
  grid(rows : 5,
    canvas(length: 25pt, {
      import draw: *
      let T = (0, 0)
      let A = (1, 0)
      let B = (0, 1)
      let C = (-1, 0)
      let D = (0, -1)
      let locs = (T, A, B, C, D)
      let mid = (A, B, C, D).map(v => (v.at(0) * 0.5, v.at(1) * 0.5))
      let e = ((0, 1), (0, 2), (0, 3), (0, 4))
      let e2 = ((0, 3), (0, 4))
      let c = ((0, $T_(1 1 k l)$), (1, $A_(1 *)$), (2, $B_(1 *)$), (3, $C_(k *)$), (4, $D_(l *)$))
      let s = 2
      let indices = (" ", " ", $k$, $l$)
      show-graph-black(locs.map(v => (v.at(0) * (s + 0.8), v.at(1) * (s + 0.8))), e2, radius:0.0)
      show-graph-black(locs.map(v => (v.at(0) * s, v.at(1) * s)), e2, radius:0.5)
      // show-graph-content(locs.map(v => (v.at(0) * s, v.at(1) * s)), e, c, radius:0.5, fontsize: 10pt)
      // content((T.at(0) * s, T.at(1) * s), text(10pt, fill:white)[$T_1$])
      for (i, v) in mid.map(v => (v.at(0) * (s), v.at(1) * (s))).enumerate() {
        circle(v, radius:0.3, fill:white, stroke:none)
        content(v, text(10pt, black)[#indices.at(i)])
      }
      line((B.at(0) * s, B.at(1) * s), (B.at(0) * s, B.at(1) * s + 1))
      line((A.at(0) * s, A.at(1) * s), (A.at(0) * s + 1, A.at(1) * s))
    }),
    v(10pt),
    text(20pt)[$+$],
    v(10pt),
    canvas(length: 25pt, {
      import draw: *
      let T = (0, 0)
      let A = (1, 0)
      let B = (0, 1)
      let C = (-1, 0)
      let D = (0, -1)
      let locs = (T, A, B, C, D)
      let mid = (A, B, C, D).map(v => (v.at(0) * 0.5, v.at(1) * 0.5))
      let e = ((0, 1), (0, 2), (0, 3), (0, 4))
      let e2 = ((0, 1), (0, 2))
      let c = ((0, $T_(i j 0 0)$), (1, $A_(i *)$), (2, $B_(j *)$), (3, $C_(0 *)$), (4, $D_(0 *)$))
      let s = 2
      let indices = ($i$, $j$, " ", " ")
      show-graph-black(locs.map(v => (v.at(0) * (s + 0.8), v.at(1) * (s + 0.8))), e2, radius:0.0)
      show-graph-black(locs.map(v => (v.at(0) * s, v.at(1) * s)), e2, radius:0.5)
      // show-graph-content(locs.map(v => (v.at(0) * s, v.at(1) * s)), e, c, radius:0.5, fontsize: 10pt)
      // content((T.at(0) * s, T.at(1) * s), text(10pt, fill:white)[$T_2$])
      for (i, v) in mid.map(v => (v.at(0) * (s), v.at(1) * (s))).enumerate() {
        circle(v, radius:0.3, fill:white, stroke:none)
        content(v, text(10pt, black)[#indices.at(i)])
      }
      line((C.at(0) * s, C.at(1) * s), (C.at(0) * s - 1, C.at(1) * s))
      line((D.at(0) * s, D.at(1) * s), (D.at(0) * s, D.at(1) * s - 1))
    })
  )
)
)

// Such a contraction can be viewed as a branching problem, one can use the optimal branching algorithm to find the optimal way. 

Such sparsity is common in many problems, including *probabilistic inference*, *combinatorial optimization*, and *quantum circuit simulations* @Markov2008.

== Summary

===

#timecounter(1)

A new method to automatically discover the optimal branching rules is proposed, by combining the tensor network method and the branching algorithm.

// With this method, we achieved an average complexity of $O(1.0441^n)$ on random 3-regular graphs, which outperforms the SOTA ($O(1.0487^n)$).

Advantages:
- generate the branching rules automatically without human effort
- fully utlize the information of the sub-graph
- the sub-graph can be selected flexibly
- can be applied to different problems, not only the MIS problem

Disadvantages:
- solving the rule can be computationally expensive
- cannot capture the rules need graph rewriting

= Summary and Outlook

== Fast Summation Algorithms
#timecounter(1)


=== Publications

- *X. Gao*, S. Jiang, J. Liang, Z. Xu, and Q. Zhou, A fast spectral sum-of-Gaussians method for electrostatic summation in quasi-2D systems, Arxiv:2412.04595 (2024)
- Z. Gan, *X. Gao*, J. Liang, and Z. Xu, Fast algorithm for quasi-2D Coulomb systems, Arxiv:2403.01521 (2024). Accepted by Journal of Computational Physics.
- Z. Gan, *X. Gao*, J. Liang, and Z. Xu, Random batch Ewald method for dielectrically confined Coulomb systems, Arxiv:2405.06333 (2024)
- *X. Gao* and Z. Gan, Broken symmetries in quasi-2D charged systems via negative dielectric confinement, The Journal of Chemical Physics 161, (2024)

=== Software packages

// - *ExTinyMD.jl*: A simple framework for MD simulation.
- *ChebParticleMesh.jl*#footnote(text(12pt)[#link("https://github.com/HPMolSim/ChebParticleMesh.jl")],): Toolkits for smooth particle mesh (type-1 and type-2 NUFFT).
- *FastSpecSoG.jl*#footnote(text(12pt)[#link("https://github.com/HPMolSim/FastSpecSoG.jl")],): Implementation of the fast spectral SOG method.
- *EwaldSummations.jl*#footnote(text(12pt)[#link("https://github.com/HPMolSim/EwaldSummations.jl")],): Various Ewald summation methods with parallelization.

#pagebreak()

== Tensor Network Algorithms
#timecounter(1)

=== Publications

- *X. Gao*, Y.-J. Wang, P. Zhang, and J.-G. Liu, Automated discovery of branching rules with optimal complexity for the maximum independent set problem, Arxiv:2412.07685 (2024)
- *X. Gao*, X. Li, and J. Liu, Programming guide for solving constraint satisfaction problems with tensor networks, Arxiv:2501.00227 (2024)
- M. Roa-Villescas, *X. Gao*, S. Stuijk, H. Corporaal, and J.-G. Liu, Probabilistic inference in the era of tensor networks and differential programming, Physical Review Research 6, 33261 (2024)

=== Software packages

- *CuTropicalGEMM.jl*#footnote(text(12pt)[#link("https://github.com/TensorBFS/CuTropicalGEMM.jl")],): Custom GPU kernel for tropical matrix multiplication (supported by OSPP 2023).
- *TreeWidthSolver.jl*#footnote(text(12pt)[#link("https://github.com/ArrogantGao/TreeWidthSolver.jl")],): Solving the treewidth problem (supported by GSoC 2024).
// - *OMEinsumContractionOrders.jl*#footnote(text(12pt)[#link("https://github.com/TensorBFS/OMEinsumContractionOrders.jl")],): Optimizing the tensor network contraction order.
- *OptimalBranching.jl*#footnote(text(12pt)[#link("https://github.com/OptimalBranching/OptimalBranching.jl")],): Implementation of the optimal branching algorithm.

#pagebreak()

== Future Research Plans
#timecounter(1)

=== Fast Summation Algorithms

- Efficient methods based on the DMK framework @greengard2023dual
- GPU acceleration for the fast algorithms

\

=== Tensor Network Algorithms

// With the advanced tensor network techniques, I am interested in exploring more flexible tensor network structures as ansatz, which may yield more accurate representations of quantum many-body states.

// I am also interested in developing efficient algorithms for evolving these states by integrating them with the Dirac–Frenkel/McLachlan variational principle @raab2000dirac, the automatic differentiation and proper pre-conditioning @ganahl2017continuous.

- Branching based sparse tensor network contraction
- More flexible quantum many-body ansatz
// - proper pre-conditioning @ganahl2017continuous


// I am keen on developing efficient summation methods for long-range interacting systems, such as those described by Poisson’s equation and the Helmholtz equation. For example, I am interested in extending the recently introduced dual-space multilevel kernel-splitting (DMK) @jiang2023dmk framework to a variety of long-range interactions and systems with periodic/partially periodic boundary conditions.

#pagebreak()

== Acknowledgements
#timecounter(1)

#let figsize = 130pt

#align(center,
grid(columns: 5, 
  grid(rows : 3, image("photos/zechenggan.jpeg", width : figsize), "" ,  text[Prof. Zecheng Gan \ HKUST(GZ)]),
  h(50pt), 
  grid(rows : 3, image("photos/zhenlixu.jpg", width : figsize), "" ,  text[Prof. Zhenli Xu \ SJTU]),
  h(50pt), 
  grid(rows : 3, image("photos/shidongjiang.jpeg", width : figsize), "" ,  text[Prof. Shidong Jiang \ CCM]),
  )
)

#align(center,
grid(columns: 3, 
  grid(rows : 3, image("photos/jinguoliu.jpeg", width : figsize), "" ,  text[Prof. Jinguo Liu, \ HKUST(GZ)]),
  h(60pt), 
  grid(rows : 3, image("photos/panzhang.jpeg", width : figsize), "" ,  text[Prof. Pan Zhang, \ ITP, CAS]),
  )
)

Also thank Jiuyang Liang (SJTU & CCM), Qi Zhou (SJTU), Martin (TU/e), and Yijia Wang (ITP, CAS) for collaboration.

#focus-slide[
  Thank you for your attention!
]

#show: appendix

= Appendix <touying:unoutlined>

#pagebreak()

== The U-series and its derivative

#figure(
  image("figs/nearfield.svg", width: 800pt),
  caption: [The U-series and its derivative, $r_c = 10.0$.],
)

#pagebreak()

== SOG parameters

#align(center, image("figs/sog_table.png", width: 500pt))

#figure(
  image("figs/sog_totalerror.png", width: 500pt),
  caption: [U-series parameters and the error],
)

#pagebreak()

== Window function

#figure(
  image("figs/sog_3dnufft.png", width: 800pt),
  caption: [Different window functions.],
)

#pagebreak()

== Zero-padding

#figure(
  image("figs/sog_zeropadding.png", width: 800pt),
  caption: [Zero-padding.],
)

#pagebreak()

== Chebyshev interpolation

#figure(
  image("figs/sog_cheb.png", width: 800pt),
  caption: [Accuracy of the Fourier-Chebyshev solver.],
)

#pagebreak()

== Strongly confined systems

#figure(
  image("figs/sog_thin.png", width: 500pt),
  caption: [Strongly confined systems.],
)

#pagebreak()

== Einsum notation

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

== Contraction order

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

== Optimizing the contraction order

Finding the optimal contraction order is a NP-hard#footnote(text(12pt)[I.L. Markov, Y. Shi, SIAM J. Comput. 38, 963–981 (2008).],) problem! 
// However, it is lucky that we are always satisfied with a good enough solution, thus heuristic methods are sufficient in practice.

In the past few years, tools have been developed to optimize the contraction order:
- OMEinsumContractionOrder.jl#footnote(text(12pt)[#link("https://github.com/TensorBFS/OMEinsumContractionOrder.jl")],) in Julia
- Cotengra#footnote(text(12pt)[#link("https://github.com/jcmgray/cotengra")],) in Python

#figure(
  image("figs/tn_order.png", width: 280pt),
  caption: [Comparison of different contraction orders.],
)

#pagebreak()

== Tropical Tensor Network

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

== Tensor Network for Maximum Independent Set Problem

A tropical TN can be used to solve the MIS problem, a simple example is shown below:

#align(center, 
  grid(columns: 5,
  align(center + horizon, align(center, canvas({
  import draw: *
  for (loc, name, color) in (((0, 0), "a", black), ((0, 3), "b", black), ((3, 3), "d", black), ((3, 0), "c", black), ((5, 4.5), "e", black)) {
    circle(loc, radius:0.5, name: name, stroke: color)
    content(loc, [$#name$])
  }
  for (a, b) in (
    ("a", "b"),
    ("b", "d"),
    ("a", "c"),
    ("c", "d"),
    ("b", "c"),
    ("d", "e"),
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

== Solving the set covering problem via Mixed Integer Linear Programming

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

== Solving the set covering problem via Mixed Integer Linear Programming

#figure(
  image("figs/ob_iterate.png", width: 500pt),
  caption: [#text(15pt)[An example of solving the set covering problem via iterative MIP.]],
)

#pagebreak()

== The branching algorithm used in benchmark

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

== Benchmarks for the worst-case complexity

#figure(
  image("figs/ob_benchmark_worst.png", width: 500pt),
  caption: [#text(15pt)[The worst-case complexity of the proposed method on random graphs.]],
)

#pagebreak()

== References

#bibliography("ref.bib", style: "apa", title: none)
