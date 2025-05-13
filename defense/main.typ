#import "@preview/cetz:0.2.2": canvas, draw, tree, plot
#import "@preview/subpar:0.2.0"

#import "@preview/pinit:0.1.3": *
#import "@preview/ctheorems:1.1.3": *

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
  // #place(top + right,text(16pt, red)[#context globalvars.get()min])
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
    title: [Confined Quasi-2D Coulomb Systems: Theory, Algorithms, and Applications],
    // subtitle: [],
    author: text(23pt)[Xuanzhao Gao \ Supervisor: Prof. Zecheng Gan \ Co-supervisor: Prof. Jinguo Liu, Prof. Yang Xiang],
    institution: text(20pt)[Hong Kong University of Science and Technology (Guangzhou)],
    date: text(23pt)[2025-05-27],
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


= Background

== Quasi-2D charged systems

#timecounter(1)

Quasi-2D systems @mazars2011long are at the macroscopic scale in $x y$, but microscopic in $z$, so that are always modeled as doubly periodic in numerical simulations.
// Q2D systems are widely exist in nature and engineering, for example, cell membranes and electrolyte near surfaces.

#figure(
  image("figs/Q2D.png", width: 400pt),
  caption: [Illustration of a quasi-2D charged system.],
)

Coulomb interaction plays a key role in nature, leading to effect such as ion transportation and self-assembly. 

However, the Coulomb interaction decays as $r^(-1)$ in 3D, so that it is long ranged and singular at $r=0$, which make such simulation computationally expensive.


== Coulomb interaction

#timecounter(1)

Coulomb interaction plays a key role in nature, leading to effect such as ion transportation and self-assembly.

#figure(
  image("figs/self-assembly.png", width: 450pt),
  caption: [Self-assembly of nanoparticles via Coulomb interaction @luijten2014.
],)

However, the Coulomb interaction decays as $r^(-1)$ in 3D, so that is long ranged and singular at $r=0$, which make such simulation computationally expensive.
Complexity of a direct sum of the Coulomb interaction in doubly periodic systems is about $O(N^2 epsilon^(-1/3))$.

== Dielectric confinements

Polarizable surfaces also play a key role in nature, which arise naturally due to dielectric mismatch between different materials, leading to ion transport in nanochannels, pattern formation and self-assembly of colloidal and polymer monolayers.

#figure(
  image("figs/pattern.png", width: 350pt),
  caption: [Pattern formation in 2D dipolar systems due to dielectric mismatch @wang2019dielectric.
],)





== Algorithms for Q2D charged systems
#timecounter(1)

Methods have been developed to accelerate the Coulomb interaction in Q2D systems.

The very first method is the Ewald2D @parry1975electrostatic method based on the Ewald splitting of the Coulomb kernel. It is accurate but with $O(N^2)$ complexity.

To reduce the complexity, most methods rely on the following three strategies:

- *Fourier spectral method* (Lindbo & Tornberg, 2011; 2012; Nestler et al., 2015; Shamshirgar & Tornberg, 2017; Shamshirgar et al., 2021; Maxian et al., 2021): based on Ewald splitting and fast Fourier transform (FFT), with $O(N log N)$ complexity.

// @lindbo2011spectral @lindbo2012fast, @nestler2015fast, @shamshirgar2017spectral, @shamshirgar2021fast, @maxian2021fast

- *Fast multipole methods* (Greengard, 1987; Greengard & Rokhlin, 1987; Berman & Greengard, 1994; Yan & Shelley, 2018; Liang et al., 2020): accelerated by hierarchical low-rank compression, adaptive and with $O(N)$ complexity.

// @greengard1987fast, @yan2018flexibly, @liang2020harmonic, @liang2022hsma, @jiang2023jcp, @berman1994renormalization

- *Random batch Ewald* (Jin et al., 2021; Liang et al., 2022; Gan et al., 2024): based on Ewald splitting and random batch sampling, stochastic and with $O(N)$ complexity, efficient parallelization.

// @Jin2020SISC, @gan2024fast, @gan2024random, @liang2022superscalability

== Algorithms for Q2D charged systems
#timecounter(1)

For doubly periodic systems, one major challenge is the large prefactor in $O(N)$ or $O(N log N)$ compared to 3D-PBC solvers @mazars2011long, especially when the system is strongly confined in the $z$ direction, i.e., $H << L_x, L_y$.

- For the FFT based methods, *huge zero-padding* is required.

- For the FMM based methods, *more near field contributions* is needed.

// Some recently developed methods offer potential solutions to this challenge, including:
// - Anisotropic truncation kernel method @greengard2018anisotropic
// - Periodic FMM (Pei, Askham, Greengard & Jiang, 2023)
// // @jiang2023jcp
// - Dual-space multilevel kernel-splitting method @greengard2023dual

// However, these methods have not yet been extended to handle quasi-2D systems.

Another challenge is the polarization effect, various strategies have been developed in recent years, either by introducing image charges @yuan2021particle @liang2020harmonic or numerically solving the Poisson equation with interface conditions @nguyen2019incorporating @maxian2021fast @ma2021modified.
These methods are also accelerated by FFT or FMM to reach a complexity of $O(N log N)$ or $O(N)$.

However, the computational cost significantly increases compared to the homogeneous case, especially for strongly confined systems.




= Theoretical Analysis

== Ewald summations revisited

Ewald2D.

== Image charge method

In Q2D systems, dielectric permittivity is a function of $z$:

$
  epsilon (z) = cases(
  epsilon_u "if" z > H,
  epsilon_c "if" 0 < z < H,
  epsilon_d "if" z < 0,
)
$



= Random Batch Ewald2D Method

== The Algorithm

rbe2d

= Applications in MD Simulations

== All-atom simulations of SPC/E water

water

= Other Works

== Sum-of-Exponential Ewald2D method

SOEwald2D.

== Quasi Ewald method

QEM.

== Fast spectral sum-of-Gaussians method

FSSOG.

== Broken Symmetries in Q2D Systems

ssb

#pagebreak()

= Summary and Outlook

== Summary

During my PhD, I focused on confined quasi-2D charged systems.

I have developed a series of fast summation algorithms for quasi-2D charged systems, including:


== Publications
#timecounter(1)

#text(18pt)[
- *X. Gao*, Z. Gan, and Y. Li, Efficient particle-based simulations of Coulomb systems under dielectric nanoconfinement, (2025)
- *X. Gao*, X. Li, and J.-G. Liu, Programming guide for solving constraint satisfaction problems with tensor networks, #emph[Chinese Physics B] 34, 050201 (2025)
- *X. Gao*, Q. Zhou, Z. Gan, and J. Liang, Accurate error estimates and optimal parameter selection in Ewald summation for dielectrically confined Coulomb systems, Arxiv:2503.18126 (2025)
- Z. Gan, *X. Gao*, J. Liang, and Z. Xu, Random batch Ewald method for dielectrically confined Coulomb systems, Arxiv:2405.06333 (2025). Accepted by #emph[SIAM Journal on Scientific Computing]
- Z. Gan, *X. Gao*, J. Liang, and Z. Xu, Fast algorithm for quasi-2D Coulomb systems. #emph[Journal of Computational Physics] 113733, (2025)
- *X. Gao*, Y.-J. Wang, P. Zhang, and J.-G. Liu, Automated discovery of branching rules with optimal complexity for the maximum independent set problem, Arxiv:2412.07685 (2024)
- *X. Gao*, S. Jiang, J. Liang, Z. Xu, and Q. Zhou, A fast spectral sum-of-Gaussians method for electrostatic summation in quasi-2D systems, Arxiv:2412.04595 (2024)
- M. Roa-Villescas, *X. Gao*, S. Stuijk, H. Corporaal, and J.-G. Liu, Probabilistic inference in the era of tensor networks and differential programming, #emph[Physical Review Research] 6, 33261 (2024)
- *X. Gao* and Z. Gan, Broken symmetries in quasi-2D charged systems via negative dielectric confinement, #emph[The Journal of Chemical Physics] 161, (2024)
]

#pagebreak()

== Software packages

=== My packages

- *ChebParticleMesh.jl*#footnote(text(12pt)[#link("https://github.com/HPMolSim/ChebParticleMesh.jl")],): Toolkits for particle mesh methods (type-1 and type-2 NUFFT).
- *CuTropicalGEMM.jl*#footnote(text(12pt)[#link("https://github.com/TensorBFS/CuTropicalGEMM.jl")],): Custom GPU kernel for tropical matrix multiplication.
- *TreeWidthSolver.jl*#footnote(text(12pt)[#link("https://github.com/ArrogantGao/TreeWidthSolver.jl")],): Solving the treewidth problem (supported by GSoC 2024).
- *FastSpecSoG.jl*#footnote(text(12pt)[#link("https://github.com/HPMolSim/FastSpecSoG.jl")],): Implementation of the fast spectral SOG method.
- *EwaldSummations.jl*#footnote(text(12pt)[#link("https://github.com/HPMolSim/EwaldSummations.jl")],): Various Ewald summation methods with parallelization.
- *OptimalBranching.jl*#footnote(text(12pt)[#link("https://github.com/OptimalBranching/OptimalBranching.jl")],): Implementation of the optimal branching algorithm.

=== Contributions to popular packages

- *OMEinsum.jl*#footnote(text(12pt)[#link("https://github.com/under-Peter/OMEinsum.jl")],) (185 stars) and its backend *OMEinsumContractionOrders.jl*#footnote(text(12pt)[#link("https://github.com/TensorBFS/OMEinsumContractionOrders.jl")],): Optimizing the tensor network contraction order and contracting the tensor network.

#pagebreak()

// == Future Research Plans
// #timecounter(1)

// === Fast Summation Algorithms


// #pagebreak()

= Acknowledgements

#pagebreak()

#let figsize = 130pt

#align(center,
grid(columns: 5, gutter: 15pt,
  grid(rows : 3, image("photos/zechenggan.jpeg", width : figsize), "" ,  text[Prof. Zecheng Gan \ HKUST(GZ)]),
  h(30pt), 
  grid(rows : 3, image("photos/jinguoliu.jpeg", width : figsize), "" ,  text[Prof. Jinguo Liu, \ HKUST(GZ)]),
  h(30pt), 
  grid(rows : 3, image("photos/xiangyang.jpg", width : figsize), "" ,  text[Prof. Yang Xiang, \ HKUST]),
  grid(rows : 3, image("photos/zhenlixu.jpg", width : figsize), "" ,  text[Prof. Zhenli Xu \ SJTU]),
  h(30pt), 
  grid(rows : 3, image("photos/shidongjiang.jpeg", width : figsize), "" ,  text[Prof. Shidong Jiang \ Flatiron Institute]),
  h(30pt), 
  grid(rows : 3, image("photos/panzhang.jpeg", width : figsize), "" ,  text[Prof. Pan Zhang, \ ITP, CAS]),
  ),
)

Also thank Jiuyang Liang (SJTU & FI), Qi Zhou (SJTU) and Yijia Wang (ITP, CAS) for collaboration.


#focus-slide[
  Thank you for your attention!
]

#show: appendix

= Appendix <touying:unoutlined>

== Supplementary results

This is a result.

#pagebreak()

== References

#bibliography("ref.bib", style: "apa", title: none)
