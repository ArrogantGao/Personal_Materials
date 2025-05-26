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
  // #place(bottom + right,text(16pt, red)[#context globalvars.get()min])
]
#let clip(image, top: 0pt, bottom: 0pt, left: 0pt, right: 0pt) = {
  box(clip: true, image, inset: (top: -top, right: -right, left: -left, bottom: -bottom))
}

#let leftright(a, b) = align(center, grid(columns: 2, gutter: 30pt, align(left, box(width: 350pt)[#a]), align(left, box(width: 350pt)[#b])))
#let leftrightw1w2(a, b, w1, w2) = align(center, grid(columns: 2, gutter: 30pt, align(left, box(width: w1)[#a]), align(left, box(width: w2)[#b])))

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
    author: text(23pt)[GAO, Xuanzhao \ Supervisor: Prof. GAN, Zecheng \ Co-supervisor: Prof. LIU, Jinguo, Prof. XIANG, Yang],
    institution: text(20pt)[Advanced Materias Thrust, Function Hub \ Hong Kong University of Science and Technology (Guangzhou)],
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

#pagebreak()

== Self-introduction
#timecounter(1)

=== Education background

- *2021 - now*: PhD in Advanced Materials Thust, Function Hub, Hong Kong University of Science and Technology (Guangzhou). Majoring in Applied Mathematics.
- *2017 - 2021*: Bachelor in School of the Gifted Young, University of Science and Technology of China. Majoring in Physics and Computer Science.

=== Research interests

- Fast summation algorithms for quasi-2D charged systems.
- Tensor network methods for combinatorial optimization problems.

= Outline <touying:hidden>

#outline(title: none, indent: 1em, depth: 1)


= Background

== Molecular dynamics simulation
#timecounter(2)

Molecular dynamics simulation is a widely used method based on theoretical modeling and computer simulation, which plays a key role in many fields, including materials science, biophysics, and drug design.

#figure(
  grid(columns: 3, 
    image("figs/pnas.2102516118fig01.jpg", width: 450pt),
    h(50pt),
    image("figs/MD_workflow.png", width: 200pt),
  ),
  caption: [#text(15pt)[(Left) Study of the structure of SARS-CoV-2 replication-transcription complex by molecular dynamics simulation @brandon2021structural. (Right) Illustration of the workflow of molecular dynamics simulation.]],
)


== Confined quasi-2D Coulomb systems

#timecounter(2)

Quasi-2D systems @mazars2011long are at the macroscopic scale in $x y$, but microscopic in $z$, which are widely exist in nature and engineering, for example, cell membranes, electrolyte near surfaces and ultrathin polymer films.
They are always modeled as infinite planar layers, and are treated as doubly periodic in numerical simulations.

#figure(
  image("figs/Q2D.png", width: 400pt),
  caption: [#text(15pt)[Illustration of a quasi-2D charged system.]],
)

We consider Q2D systems with charged particles, which are modeled as point charges.
These charges interact with each other via Coulomb interaction, which plays a key role in nature, leading to effect such as self-assembly. 

// However, the Coulomb interaction decays as $r^(-1)$ in 3D, so that it is long ranged and singular at $r=0$, which make such simulation computationally expensive.


// // == Coulomb interaction
// #pagebreak()

// #timecounter(1)

// Coulomb interaction plays a key role in nature, leading to effect such as ion transportation and self-assembly.

// #figure(
//   image("figs/self-assembly.png", width: 450pt),
//   caption: [#text(15pt)[Self-assembly of nanoparticles via Coulomb interaction @luijten2014.]],
// )


// == Dielectric confinements
#pagebreak()

#timecounter(1)

Confinements means the particles are confined by substrates below and above.

Polarization arises naturally due to dielectric mismatch between different materials, leading to effects such as modulation of ion transport @antila2018dielectric, and pattern formation @wang2019dielectric.

#leftright(
  figure(
    image("figs/ion_transport.png", width: 350pt),
    caption: [#text(15pt)[Dielectric modulation of ion transport near interfaces @antila2018dielectric.]],
  ),
  figure(
    image("figs/pattern.png", width: 350pt),
    caption: [#text(15pt)[Pattern formation in 2D dipolar systems due to dielectric mismatch @wang2019dielectric.]],
  ),
)

== Algorithms for Q2D charged systems
#timecounter(1)

The bottleneck of simulating Q2D system is the Coulomb interaction.
The electrostatic interaction in Q2D systems is given as follows:
$
  U = 1 / 2 sum_(i,j=1)^N sum_(bold(m)) q_i q_j G(bold(r)_i, bold(r)_j + L_bold(m)),
$
where $G(bold(r), bold(r)')$ is the Green's function, $L_bold(m) = (L_x m_x, L_y m_y, 0)$.

If the system is homogeneous, the Green's function is given by:
$
  G(bold(r), bold(r)') = 1 / abs(bold(r) - bold(r)')
$
It decays slowly and is singular at $r=0$, which make such summation computationally expensive.

#pagebreak()

// Coulomb interaction 

Methods have been developed to accelerate the calculation of Coulomb interaction in Q2D systems.

The very first method is the Ewald2D @parry1975electrostatic method based on the Ewald splitting of the Coulomb kernel. It is accurate but with *$O(N^2)$* complexity.

To reduce the complexity, most methods rely on the following two strategies:

- *Fourier spectral method* @maxian2021fast @yuan2021particle: based on Ewald splitting and fast Fourier transform (FFT), with *$O(N log N)$* complexity. Example: PPPM2D, ICM-PPPM.

// Lindbo & Tornberg, 2011; 2012; Nestler et al., 2015; Shamshirgar & Tornberg, 2017; Shamshirgar et al., 2021; Maxian et al., 2021
// @lindbo2011spectral @lindbo2012fast, @nestler2015fast, @shamshirgar2017spectral, @shamshirgar2021fast, @maxian2021fast

- *Fast multipole methods* @greengard1987fast @liang2020harmonic: accelerated by hierarchical low-rank compression, adaptive and with *$O(N)$* complexity. Example: HSMA, periodic FMM.

// (Greengard, 1987; Greengard & Rokhlin, 1987; Berman & Greengard, 1994; Yan & Shelley, 2018; Liang et al., 2020)
// @greengard1987fast, @yan2018flexibly, @liang2020harmonic, @liang2022hsma, @jiang2023jcp, @berman1994renormalization

// - *Random batch Ewald* @Jin2020SISC: based on Ewald splitting and random batch sampling, stochastic and with $O(N)$ complexity, efficient parallelization.

// @Jin2020SISC, @gan2024fast, @gan2024random, @liang2022superscalability

== Algorithms for Q2D charged systems
#timecounter(1)

For doubly periodic systems, there are two major challenges:

=== Strongly confinements

Confinements leads to large prefactors compared to 3D-PBC solvers @mazars2011long, especially when $H << L_x, L_y$:
- For the FFT based methods, *huge zero-padding* is required.
- For the FMM based methods, *more near field contributions* is needed.

=== Polarization effect

Various strategies have been developed in recent years, by:
- Introducing image charges @yuan2021particle @liang2020harmonic.
- Numerically solving the Poisson's equation with interface conditions @nguyen2019incorporating @maxian2021fast @ma2021modified.
Both of them significantly increase the computational cost compared to the homogeneous case.

// Main target of our work is to develop efficient methods for quasi-2D systems, overcomes the challenges mentioned above.

= Theoretical Analysis

== Ewald summations revisited

#timecounter(1)

Key point: handle short-range and long-range contributions separately.

$
  1 / r = underbrace((1 - "erf"(alpha r)) / r, "short-range") + underbrace("erf"(alpha r) / r, "long-range") = "erfc"(alpha r) / r + "erf"(alpha r) / r
$

#align(center,
  image("figs/erf.png", width: 400pt),
)

These kernels decay rapidly in real space and Fourier space, respectively.

#pagebreak()

#timecounter(1)

The short-range part is summed in the real space, and the long-range part is summed in the Fourier space.
The total interaction energy can be computed as follows:
$
  U_s = 1 / 2 sum_(i,j=1)^N sum_(bold(m)) q_i q_j "erfc"(alpha |bold(r)_i - bold(r)_j + L_bold(m)|) / abs(bold(r)_i - bold(r)_j + L_bold(m))
$
$
  U_l = pi / (2 L_x L_y) sum_(i,j=1)^N q_i q_j sum_(bold(h) != bold(0)) e^(i bold(h) dot (bold(r)_i - bold(r)_j)) / h cal(G)_alpha (bold(h), z_i - z_j) - alpha / sqrt(pi) sum_(i=1)^N q_i^2 + cal(J)_0
$
where $bold(m) = (m_x, m_y)$ and $bold(h) = ((2 pi m_x) / L_x, (2 pi m_y) / L_y)$ is the reciprocal lattice vector.
$
  cal(J)_0 = -pi / (L_x L_y) sum_(i,j=1)^N q_i q_j ("erf"(alpha z) + alpha / sqrt(pi) e^(- alpha^2 z^2)) \
  cal(G)_alpha (bold(h), z) = xi^+ (bold(h), z) + xi^- (bold(h), z) = e^(-h z) "erfc"(h / (2 alpha) + alpha z) + e^(h z) "erfc"(h / (2 alpha) - alpha z)
$
where $z = |z_i - z_j|$.
The summaitons can be truncated as $r_c = s / alpha$, and $k_c = 2 s alpha$.

The resulting method is the so-called Ewald2D method, have a complexity of $O(N^2)$ due to the *unseparable double summation* over $i$ and $j$ in $U_l$.

#pagebreak()

#timecounter(2)

One way to accelerate the summation is transforming the system as triply periodic one by adding vacuum layer in the $z$ direction.

#align(center,
  image("figs/vacuum_layer.png", width: 200pt),
)

The summation is reformulated as:
$
  U_l = underbrace((2 pi) / (L_x L_y L_z) sum_(bold(k) != bold(0)) e^(- k^2 / (4 alpha^2)) / k^2 (sum_(i = 1)^N q_i e^(i bold(k) dot bold(r)_i))^2 - alpha / sqrt(pi) sum_(i=1)^N q_i^2, "Ewald3D summation") + underbrace(U_("YB") , #text[$O(N)$]) + underbrace(U_("ELC"), #text[$O(N^2)$]) + underbrace(U_("Trap") , "Error")
$
where $bold(k) = 2 pi (m_x / L_x, m_y / L_y, m_z / L_z)$ is the reciprocal lattice vector.
Complexity of EwaldELC method is *$O(N^(1.5))$*.
// and can easily be accelerated by FFT or FMM to reach linear complexity.


== Image charge method revisited

#timecounter(1)

In dielectrically confined Q2D systems, we assume that the dielectric permittivity is uniform in each region, then can be written as a function of $z$:
$
  epsilon (z) = cases(
  epsilon_u "if" z > H,
  epsilon_c "if" 0 < z < H,
  epsilon_d "if" z < 0,
)
$

The governing equation of the Green's function $G(bold(r), bold(r)')$:
$
  cases(
    gradient_bold(r) (epsilon(bold(r)) gradient_bold(r) G(bold(r), bold(r)')) = - delta(bold(r) - bold(r)') & "if" bold(r) in R^3,
    G(bold(r), bold(r)') |_- = G(bold(r), bold(r)') |_+ & "if" bold(r) in partial Omega_c,
    epsilon_c gradient_bold(r) G(bold(r), bold(r)') |_- = epsilon_u gradient_bold(r) G(bold(r), bold(r)') |_+ & "if" bold(r) in partial Omega_c \u{2229} partial Omega_u,
    epsilon_c gradient_bold(r) G(bold(r), bold(r)') |_+ = epsilon_d gradient_bold(r) G(bold(r), bold(r)') |_- & "if" bold(r) in partial Omega_c \u{2229} partial Omega_d,
    G(bold(r), bold(r)') |_- = G(bold(r), bold(r)') |_+ & "if" bold(r) in partial Omega_u,
  )
$

#leftright(
  text[
    Define the dielectric reflection rate:
    $
      gamma_u = (epsilon_c - epsilon_u) / (epsilon_c + epsilon_u), " and "
      gamma_d = (epsilon_c - epsilon_d) / (epsilon_c + epsilon_d),
    $
  ],
  align(center,
    image("figs/gamma.png", width:100pt),
  ),
)


#pagebreak()

#timecounter(1)

We get the so-called *image charge series* for the polarization potential:
$
  G(bold(r), bold(r)') = 1 / (4 pi epsilon_c) [1 / abs(bold(r) - bold(r)') + sum_(l=1)^infinity (gamma_+^(l) / abs(bold(r) - bold(r)_(+)^(l) ') + gamma_-^(l) / abs(bold(r) - bold(r)_(-)^(l) '))]
$
where $gamma_+^(l) = gamma_d^(ceil(l/2)) gamma_u^(floor(l/2))$ and $gamma_-^(l) = gamma_d^(floor(l/2)) gamma_u^(ceil(l/2))$, and $bold(r)_(+)^(l) = (x, y, z_(+)^(l))$ and $bold(r)_(-)^(l) = (x, y, z_(-)^(l))$ are the positions of the $l$-th level image charges, and 
$
  z_+^(l) = (-1)^l z + 2 ceil(l/2) H,
  z_-^(l) = (-1)^l z - 2 floor(l/2) H.
$
Assume that $epsilon_u, epsilon_d, epsilon_c$ are all positive, then $gamma_u, gamma_d in (-1, 1)$, the series is decaying exponentially, and can be truncated at the $M$-th level.

#align(center,
  image("figs/icm_system.png", width: 400pt)
)


$N$ particles inhomogenous system $arrow$ $(2M + 1) N$ particles homogenous system

== ICM-EwaldELC method

#timecounter(1)

A combination of image charge method and EwaldELC method is called ICM-EwaldELC method, which is a standard technique for dielectrically confined systems.

- *$L_z$*: thickness of the padded 3D system.

- *$M$*: number of layers of image charges.

We noticed that the total error is quite complex for the ICM-EwaldELC method:

#figure(
  image("figs/error_yuan.png", width: 500pt),
  caption: [#text(15pt)[Relative precision in the electrostatic energy computed by ICM-PPPM as a function of the number of reciprocal reflections and the slab factor for a 2:1 salt confined between dielectric interfaces @yuan2021particle.]],
)

A rigorous error analysis is thus very important for the effective application.

== Error estimation

#timecounter(1)

In our work @icm_error, we show that the error of ICM-EwaldELC method is given as follows:
$
  cal(E) ~  & underbrace(O(e^(-s^2) / (s^2)), "Ewald sum") + underbrace(O(e^( - alpha^2 (L_z - H)^2 ) ), "Trapezoidal rule") + underbrace(O(abs(gamma_u gamma_d)^(floor((M + 1) / 2)) e^( - (4 pi H floor((M + 1) / 2)) / max(L_x, L_y) )), "Image charge") \ 
  & + underbrace(O(e^( - (2 pi (L_z - H)) / max(L_x, L_y) ) + sum_(l=1)^M (abs(gamma_d^(ceil(l/2)) gamma_u^(floor(l/2))) + abs(gamma_d^(floor(l/2)) gamma_u^(ceil(l/2)))) e^( - (2 pi (L_z - (l + 1) H)) / max(L_x, L_y) )), "ELC term")
$
with the assumption that $abs(gamma_u) <= 1$ and $abs(gamma_d) <= 1$.

#pagebreak()

=== Image charge series

#timecounter(1)

Error induced by truncation of the image charge series decays exponentially with the layer number $M$, the speed depends on $gamma$ and aspect ratio $L_x / H$.
$
  cal(E) ~ O(abs(gamma_u gamma_d)^(floor((M + 1) / 2)) e^( - (4 pi H floor((M + 1) / 2)) / max(L_x, L_y) ))
$

#figure(
  image("figs/icm_image_charge_error.png", width: 550pt),
  caption: [#text(15pt)[Error of image charge series as a function of the layer number $M$ for different $gamma$ and aspect ratio $L_x / H$, the dashed lines are fitted according to the theoretical prediction.]],
)

#pagebreak()

=== ELC term

#timecounter(1)

First consider the case when $gamma_u = gamma_d = 0$, i.e., the system is homogeneous, then
$
  cal(E) ~ O(e^( - (2 pi (L_z - H)) / max(L_x, L_y)))
$
we thus define the padding ratio $R = (L_z - H) / max(L_x, L_y)$, and ELC term decays exponentially with $R$.

#figure(
  image("figs/elc_error_force.svg", width: 300pt),
  caption: [#text(15pt)[Error of ELC term as a function of the padding ratio $R$ for different aspect ratio.]],
)

#pagebreak()

#timecounter(1)

When $gamma_u, gamma_d != 0$, there are two cases, depending on the value of $abs(gamma_u gamma_d e^((4 pi H) / max(L_x, L_y)))$:

If $abs(gamma_u gamma_d e^((4 pi H) / max(L_x, L_y))) <= 1$, the ELC term is independent of $M$:
$
  cal(E) ~ O((gamma_u e^((2 pi H) / max(L_x, L_y)) + gamma_d e^((2 pi H) / max(L_x, L_y)) + 2)e^( - (2 pi (L_z - H)) / max(L_x, L_y)))
$

#figure(
  image("figs/error_icm_pad_gamma_0.6.svg", width: 600pt),
  caption: [#text(15pt)[Error of ELC term as a function of $M$ and $R$. We fix $L_x = L_y = 10$, $H = 0.5$ and $gamma_u = gamma_d = 0.6$ so that $abs(gamma_u gamma_d e^((4 pi H) / max(L_x, L_y))) < 1$]],
)

#pagebreak()

#timecounter(1)

Interestingly, if $abs(gamma_u gamma_d e^((4 pi H) / max(L_x, L_y))) > 1$, the ELC term also decays exponentially to $R$, but grows exponentially to $M$:
$
  cal(E) ~ O(abs(gamma_u gamma_d e^((4 pi H) / max(L_x, L_y)))^(M / 2) e^(- (2 pi (L_z - H)) / max(L_x, L_y)))
$

#figure(
  image("figs/error_icm_pad_gamma_1.svg", width: 600pt),
  caption: [#text(15pt)[Error of ELC term as a function of $M$ and $R$. We fix $L_x = L_y = 10$, $H = 0.5$ and $gamma_u = gamma_d = 1$ so that $abs(gamma_u gamma_d e^((4 pi H) / max(L_x, L_y))) > 1$]],
)

// #pagebreak()

// Our results well explain the behavior of the error of ICM-EwaldELC method in the literature @yuan2021particle.

// #figure(
//   image("figs/icm_error_yuan.png", width: 700pt),
//   caption: [#text(15pt)[Repeat of the results in @yuan2021particle with fitted lines according to our theoretical prediction, where $L_x = L_y = 15$, $H = 5$. The system is a 2:1 electrolyte with $39$ particles.]],
// )

== Parameter selection

#timecounter(2)

Given the system size $L_x, L_y, H$, the dielectric constants $epsilon_u, epsilon_d, epsilon_c$, and desired precision $epsilon$, we can select the parameters as follows:

*Step 1*:
Select the layer number $M$:
$
  M ~ (2 log(epsilon) - (4 pi H) / max(L_x, L_y) - log(abs(gamma_u gamma_d))) / (log(abs(gamma_u gamma_d)) - (4 pi H) / max(L_x, L_y))
$

*Step 2*:
Select height of the padded 3D system $L_z$:
- Case 1: $abs(gamma_u gamma_d e^((4 pi H) / max(L_x, L_y))) < 1$
$
  L_z >= H + (max(L_x, L_y) / (2 pi)) (log(1 / epsilon) + log(abs(gamma_u + gamma_d + e^(-2 pi H / max(L_x, L_y)))))
$
- Case 2: $abs(gamma_u gamma_d e^((4 pi H) / max(L_x, L_y))) >= 1$
$
  L_z >= (M + 1) H + (max(L_x, L_y) / (2 pi)) (log(1 / epsilon) + log(abs(gamma_u gamma_d)))
$ 

*Step 3*:
Select $alpha$ and $s$ for Ewald splitting: $alpha >= sqrt(- log(epsilon)) / (L_z - H)$, $e^(-s^2) / s^2 <= epsilon$

#pagebreak()

#timecounter(1)

We applied our method for system with $L_x = L_y = 10$, $H = 1$, $gamma_u = gamma_d = 0.6$ and $1$.

#align(center,
  grid(columns:2, gutter: 15pt, 
    image("figs/icm_parameter_table.png", width: 200pt),
    figure(
      image("figs/icm_parameter_benchmark.png", width: 450pt),
      caption: [#text(15pt)[Numerical tests for the parameter selection. Dashed lines are numerically error, dots are pre-selected parameters.]],
    ),
  )
)

= Algorithms for Q2D Systems

== Random batch Ewald2D method

#timecounter(2)

=== Overview

- For the dielectrically confined Q2D systems.
- Accelerated by random batch Ewald method, with linear complexity.
- Reduce the complexity of ICM-EwaldELC based method from $O(M N)$ to $O(M + N)$.
- Parallelizable and scalable.

=== Basic idea

Recall that in EwaldELC, long-range energy of the Q2D system is approximated by that of a padded 3D system.
$
  U_("Q2D")^l approx underbrace((2 pi) / (L_x L_y L_z) sum_(bold(k) != bold(0)) e^(- k^2 / (4 alpha^2)) / k^2 (sum_(i = 1)^N q_i e^(i bold(k) dot bold(r)_i))^2 - alpha / sqrt(pi) sum_(i=1)^N q_i^2, #text[$U_("3D")^l$]) + U_("YB")
$
// A direct forward idea is to calculate $U_("3D")^l$ via random batch Ewald, resulting in an $O(N)$ algorithm.
By accelerating $U_("3D")^l$, the complexity of EwaldELC is reduced. For example, PPPM2D accelerates via FFT and has a complexity of $O(N log N)$. However, it suffers from the huge zero padding for strongly confined systems and the communication overhead of FFT when using multiple CPUs.

Our idea is to accelerate $U_("3D")^l$ via random batch Ewald method.

#pagebreak()

#timecounter(1)

Random batch Ewald method @Jin2020SISC is a stochastic method proposed to accelerate the Ewald3D summation.
The main idea is to sample the long-range part in the Fourier space using importance sampling instead of summing over all the Fourier modes.
$
  sum_(bold(k) != bold(0)) e^(- k^2 / (4 alpha^2)) / k^2 (sum_(i = 1)^N q_i e^(i k r_i))^2 approx (sum_(bold(k) != bold(0)) e^(- k^2 / (4 alpha^2))) / P sum_(k in cal(K)_P) 1 / k^2 (sum_(i = 1)^N q_i e^(i k r_i))^2,
$
where $cal(K)_P$ is the set of $P$ Fourier modes sampled from the distribution $cal(P) (k) ~ e^(- k^2 / (4 alpha^2))$.
$P$ is proven to be an $O(1)$ constant so that the total complexity is $O(P N)$.

#figure(
  grid(columns:2, gutter: 25pt, 
    image("figs/rbe_accuracy.png", height: 220pt),
    image("figs/rbe_scale.png", height: 220pt),
  ),
  caption: [#text(15pt)[Accuracy and scalability of random batch Ewald method @liang2022superscalability.]],
)

#pagebreak()

#timecounter(1)

RBE is not sensitive to the aspect ratio of the system.

#figure(
  image("figs/rbe2d_accuracy.png", width: 600pt),
  caption: [#text(15pt)[Accuracy of RBE2D method for simulating electrolytes.]],
)

#pagebreak()

#timecounter(1)

=== RBE2D for dielectrically confined systems

Combined with the image charge method, RBE2D can be applied to the dielectrically confined systems:
$
  U_("3D")^M = (2 pi) / (L_x L_y L_z) sum_(bold(k) != bold(0)) e^(- k^2 / (4 alpha^2)) / k^2 rho_k overline(rho)_k^M - alpha / sqrt(pi) sum_(i=1)^N q_i^2
$
where
$
  rho_k = sum_(i = 1)^N q_i e^( - i k r_i) \ overline(rho)_k^M = sum_(j = 1)^N q_j (e^(i k r_j) + sum_(l = 1)^M (gamma_+^l e^(i k r_(j +)^l) + gamma_-^l e^(i k r_(j -)^l)))
$

However, cost of evaluating $overline(rho)_k^M$ is $O(M N)$, which is too expensive.

#pagebreak()

#timecounter(1)

=== From $O(M N)$ to $O(M + N)$

Recall the location of image charges:
$
  r_(+)^(l) = (x, y, (-1)^l z + 2 ceil(l/2) H) " and " r_(-)^(l) = (x, y, (-1)^l z - 2 floor(l/2) H)
$
Notice that
$
  overline(rho)_k^M = sum_(j = 1)^N q_j e^(i bold(k) dot bold(r)_j) (1 + e^(-2 i k_z z_j) Y_("odd")(k_z) + Y_("even")(k_z)) \
  Y_("odd")(k_z) = sum_(l = 1, "odd")^M (gamma_+^l e^(i k_z (l + 1) H) + gamma_-^l e^( - i k_z (l - 1) H)) \
  Y_("even")(k_z) = sum_(l = 1, "even")^M (gamma_+^l e^(i k_z l H) + gamma_-^l e^( - i k_z (l - 1) H))
$
Then for each Fourier mode $k$, one can compute $Y_("odd")(k_z)$ and $Y_("even")(k_z)$ in $O(M)$ time, and then $overline(rho)_k^M$ in $O(N)$ time. The overall complexity is $O(M + N)$.


#pagebreak()

#timecounter(1)

#align(center,
  image("figs/rbe2d_table.png", height: 110pt),
)

#figure(
  image("figs/rbe2d_dielectric.png", width: 500pt),
  caption: [#text(15pt)[Accuracy of RBE2D method for dielectrically confined electrolytes, $P = 100$.]],
)

#pagebreak()

=== CPU performance

#timecounter(1)

The SPC/E bulk water systems are simulated, where the system dimensions are set as $L_x = L_y = H$, and the system size changes as one varies $N$ with the density of water being fixed at $1 "g/cm"^3$.

#figure(
  image("figs/rbe2d_runtime.png", width: 700pt),
  caption: [#text(15pt)[CPU performance of RBE2D method, compared with PPPM. (Left) number of CPU cores is fixed as 64, (Right) number of atoms is fixed as 139968, inset shows the strong parallel efficiency.]],
)

#pagebreak()

#timecounter(1)

=== Advantages

- Mesh free, only use simple data structure.
- Not sensitive to the aspect ratio of the system.
- Reduce the ICM complexity from $O(M N)$ to $O(M + N)$.
- Parallelizable and scalable.

=== Disadvantages

- Stochastic, replying on law of large numbers.
- Limited to planar systems with uniform dielectric constant.
- Strong convergence and ergodicity remaining an open question.

== Other methods

#timecounter(2)

=== Sum-of-Exponential Ewald2D method

- For the Q2D systems in homogeneous media.

- Reduce the complexity of double summation to $O(N)$ via SOE approximation.

- Accelerated by random batch Ewald method, with linear complexity.

#linebreak()

=== Quasi-Ewald method

- For the negatively confined Q2D systems.

- Remove the divergence due to negative dielectric constant via singularity subtractions.


// == Short summary

// #timecounter(1)

// In summary, we developed a series of methods for simulating quasi-2D charged systems, including:

// - SOEwald2D method for non-dielectrically confined systems.
// - RBE2D method for dielectrically confined systems.
// - QEM method for negatively confined systems.


= Applications in MD Simulations

== All-atom simulations of SPC/E water

#timecounter(1)

=== SPC/E water in homogeneous media

RBE2D is applied to the all-atom simulations of SPC/E water in homogeneous media, where $L_x = L_y = H = 55.9 angstrom$ and consists of 17496 atoms.

#figure(
  image("figs/spce.png", width: 550pt),
  caption: [#text(15pt)[All-atom simulations of SPC/E water in homogeneous media.]],
)

#pagebreak()

=== SPC/E water in dielectrically confined systems

#timecounter(1)

We applied our method to the simulations of SPC/E water confined by slabs with different dielectric constants, where $gamma_u = gamma_d = -0.5$, the system consists of 53367 atoms.

#figure(
  image("figs/spce_dielectric.png", width: 500pt),
  caption: [#text(15pt)[All-atom simulations of SPC/E water under dielectric confinements.]],
)

== Simulations of electrolytes

=== Electrolytes in homogeneous media

We applied the SOEwald2D method for the simulations of electrolytes in homogeneous media, with external electric field in $z$ direction.

#figure(
  image("figs/nonconfined.png", width: 600pt),
  caption: [#text(15pt)[Results of SOEWald2D method for simulations of electrolytes in homogeneous media.]],
)

#pagebreak()

=== Electrolytes in dielectrically confined systems

We studied the electrolytes confined by slabs with different dielectric constants.

#figure(
  image("figs/confined_sym.png", width: 650pt),
  caption: [#text(15pt)[Concentration of ions in dielectrically confined 3:1 electrolyte systems with symmetric dielectric interfaces, (a) $gamma = -0.95$, (b) $gamma = 0.95$.]],
)

#pagebreak()

=== Electrostatic interactions under negative dielectric confinements

With the quasi-Ewald method, we examined the force between two charged particles in the presence of negative dielectric confinements.
The force is oscillatory as a function of the distance between the two particles, and like-charge attraction is observed.

#figure(
  image("figs/qem_force.png", width: 700pt),
  caption: [#text(15pt)[Force between two charged particles in the presence of negative dielectric confinements.]],
)

#pagebreak()

=== Spontaneous symmetry broken via negative dielectric confinements

By simulating 1:1 electrolyte systems under negative dielectric confinements, we observed spontaneous symmetry broken solely via dielectric confinements.
The charged particles gather into checkerboard patterns, and the system exhibits a broken symmetry.

#figure(
  image("figs/ssb.png", width: 700pt),
  caption: [#text(15pt)[Spontaneous symmetry broken via dielectric confinements.]],
)

// #pagebreak()

// #figure(
//   image("figs/ssb_force.png", width: 600pt),
//   caption: [#text(15pt)[Detailed force analysis of the symmetry broken phase.]],
// )

= Summary and Outlook

== Summary and outlook

#timecounter(1)

During my PhD, I focused on confined quasi-2D charged systems, including:

- Theoretical analysis of Ewald summation for dielectric-confined planar systems.

- Development of novel fast algorithms for simulating quasi-2D Coulomb systems, for non-dielectrically confined systems, dielectrically confined systems, and negatively confined systems, which overcomes the difficulty of strongly confinement and polarization effects.

- Applications in MD simulations, including all-atom simulations of SPC/E water, and the observation of spontaneous symmetry broken solely via dielectric confinements for the first time.


=== Future works
#timecounter(1)

- Extend our methods to more complex systems, including systems with complex geometries such as curved surfaces and polarizable particles.

- Extend our methods to different interaction kernels, such as the dipolar interactions, the Stokes interactions and the Yukawa interactions.

- Application to the simulation of real world systems, such as biomolecular systems and colloidal systems.

== List of publications

#text(16pt)[
1. *X. Gao*, Z. Gan, and Y. Li, Efficient particle-based simulations of Coulomb systems under dielectric nanoconfinement, (2025), under preparation

2. *X. Gao*, Q. Zhou, Z. Gan, and J. Liang, Accurate error estimates and optimal parameter selection in Ewald summation for dielectrically confined Coulomb systems, Arxiv:2503.18126 (2025). Accepted by #emph[#text(blue, "Journal of Chemical Theory and Computation")]

3. *X. Gao*, X. Li, and J.-G. Liu, Programming guide for solving constraint satisfaction problems with tensor networks, #emph[#text(blue, "Chinese Physics B")] 34, 050201 (2025)

4. (Alphabetical order) Z. Gan, *X. Gao*, J. Liang, and Z. Xu, Random batch Ewald method for dielectrically confined Coulomb systems, Arxiv:2405.06333 (2025). Accepted by #emph[#text(blue, "SIAM Journal on Scientific Computing")]

5. (Alphabetical order) Z. Gan, *X. Gao*, J. Liang, and Z. Xu, Fast algorithm for quasi-2D Coulomb systems. #emph[#text(blue, "Journal of Computational Physics")] 113733, (2025)

6. *X. Gao*, Y.-J. Wang, P. Zhang, and J.-G. Liu, Automated discovery of branching rules with optimal complexity for the maximum independent set problem, Arxiv:2412.07685 (2024)

7. (Alphabetical order) *X. Gao*, S. Jiang, J. Liang, Z. Xu, and Q. Zhou, A fast spectral sum-of-Gaussians method for electrostatic summation in quasi-2D systems, Arxiv:2412.04595 (2024)

8. M. Roa-Villescas, *X. Gao*, S. Stuijk, H. Corporaal, and J.-G. Liu, Probabilistic inference in the era of tensor networks and differential programming, #emph[#text(blue, "Physical Review Research")] 6, 33261 (2024)

9. *X. Gao* and Z. Gan, Broken symmetries in quasi-2D charged systems via negative dielectric confinement, #emph[#text(blue, "The Journal of Chemical Physics")] 161, (2024)
]

#pagebreak()

== Software packages

#timecounter(1)

- *ChebParticleMesh.jl*#footnote(text(12pt)[#link("https://github.com/HPMolSim/ChebParticleMesh.jl")],): Toolkits for particle mesh methods (type-1 and type-2 NUFFT).
- *CuTropicalGEMM.jl*#footnote(text(12pt)[#link("https://github.com/TensorBFS/CuTropicalGEMM.jl")],): Custom GPU kernel for tropical matrix multiplication.
- *TreeWidthSolver.jl*#footnote(text(12pt)[#link("https://github.com/ArrogantGao/TreeWidthSolver.jl")],): Solving the treewidth problem (supported by GSoC 2024).
- *FastSpecSoG.jl*#footnote(text(12pt)[#link("https://github.com/HPMolSim/FastSpecSoG.jl")],): Implementation of the fast spectral SOG method.
- *EwaldSummations.jl*#footnote(text(12pt)[#link("https://github.com/HPMolSim/EwaldSummations.jl")],): Various Ewald summation methods with parallelization.
- *OptimalBranching.jl*#footnote(text(12pt)[#link("https://github.com/OptimalBranching/OptimalBranching.jl")],): Implementation of the optimal branching algorithm.
- *ExTinyMD.jl*#footnote(text(12pt)[#link("https://github.com/HPMolSim/ExTinyMD.jl")],): A simple MD simulation package.
- *SoEwald2D.jl*#footnote(text(12pt)[#link("https://github.com/HPMolSim/SoEwald2D.jl")],): Implementation of the sum-of-exponential Ewald2D method.
- *QuasiEwald.jl*#footnote(text(12pt)[#link("https://github.com/HPMolSim/QuasiEwald.jl")],): Implementation of the quasi-Ewald method.

// === Contributions to popular packages

// - *OMEinsum.jl*#footnote(text(12pt)[#link("https://github.com/under-Peter/OMEinsum.jl")],) (185 stars) and its backend *OMEinsumContractionOrders.jl*#footnote(text(12pt)[#link("https://github.com/TensorBFS/OMEinsumContractionOrders.jl")],): Optimizing the tensor network contraction order and contracting the tensor network.

// = Acknowledgements

#pagebreak()

== Acknowledgements

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

== Performance of RBE2D method

#pagebreak()

#figure(
  grid(rows:2, gutter: 10pt,
    image("figs/spce_time.png", width: 700pt),
    image("figs/spce_time_dielectric.png", width: 700pt)
  ),
  caption: [#text(15pt)[CPU performance of RBE2D method, compared with PPPM, for SPC/E water without/with dielectric confinement, with different number of cores.]],
)


== Sum-of-Exponential Ewald2D method

=== Sum-of-Exponential approximation

Approximating an arbitrary function $f(x)$ by a sum of exponential functions:
$
  abs(f(x) - sum_(l=1)^M w_l e^(- s_l abs(x))) < epsilon
$
For Gaussian kernel, VPMR method @AAMM-13-1126 is shown to be an optimal choice:

#figure(
  image("figs/soe_vpmr.png", width: 400pt),
  caption: [#text(15pt)[Approximation of Gaussian kernel by VPMR method.]],
)

#pagebreak()

=== SOE + Ewald2D

Main problem of Ewald2D method: double summation, resulting in $O(N^2)$ complexity.
$
  sum_(i, j) q_i q_j e^(i h rho_(i j)) (e^(-h z) "erfc"(h / (2 alpha) + alpha z) + e^(h z) "erfc"(h / (2 alpha) - alpha z))
$
where $z = abs(z_i - z_j)$. Notice that
$
  xi^+ = e^( - h z) "erfc"(h / (2 alpha) + alpha z) =  (2 alpha) / sqrt(pi) e^(-h^2 / (4 alpha^2)) e^(-h z) integral_(-z)^infinity e^(- alpha^2 t^2) e^(-h t) d t\
  xi^- = e^(h z) "erfc"(h / (2 alpha) - alpha z) =  (2 alpha) / sqrt(pi) e^(-h^2 / (4 alpha^2)) e^(h z) integral_z^infinity e^(- alpha^2 t^2) e^(-h t) d t
  // e^(-alpha^2 t^2) approx sum_(l=1)^M w_l e^(- s_l alpha abs(t))
$
Then the integral can be calculated analytically:
$
  xi^+ approx (2 alpha) / sqrt(pi) e^(-h^2 / (4 alpha^2)) e^(-h z) integral_(-z)^infinity sum_(l=1)^M w_l e^(- s_l alpha abs(t)) e^(-h t) d t = (2 alpha) / sqrt(pi) e^(-h^2 / (4 alpha^2)) sum_(l=1)^M w_l e^(- s_l alpha z) / (alpha s_l + h) \
  xi^- approx (2 alpha) / sqrt(pi) e^(-h^2 / (4 alpha^2)) e^(h z) integral_z^infinity sum_(l=1)^M w_l e^(- s_l alpha abs(t)) e^(-h t) d t = (2 alpha) / sqrt(pi) e^(-h^2 / (4 alpha^2)) sum_(l=1)^M w_l (-e^(- s_l alpha z) / (alpha s_l - h) + (2 alpha s_l e^(-h z)) / ((alpha s_l)^2 - h^2))
$

#pagebreak()

#figure(
  image("figs/soe_xi.png", width: 600pt),
  caption: [#text(15pt)[Absolute error of approximating $xi^+$ and $xi^-$ as a function of $h$ and $z$.]],
)

#pagebreak()

=== Summing up exponentials

Now we have the following double summation:
$
  S = sum_(i, j) q_i q_j e^(-h abs(z_i - z_j))
$
Direct sum still leads to $O(N^2)$ complexity.
To reduce the complexity, we first sort the particles in $z$ direction, so that
$
  0 < z_1 < z_2 < ... < z_N < H
$
Then we can rewrite the summation as follows:
$
  S = 2 sum_(i=1)^N sum_(j=i + 1)^N q_i q_j e^(-h (z_j - z_i)) + sum_(i=1)^N q_i^2 = 2 sum_(i=1)^N q_i e^(h z_i) underbrace(sum_(j=i + 1)^N q_j e^(-h z_j), #text[$A(i)$]) + sum_(i=1)^N q_i^2 \
$
where $A(i)$ can be computed as follows:
$
  A(i) = A(i + 1) + q_(i + 1) e^(-h z_(i + 1)), A(N) = 0
$
Thus it only takes $N$ steps to compute $A(i)$, and another $N$ steps to compute $S$ with $A$.

#pagebreak()

=== SOEwald2D

Summing up all Fourier modes directly, the complexity is *$O(N^(1.4))$*.

#figure(
  image("figs/soe_accuracy_energy.png", width: 600pt),
  caption: [#text(15pt)[Accuracy of SOEwald2D method for the energy.]],
)

#pagebreak()

#figure(
  image("figs/soe_accuracy_force.png", width: 700pt),
  caption: [#text(15pt)[Accuracy of SOEwald2D method for the force.]],
)

#pagebreak()

=== SOEwald2D + RBE = RBSE2D

Recall the energy expression:
$
  U_l^h & = sum_(h) e^(-h^2 / (4 alpha^2)) underbrace((sqrt(pi) alpha) /  (L_x L_y) sum_(i, j) q_i q_j e^(i h rho_(i j)) / h sum_(l=1)^M w_l (e^(- s_l alpha z) / (alpha s_l + h) - e^(- s_l alpha z) / (alpha s_l - h) + (2 alpha s_l e^(-h z)) / ((alpha s_l)^2 - h^2)), "Linear complexity") \ 
  & = sum_h e^(-h^2 / (4 alpha^2)) f(h)
$
Use importance sampling instead of direct sum:
$
  sum_h e^(-h^2 / (4 alpha^2)) f(h) approx (sum_h e^(-h^2 / (4 alpha^2)))/P sum_(h in cal(K)_P) f(h)
$
where $cal(K)_P$ is the set of $P$ Fourier modes sampled from $cal(P)(h) ~ e^(-h^2 / (4 alpha^2))$, where $P$ is a $O(1)$ constant.

The resulting stochastic method is called RBSE2D method, with *$O(N)$* complexity.

#pagebreak()

=== Numerical results

#figure(
  image("figs/soe_md.png", width: 500pt),
  caption: [#text(15pt)[Results of coarse-grained MD simulation of 1 : 1 electrolytes in the NVT ensemble with RBSE2D method.]],
)

#pagebreak()

#figure(
  image("figs/soe_runtime.png", width: 450pt),
  caption: [#text(15pt)[Time cost of SOEwald2D and RBSE2D methods for different system sizes.]],
)

== Quasi-Ewald method

=== Negatively confined systems

In previous works, we assume that $abs(gamma) <= 1$, which requires all dielectric constants are positive.
In recent years, metamaterials characterized by negative dielectric constants have been widely studied @cheng2017tunable @xu2020polyaniline @xie2022recent.

Assume a system with $epsilon_c >0$, but confined by a dielectric medium with $epsilon_u, epsilon_d < 0$, then the system is *negatively confined*, with $abs(gamma) > 1$.

Recall that the image charge of the $l$-th level:
$
  q_+^l = q gamma_d^(ceil(l/2)) gamma_u^(floor(l/2)) ", and " 
  q_-^l = q gamma_d^(floor(l/2)) gamma_u^(ceil(l/2)) 
$
the image charge series is divergent, resulting in *failure of the image charge method*.

#pagebreak()

=== Green's function in Fourier space

Another way is to solve the Poisson equation in Fourier space.

Applying Fourier transform to the Poisson equation in $x$ and $y$ directions, then
$
  gradient^2 G(r, r') = delta(r - r') arrow integral_(bb(R)^2) gradient^2 G(r, r') e^( - i h rho) d rho = integral_(bb(R)^2) delta(r - r') e^( - i h rho) d rho
$
is transformed to
$
  partial_z^2 hat(G)(h, z, r') - h^2 hat(G)(h, z, r') = e^(-i h rho') delta(z - z') \
  G(r, r') = 1 / (4 pi^2) integral_(bb(R)^2) hat(G)(h, z, r') e^(i h rho) d h
$

Solution to the above equation is
$
  hat(G)(h, z, r') = e^(-i h rho') (e^(-h |z - z'|) + gamma_d e^(-h (z + z')) + gamma_u e^(-h (2H - z - z')) + gamma_d gamma_u e^(-h (2H - |z - z'|))) / (2h(1 - gamma_u gamma_d e^(-2 h H)))
$
Notice that when $h_c = ln(gamma_u gamma_d) / (2H)$, $hat(G)$ diverges.

#pagebreak()

=== Singular substraction

With $hat(G)$, we have
$
  G(r, r') = & 1 / (4 pi^2) integral_(bb(R)^2) hat(G)(h, z, r') e^(i h rho) d h \
  = & integral_(bb(R)^2) (f(h, z, z')) / (8 pi^2 h(1 - gamma_u gamma_d e^(-2 h H))) e^(i h (rho - rho')) d h = 1 / (4 pi) integral_(bb(R)) (f(h, z, z') J_0(h Delta rho)) / (1 - gamma_u gamma_d e^(-2 h H)) d h,
$
which is a singular integral diverging at $h = h_c$, and
$
  f(h, z, z') = e^(-h |z - z'|) + gamma_d e^(-h (z + z')) + gamma_u e^(-h (2H - z - z')) + gamma_d gamma_u e^(-h (2H - |z - z'|)).
$

Notice that
$
  lim_(d h arrow 0) 1 / (1 - gamma_u gamma_d e^(-2 (h_c + d h) H)) = lim_(d h arrow 0) 1 / (1 - e^(2 H d h)) = 1 / (2 H d h)
$
thus $h_c$ is a 1st order pole.

#pagebreak()

Split the integral into two parts:
$
  I_1 = 1 / (4 pi) integral_(bb(R)) ((f(h, z, z') J_0(h Delta rho) - f(h_c, z, z') J_0(h_c Delta rho)) / (1 - e^(-2 H (h - h_c))) + f(h_c, z, z') J_0(h_c Delta rho)) d h \
  I_2 = (f(h_c, z, z') J_0(h_c Delta rho)) / (4 pi) integral_(bb(R)) (e^(-2 H (h - h_c))) / (1 - e^(-2 H (h - h_c))) d h =  (f(h_c, z, z') J_0(h_c Delta rho)) / (8 pi H) ln(gamma_u gamma_d - 1)
$

#figure(
  image("figs/singular_integral.png", width: 350pt),
  caption: [#text(15pt)[Integrand of $I_1$ and $I_2$ as function of $h$, with $h_c = 1$.]],
)

#pagebreak()

With the above two integrals, we can compute the force between two charges in the system confined by dielectric slabs with different $gamma_u = gamma_d = gamma$.
When $abs(gamma) < 1$, the results are compared against that of the ICM.

#figure(
  image("figs/qem_force.png", width: 700pt),
  caption: [#text(15pt)[(a) Force between two charges in systems confined by dielectric slabs with different $gamma_u = gamma_d = gamma$, (b) the corresponding field lines.]],
)

#pagebreak()

=== Quasi-Ewald splitting

Spliting short-range and long-range contributions is needed for simulation:
$
  delta(bold(r)) = sigma_s (bold(r)) + sigma_l (bold(r)) \ 
  G (bold(r), bold(r)') = 1 / (4 pi^2) integral.triple e^(i bold(h) dot bold(rho)) hat(G)(bold(h), z, z'') hat(sigma) (bold(h), z'' - z') d bold(h) d z'' 
$
To simplify the computation, we introduce the quasi-Ewald splitting:
$
  delta(bold(r)) = underbrace(delta(bold(r)) - e^(-alpha^2 bold(rho)^2) delta(z), #text[$sigma_s$]) + underbrace(e^(-alpha^2 bold(rho)^2) delta(z), #text[$sigma_l$])
$
so that the integral is simplified.

#figure(
  image("figs/qem_splitting.png", width: 500pt),
  caption: [#text(15pt)[Illustration of the quasi-Ewald splitting.]],
)

#pagebreak()

The accuracy of the quasi-Ewald method is compared with that of the ICM-Ewald2D method.

#figure(
  image("figs/qem_accuracy.png", width: 450pt),
  caption: [#text(15pt)[Accuracy of the quasi-Ewald method for different $gamma_u = gamma_d = gamma$, box size is fixed to $(100, 100, 50)$, $N = 100$, $E$ is the parameter for accuracy.]],
)

#pagebreak()

Further combined with random batch sampling, we obtain a method with $O(N)$ complexity.

#figure(
  grid(columns: 2, gutter: 10pt,
    image("figs/qem_md.png", width: 400pt),
    image("figs/qem_time.png", width: 330pt),
  ),
  caption: [#text(15pt)[(Left) Distributions of ion density in $z$ for symmetric electrolytes containing 218 cations and 218 anions. (Right) time cost of QEM method for different system sizes.]],
)


#pagebreak()

== References

#bibliography("ref.bib", style: "apa", title: none)
