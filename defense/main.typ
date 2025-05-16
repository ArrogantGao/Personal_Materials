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

#pagebreak()

== Self-introduction

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

Molecular dynamics simulation is a widely used method to study the motion of atoms and molecules, which plays a key role in many fields, including materials science, biophysics, and drug design.

#figure(
  grid(columns: 3, 
    image("figs/pnas.2102516118fig01.jpg", width: 450pt),
    h(50pt),
    image("figs/MD_workflow.png", width: 200pt),
  ),
  caption: [#text(15pt)[(Left) Study of the structure of SARS-CoV-2 replication-transcription complex by molecular dynamics simulation @brandon2021structural. (Right) Illustration of the workflow of molecular dynamics simulation.]],
)


== Quasi-2D charged systems

#timecounter(1)

Quasi-2D systems @mazars2011long are at the macroscopic scale in $x y$, but microscopic in $z$, so that are always modeled as doubly periodic in numerical simulations.
// Q2D systems are widely exist in nature and engineering, for example, cell membranes and electrolyte near surfaces.

#figure(
  image("figs/Q2D.png", width: 400pt),
  caption: [#text(15pt)[Illustration of a quasi-2D charged system.]],
)

Coulomb interaction plays a key role in nature, leading to effect such as ion transportation and self-assembly. 

However, the Coulomb interaction decays as $r^(-1)$ in 3D, so that it is long ranged and singular at $r=0$, which make such simulation computationally expensive.


== Coulomb interaction

#timecounter(1)

Coulomb interaction plays a key role in nature, leading to effect such as ion transportation and self-assembly.

#figure(
  image("figs/self-assembly.png", width: 450pt),
  caption: [#text(15pt)[Self-assembly of nanoparticles via Coulomb interaction @luijten2014.]],
)

However, the Coulomb interaction decays as $r^(-1)$ in 3D, so that is long ranged and singular at $r=0$, which make such simulation computationally expensive.
Complexity of a direct sum of the Coulomb interaction in doubly periodic systems is about $O(N^2 epsilon^(-1/3))$.

== Dielectric confinements

Polarizable surfaces also play a key role in nature, which arise naturally due to dielectric mismatch between different materials, leading to ion transport in nanochannels, pattern formation and self-assembly of colloidal and polymer monolayers.

#figure(
  image("figs/pattern.png", width: 400pt),
  caption: [#text(15pt)[Pattern formation in 2D dipolar systems due to dielectric mismatch @wang2019dielectric.]],
)





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

Another challenge is the polarization effect, various strategies have been developed in recent years, by:
- Introducing image charges @yuan2021particle @liang2020harmonic.
- Numerically solving the Poisson equation with interface conditions @nguyen2019incorporating @maxian2021fast @ma2021modified.
These methods are also accelerated by FFT or FMM to reach a complexity of $O(N log N)$ or $O(N)$.
However, the computational cost significantly increases compared to the homogeneous case, especially for strongly confined systems.
Rugious error analysis and parameter selection are also in lack.

== Ewald summations revisited

Key point: handle short-range and long-range contributions separately.

$
  1 / r = underbrace((1 - "erf"(alpha r)) / r, "short-range") + underbrace("erf"(alpha r) / r, "long-range") = "erfc"(alpha r) / r + "erf"(alpha r) / r
$

#align(center,
  image("figs/erf.png", width: 400pt),
)

#pagebreak()

The short-range part is summed in the real space, and the long-range part is summed in the Fourier space.
The total interaction energy can be computed as follows:
$
  U_s = 1 / 2 sum_(i,j=1)^N sum_(bold(m)) q_i q_j "erfc"(alpha |r_i - r_j + L_bold(m)|) / abs(r_i - r_j + L_bold(m))
$
$
  U_l = pi / (2 L_x L_y) sum_(i,j=1)^N q_i q_j sum_(bold(h) != bold(0)) e^(i bold(h) dot (r_i - r_j)) / h G_alpha (bold(h), z_i - z_j) - alpha / sqrt(pi) sum_(i=1)^N q_i^2 + J_0
$
where $bold(h) = ((2 pi m_x) / L_x, (2 pi m_y) / L_y)$ is the reciprocal lattice vector, and
$
  J_0 = -pi / (L_x L_y) sum_(i,j=1)^N q_i q_j ("erf"(alpha |z_i - z_j|) + alpha / sqrt(pi) e^(- alpha^2 (z_i - z_j)^2)) \
  G_alpha (bold(h), z) = xi^+ (bold(h), z) + xi^- (bold(h), z) = e^(-h z) "erfc"(h / (2 alpha) + alpha z) + e^(h z) "erfc"(h / (2 alpha) - alpha z)
$

The resulting method is the so-called Ewald2D method with *$O(N^2 log(epsilon))$* complexity.

== Ewald summations revisited

One way to accelerate the summation is transforming the system as triply periodic one by adding vacuum layer in the $z$ direction.

#align(center,
  image("figs/vacuum_layer.png", width: 200pt),
)

The summation is reformulated as:
$
  U_l = underbrace((2 pi) / (L_x L_y L_z) sum_(bold(k) != bold(0)) e^(- k^2 / (4 alpha^2)) / k^2 (sum_(i = 1)^N q_i e^(i k r_i))^2 - alpha / sqrt(pi) sum_(i=1)^N q_i^2, "Ewald3D summation") + underbrace(U_("YB") + U_("ELC") , "correction terms") + U_("Trap")
$
where $k = 2 pi (m_x / L_x, m_y / L_y, m_z / L_z)$ is the reciprocal lattice vector.
Complexity of EwaldELC method is *$O(N^(1.5))$*, and can easily be accelerated by FFT or FMM to reach linear complexity.

== Random batch Ewald revisited

Random batch Ewald method @Jin2020SISC is a stochastic method proposed to accelerate the Ewald3D summation.
The main idea is to sample the long-range part in the Fourier space using importance sampling instead of summing over all the Fourier modes.
$
  sum_(bold(k) != bold(0)) e^(- k^2 / (4 alpha^2)) / k^2 (sum_(i = 1)^N q_i e^(i k r_i))^2 approx (sum_(bold(k) != bold(0)) e^(- k^2 / (4 alpha^2))) / P sum_(k in cal(K)_P) 1 / k^2 (sum_(i = 1)^N q_i e^(i k r_i))^2,
$
where $cal(K)_P$ is the set of $P$ Fourier modes sampled from the distribution $cal(P) (k) ~ e^(- k^2 / (4 alpha^2))$.
$P$ is proven to be an $O(1)$ constant so that the total complexity is $O(N)$.

#figure(
  grid(columns:2, gutter: 25pt, 
    image("figs/rbe_accuracy.png", height: 220pt),
    image("figs/rbe_scale.png", height: 220pt),
  ),
  caption: [#text(15pt)[Accuracy and scalability of random batch Ewald method @liang2022superscalability.]],
)

== Image charge method revisited

In dielectrically confined Q2D systems, dielectric permittivity is a function of $z$:

$
  epsilon (z) = cases(
  epsilon_u "if" z > H,
  epsilon_c "if" 0 < z < H,
  epsilon_d "if" z < 0,
)
$

The governing equation of the pairwise electrostatic interaction $G(r, r')$:
$
  cases(
    gradient_r (epsilon(r) gradient_r G(r, r')) = delta(r - r') & "if" r in R^3,
    G(r, r') |_- = G(r, r') |_+ & "if" r in partial Omega_c,
    epsilon_c gradient_r G(r, r') |_- = epsilon_u gradient_r G(r, r') |_+ & "if" r in partial Omega_c \u{2229} partial Omega_u,
    epsilon_c gradient_r G(r, r') |_+ = epsilon_d gradient_r G(r, r') |_- & "if" r in partial Omega_c \u{2229} partial Omega_d,
    G(r, r') |_- = G(r, r') |_+ & "if" r in partial Omega_u,
  )
$

Define the dielectric reflection rate, $gamma_u$ and $gamma_d$:
$
  gamma_u = (epsilon_c - epsilon_u) / (epsilon_c + epsilon_u), " and "
  gamma_d = (epsilon_c - epsilon_d) / (epsilon_c + epsilon_d),
$

#pagebreak()

Assume that $epsilon_u, epsilon_d, epsilon_c$ are all positive, then $gamma_u, gamma_d in (-1, 1)$, and we get the so-called *image charge series* for the polarization potential:
$
  G(r, r') = 1 / (4 pi epsilon_c) [1 / abs(r - r') + sum_(l=1)^infinity (gamma_u^(l) / abs(r - r_+^(l) ') + gamma_d^(l) / abs(r - r_-^(l) '))]
$
where $r_(+)^(l) = (x, y, z_(+)^(l))$ and $r_(-)^(l) = (x, y, z_(-)^(l))$ are the positions of the $l$-th level image charges, and 
$
  z_+^(l) = (-1)^l z + 2 ceil(l/2) H,
  z_-^(l) = (-1)^l z - 2 floor(l/2) H.
$

#align(center,
  image("figs/icm_system.png", width: 500pt)
)

= Theoretical Analysis

== Error estimation

When combined with image charge method, the error of EwaldELC method is quite complex.

#figure(
  image("figs/error_yuan.png", width: 500pt),
  caption: [#text(15pt)[Relative precision in the electrostatic energy computed by ICM-PPPM as a function of the number of reciprocal reflections and the slab factor for a 2:1 salt confined between dielectric interfaces @yuan2021particle.]],
)

A rigorous error analysis is thus very important for the effective application of the EwaldELC method in systems with dielectric confinement.

#pagebreak()

In our work @icm_error, we show that the error of ICM-EwaldELC method is given as follows:
$
  cal(E) ~  underbrace(O(e^(-s^2) / (s^2)), "Ewald sum") + underbrace(O(abs(gamma_u gamma_d)^(floor((M + 1) / 2)) e^( - (4 pi H floor((M + 1) / 2)) / max(L_x, L_y) )), "Image charge") + underbrace(O(e^( - alpha^2 (L_z - H)^2 ) ), "Trapezoidal rule") \ + underbrace(O(e^( - (2 pi (L_z - H)) / max(L_x, L_y) ) + sum_(l=1)^M (abs(gamma_d^(ceil(l/2)) gamma_u^(floor(l/2))) + abs(gamma_d^(floor(l/2)) gamma_u^(ceil(l/2)))) e^( - (2 pi (L_z - (l + 1) H)) / max(L_x, L_y) )), "ELC term")
$
with the assumption that $abs(gamma_u) <= 1$ and $abs(gamma_d) <= 1$.

#pagebreak()

=== Image charge series

Error induced by truncation of the image charge series decays exponentially with the layer number $M$, the speed depends on $gamma$ and aspect ratio $H / max(L_x, L_y)$.
$
  cal(E) ~ O(abs(gamma_u gamma_d)^(floor((M + 1) / 2)) e^( - (4 pi H floor((M + 1) / 2)) / max(L_x, L_y) ))
$

#figure(
  image("figs/icm_image_charge_error.png", width: 550pt),
  caption: [#text(15pt)[Error of image charge series as a function of the layer number $M$ for different $gamma$ and aspect ratio $H / max(L_x, L_y)$, the dashed lines are fitted according to the theoretical prediction.]],
)

#pagebreak()

=== ELC term

First consider the case when $gamma_u = gamma_d = 0$, i.e., the system is homogeneous, then
$
  cal(E) ~ O(e^( - (2 pi (L_z - H)) / max(L_x, L_y)))
$
we thus define the padding ratio $R = (L_z - H) / max(L_x, L_y)$, and ELC term decays exponentially with $R$.

#figure(
  image("figs/icm_elc.png", width: 300pt),
  caption: [#text(15pt)[Error of ELC term as a function of the padding ratio $R$ for different aspect ratio.]],
)

#pagebreak()

If $abs(gamma_u gamma_d e^((4 pi H) / max(L_x, L_y))) <= 1$, the ELC term is a constant independent of $M$:
$
  cal(E) ~ O((gamma_u e^((2 pi H) / max(L_x, L_y)) + gamma_d e^((2 pi H) / max(L_x, L_y)) + 2)e^( - (2 pi (L_z - H)) / max(L_x, L_y)))
$

#figure(
  image("figs/icm_elc_decay.png", width: 600pt),
  caption: [#text(15pt)[Error of ELC term as a function of $M$ and $R$. We fix $L_x = L_y = 10$, $H = 0.5$ and $gamma_u = gamma_d = 0.6$ so that $abs(gamma_u gamma_d e^((4 pi H) / max(L_x, L_y))) < 1$]],
)

#pagebreak()

Interestingly, if $abs(gamma_u gamma_d e^((4 pi H) / max(L_x, L_y))) > 1$, the ELC term grows exponentially to $M$:
$
  cal(E) ~ O(abs(gamma_u gamma_d e^((4 pi H) / max(L_x, L_y)))^(M / 2) e^(- (2 pi (L_z - H)) / max(L_x, L_y)))
$

#figure(
  image("figs/icm_elc_increase.png", width: 600pt),
  caption: [#text(15pt)[Error of ELC term as a function of $M$ and $R$. We fix $L_x = L_y = 10$, $H = 0.5$ and $gamma_u = gamma_d = 1$ so that $abs(gamma_u gamma_d e^((4 pi H) / max(L_x, L_y))) > 1$]],
)

#pagebreak()

Our results well explain the behavior of the error of ICM-EwaldELC method in the literature @yuan2021particle.

#figure(
  image("figs/icm_error_yuan.png", width: 700pt),
  caption: [#text(15pt)[Repeat of the results in @yuan2021particle with fitted lines according to our theoretical prediction, where $L_x = L_y = 15$, $H = 5$. The system is a 2:1 electrolyte with $39$ particles.]],
)

== Parameter selection

Given the system size $L_x, L_y, H$, the dielectric constants $epsilon_u, epsilon_d, epsilon_c$, and desired precision $epsilon$, we can select the parameters as follows:

*Step 1*:
Select the layer number $M$:
$
  M ~ (2 log(epsilon) - (4 pi H) / max(L_x, L_y) - log(abs(gamma_u gamma_d))) / (log(abs(gamma_u gamma_d)) - (4 pi H) / max(L_x, L_y))
$

*Step 2*:
Select height of the vacuum layer $L_z$:
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

== Overview

=== Sum-of-Exponential Ewald2D method
- For the Q2D systems in homogeneous media.
- Reduce the complexity of double summation to $O(N)$ via SOE approximation.
- Numerically stable.

=== Random batch Ewald2D method
- For the dielectrically confined Q2D systems.
- Reduce the complexity of ICM-EwaldELC based method from $O(M N)$ to $O(M + N)$.
- Parallelizable and scalable.

=== Quasi-Ewald method
- For the negatively confined Q2D systems.
- Remove the divergence due to negative dielectric constant via singularity subtractions.


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

== Random batch Ewald2D method

rbe2d

== Quasi-Ewald method

qem

== Short summary

In summary, we developed a series of methods for simulating quasi-2D charged systems.

=== Advantages

- Mesh free, only use simply data structure.
- Not sensitive to the aspect ratio of the system.
- Parallelizable and scalable.

=== Disadvantages

- Stochastic, replying on law of large numbers.
- Limited to planar systems with uniform dielectric constant.
- Strong convergence and ergodicity remaining an open question.

= Applications in MD Simulations

== All-atom simulations of SPC/E water

Our method is applied to the all-atom simulations of SPC/E water confined by slabs, where $L_x = L_y = H = 55.9 angstrom$ and consists of 17496 atoms.

#figure(
  image("figs/spce.png", width: 500pt),
  caption: [#text(15pt)[All-atom simulations of SPC/E water confined by slabs.]],
)

#pagebreak()

We further applied our method to the simulations of SPC/E water confined by slabs with different dielectric constants, where $gamma_u = gamma_d = -0.5$, the system consists of 53367 atoms.

#figure(
  image("figs/spce_dielectric.png", width: 500pt),
  caption: [#text(15pt)[All-atom simulations of SPC/E water dielectrically confined by slabs.]],
)

// #pagebreak()

// #figure(
//   grid(rows:2, gutter: 15pt, 
//     image("figs/spce_time.png", width: 650pt),
//     image("figs/spce_time_dielectric.png", width: 650pt),
//   ), 
//   caption: [#text(15pt)[Time cost of all-atom simulations of SPC/E water in homogeneous media (upper) and confined by slabs with different dielectric constants (lower).]],
// )

== Broken symmetries via dielectric confinements

With the help of the quasi Ewald method, we observed spontaneous symmetry broken via dielectric confinements when $abs(gamma) > 1$.

#figure(
  image("figs/ssb.png", width: 700pt),
  caption: [#text(15pt)[Spontaneous symmetry broken via dielectric confinements.]],
)

#pagebreak()

#figure(
  image("figs/ssb_force.png", width: 600pt),
  caption: [#text(15pt)[Detailed force analysis of the symmetry broken phase.]],
)

// = Introduction to Other Works

// == Fast spectral sum-of-Gaussians method

// Based on the sum-of-Gaussian approximation of the Coulomb kernel and non-uniform fast Fourier transform, we developed a fast spectral method for the Q2D systems in homogeneous media.
// We split the Coulomb kernel into three parts:
// $
//   1 / r approx underbrace(( 1 / r - sum_(l = 0)^M w_l e^(- r^2 / s_l^2)) bb(1)_(r < r_c), "near-field")+ underbrace(sum_(l = 0)^m w_l e^(- r^2 / s_l^2), "mid-range") + underbrace(sum_(l = m + 1)^M w_l e^(- r^2 / s_l^2), "long-range")
// $

// The mid-range potential is computed by a standard Fourier spectral solver #footnote(text(12pt)[#link("https://github.com/HPMolSim/ChebParticleMesh.jl")],) with *little zero padding* ($lambda_z < 2$ for double precision).
// *No need of the kernel truncation* in the free direction due to the smoothness and separability of the Gaussian.

// The extremely smooth long-range Gaussians are interpolated on the Chebyshev proxy points in $z$, similar to that of the periodic FMM (Pei, Askham, Greengard & Jiang, 2023), and only *$O(1)$ number of Chebyshev points are required*.

// #pagebreak()


// The method #footnote(text(12pt)[#link("https://github.com/HPMolSim/FastSpecSoG.jl")],) is benchmarked on the following systems:
// - Cubic systems with fixed aspect ratio equals to $1$.
// - Strongly confined systems with fixed $L_z$, aspect ratio up to $10^(3.5)$

// #figure(
//   image("figs/sog_benchmark.png", width: 370pt),
//   caption: [#text(15pt)[Error and time cost for the SOG method in the (a,c) cubic and (b,d) strongly confined systems.]],
// )

= Summary and Outlook

== Summary

During my PhD, I focused on confined quasi-2D charged systems, including:

- Theoretical analysis of Ewald summation for dielectric-confined planar systems.

- Development of novel fast algorithms for simulating quasi-2D Coulomb systems, including both stochastic methods and spectral methods, for non-dielectrically confined systems, dielectrically confined systems, and negatively confined systems.

- Applications in MD simulations, including all-atom simulations of SPC/E water, and the observation of spontaneous symmetry broken solely via dielectric confinements for the first time.


== Future works

- Extend our methods to more complex systems, including systems with complex geometries such as curved surfaces and polarizable particles.

- Extend our methods to different interaction kernels, such as the dipolar interactions, the Stokes interactions and the Yukawa interactions.

- Application to the simulation of real world systems, such as biomolecular systems and colloidal systems.

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

// === My packages

- *ChebParticleMesh.jl*#footnote(text(12pt)[#link("https://github.com/HPMolSim/ChebParticleMesh.jl")],): Toolkits for particle mesh methods (type-1 and type-2 NUFFT).
- *CuTropicalGEMM.jl*#footnote(text(12pt)[#link("https://github.com/TensorBFS/CuTropicalGEMM.jl")],): Custom GPU kernel for tropical matrix multiplication.
- *TreeWidthSolver.jl*#footnote(text(12pt)[#link("https://github.com/ArrogantGao/TreeWidthSolver.jl")],): Solving the treewidth problem (supported by GSoC 2024).
- *FastSpecSoG.jl*#footnote(text(12pt)[#link("https://github.com/HPMolSim/FastSpecSoG.jl")],): Implementation of the fast spectral SOG method.
- *EwaldSummations.jl*#footnote(text(12pt)[#link("https://github.com/HPMolSim/EwaldSummations.jl")],): Various Ewald summation methods with parallelization.
- *OptimalBranching.jl*#footnote(text(12pt)[#link("https://github.com/OptimalBranching/OptimalBranching.jl")],): Implementation of the optimal branching algorithm.

// === Contributions to popular packages

// - *OMEinsum.jl*#footnote(text(12pt)[#link("https://github.com/under-Peter/OMEinsum.jl")],) (185 stars) and its backend *OMEinsumContractionOrders.jl*#footnote(text(12pt)[#link("https://github.com/TensorBFS/OMEinsumContractionOrders.jl")],): Optimizing the tensor network contraction order and contracting the tensor network.

= Acknowledgements

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

// = Appendix <touying:unoutlined>

// == Supplementary results

// This is a result.

// #pagebreak()

== References

#bibliography("ref.bib", style: "apa", title: none)
