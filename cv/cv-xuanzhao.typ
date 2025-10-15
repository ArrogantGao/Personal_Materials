#import "@preview/modernpro-cv:1.0.2": *

#import "@preview/fontawesome:0.5.0": *

#show: cv-single.with(
  font-type: "PT Serif",
  continue-header: "false",
  name: [Xuanzhao Gao],
  address: [Center for Computational Mathematics, Flatiron Institute, New York, NY, USA],
  lastupdated: "true",
  pagecount: "false",
  date: [2025-10-15],
  contacts: (
    (text: [xgao\@flatironinstitute.org], link: "mailto:xgao@flatironinstitute.org"),
    (text: [github.com/ArrogantGao], link: "https://github.com/ArrogantGao"),
    (text: [arrogantgao.github.io], link: "https://arrogantgao.github.io/")
  ),
)

#set par(justify: true)
#show link: set text(blue)

// about
#section[Research Interests]
I am interested in applied and computational mathematics in general, working on the following topics:
- Boundary integral equations for complex geometries.
- Fast summation algorithms for long-range interactions.
- Tensor network-based algorithms for combinatorial optimization problems.
#sectionsep
#section("Education")
#education(
  institution: [Hong Kong University of Science and Technology],
  major: [P.h.D. in Individual Interdisciplinary Program, major in Applied Mathematics. Thesis advisors: Zecheng Gan, Jinguo Liu, and Yang Xiang],
  date: "2021 -- 2025",
  location: "Hong Kong SAR, China",
)
// Advisor: #link("https://zcgan.github.io/")[Prof. Zecheng Gan]\; Co-advisor: #link("https://giggleliu.github.io/")[Prof. Jin-Guo Liu] & #link("https://www.math.hkust.edu.hk/people/faculty/profile/maxiang/")[Prof. Yang Xiang].
#subsectionsep
#education(
  institution: [The University of Science and Technology of China],
  major: [B.S. in Condensed Matter Physics \& B.S. in Computer Science. Undergraduate thesis advisor: Xianhui Chen],
  date: "2017 -- 2021",
  location: "China",
)
#sectionsep
#section("Positions")
#job(
  date: "2025 -- now",
  position: [Flatiron Research Fellow],
  institution: [Center for Computational Mathematics, Flatiron Institute, Simons Foundation, New York, NY, USA],
)
#subsectionsep
#sectionsep
#section("Publications")
#subsection("Peer-reviewed Publications")

@icm_error

@gan2024random

@complexitypaper

@gan2024fast

@roa2024probabilistic

@gao2024broken

@nanolett2020
#subsectionsep
#subsection("Manuscripts Under Review")

@fssog

@optimalbranching

#subsectionsep

#subsection("In Draft (preprint available upon request)")

// @tensorbranching

@quasiewald

#pagebreak()

#section("Software Packages")
#oneline-title-item(
  title: [#link("https://github.com/HPMolSim/ExTinyMD.jl")[ExTinyMD.jl]],
  content: [A framework for molecular dynamics simulations.],
)
#subsectionsep
#oneline-title-item(
  title: [#link("https://github.com/HPMolSim/EwaldSummations.jl")[EwaldSummations.jl]],
  content: [A comprehensive implementation of the Ewald summation method for electrostatic interactions in both triply and doubly periodic systems with and without dielectric mismatches.],
)
#subsectionsep
  #oneline-title-item(
  title: [#link("https://github.com/HPMolSim/ChebParticleMesh.jl")[ChebParticleMesh.jl]],
  content: [A suite of highly efficient tools for the widely used Particle-Mesh methods applicable to systems with arbitrary dimensions and periodicity.],
)
#subsectionsep
#oneline-title-item(
  title: [#link("https://github.com/TensorBFS/TropicalNumbers.jl")[TropicalNumbers.jl]],
  content: [A refined implementation of the tropical semiring.],
)
#subsectionsep
#oneline-title-item(
  title: [#link("https://github.com/TensorBFS/CuTropicalGEMM.jl")[CuTropicalGEMM.jl]],
  content: [A GPU-accelerated implementation of the tropical matrix multiplication.],
)
#subsectionsep
#oneline-title-item(
  title: [#link("https://github.com/ArrogantGao/TreeWidthSolver.jl")[TreeWidthSolver.jl]],
  content: [A collection of tools for calculating the exact tree width and tree decomposition of a given graph.],
)
#subsectionsep
#oneline-title-item(
  title: [#link("https://github.com/OptimalBranching/OptimalBranching.jl")[OptimalBranching.jl]],
  content: [A framework for automated discovery of optimal branching rules for the branch-and-bound algorithm.],
)
#sectionsep
#section("Open Source Projects")

#education(
  institution: [#link("https://summerofcode.withgoogle.com")[Google Summer of Code 2024]],
  major: [Contributed to the project #link("https://summerofcode.withgoogle.com/programs/2024/projects/B8qSy9dO")["Tensor network contraction order optimization and visualization"] released by the Julia Language community in GSoC 2024.],
  location: "The Julia Language",
)
#subsectionsep
#education(
  institution: [#link("https://summer-ospp.ac.cn")[Open Source Promotion Plan 2023]], 
  major: [Contributed to the project #link("https://summer-ospp.ac.cn/2023/org/prodetail/23fec0105?lang=en&list=pro")["TropicalGEMM on GPU"] released by the JuliaCN community in OSPP 2023.],
  location: "JuliaCN",
)
#sectionsep
#section("Presentations and Posters")
#education(
  institution: [The 14th CSCM Annual Conference],
  major: [A Fast Spectral Sum-of-Gaussians Method for Electrostatic Summation in Quasi-2D Systems],
  date: "August 17-21, 2025",
  location: "Invited Talk",
)
#subsectionsep
#sectionsep
#education(
  institution: [JuliaCN Meetup 2024],
  major: [#link("https://raw.githubusercontent.com/ArrogantGao/my_presentations/main/pre/treewidth.pdf")[TreeWidthSolver.jl: From Treewidth to Tensor Network Contraction Order]],
  date: "Nov 2-3, 2024",
  location: "Invited Talk",
)
#sectionsep
#education(
  institution: [SciCADE 2024],
  major: [#link("https://raw.githubusercontent.com/ArrogantGao/my_presentations/main/pre/fast_algorithm_for_q2d_coulomb_systems.pdf")[Fast Algorithm for Quasi-2D Coulomb Systems]],
  date: "July 15-19, 2024",
  location: "Contributed Talk",
)
#sectionsep
#education(
  institution: [JuliaCN Meetup 2023],
  major: [#link("https://raw.githubusercontent.com/ArrogantGao/my_presentations/main/pre/CuTropicalGEMM.pdf")[How to Implement Generic Matrix-Mul with Generic Element Types on GPU?]],
  date: "Dec 9, 2023",
  location: "Contributed Talk",
)
#sectionsep
#education(
  institution: [ICIAM 2023],
  major: [Random Batch Quasi-Ewald Method for the Simulations of Charged Particles under Dielectric Confinement],
  date: "August 20-25, 2023",
  location: "Poster",
)
#sectionsep
#section("Skills")
#oneline-title-item(
  title: [Programming Languages],
  content: [Julia (proficient), Python, C/C++, CUDA],
)
#subsectionsep
#oneline-title-item(
  title: [Languages],
  content: [Mandarin Chinese (native), English (proficient)],
)

// Keep this at the end
#show bibliography: none
#bibliography("my.bib", style: "american-physics-society")
