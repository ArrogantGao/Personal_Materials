#import "@preview/modernpro-cv:1.0.2": *

#import "@preview/fontawesome:0.5.0": *

#show: cv-single.with(
  font-type: "PT Serif",
  continue-header: "false",
  name: [Xuanzhao Gao],
  address: [Hong Kong University of Science and Technology, Hong Kong SAR, China],
  lastupdated: "true",
  pagecount: "true",
  date: [2025-07-10],
  contacts: (
    (text: [xgaobj\@connect.ust.hk], link: "mailto:xgaobj@connect.ust.hk"),
  ),
)

#show link: set text(blue)

#section[Research Interests]
I am interested in applied and computational mathematics in general, with a particular focus on developing efficient numerical algorithms for modeling and simulating complex systems. 
Specifically, I am engaged in research on fast summation algorithms tailored for long-range interactions. 
I have also concentrated on algorithms for combinatorial optimization problems.

#sectionsep

#section("Education")
#education(
  institution: [Hong Kong University of Science and Technology],
  major: [P.h.D. Majoring in Applied Mathematics, under the Individual Interdisciplinary Program (Advanced Materials) \ Supervised by #link("https://zcgan.github.io/")[Prof. Zecheng Gan]],
  date: "2021 -- now",
  location: "Hong Kong SAR, China",
)
#subsectionsep
#education(
  institution: [The University of Science and Technology of China],
  major: [B.S. Majoring in Applied Physics],
  date: "2017 -- 2021",
  location: "China",
)

#sectionsep
#section("Publications")

@gan2024random

@icm_error

@complexitypaper

@gan2024fast

@tensorbranching

@fssog

@roa2024probabilistic

@gao2024broken

@nanolett2020

#sectionsep
#section("Presentations and Posters")
#education(
  institution: [SciCADE 2024, Singapore],
  major: [Fast Algorithm for Quasi-2D Coulomb Systems],
  date: "July 15-19, 2024",
  location: "Contributed Talk",
)
#sectionsep
#education(
  institution: [ICIAM 2023, Japan],
  major: [Random Batch Quasi-Ewald Method for the Simulations of Charged Particles under Dielectric Confinement],
  date: "August 20-25, 2023",
  location: "Poster",
)

#sectionsep
#section("Skills")
#oneline-title-item(
  title: [Programming Languages],
  content: [Julia (proficient), Python, C/C++],
)
// #subsectionsep
// #oneline-title-item(
//   title: [Languages],
//   content: [Mandarin Chinese (native), English (proficient)],
// )

// Keep this at the end
#show bibliography: none
#bibliography("./my.bib")
