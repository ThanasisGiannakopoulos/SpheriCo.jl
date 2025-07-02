---
title: 'SpheriCo.jl: Spherical Collapse in classical and semiclassical gravity with Julia'
tags:
  - Julia
  - numerical relativity
  - semiclassical gravity
  - black holes
authors:
  - name: Thanasis Giannakopoulos
    orcid: 0000-0002-3055-9652
    affiliation: "1, 2"
affiliations:
 - name: School of Physics and Astronomy, University of Nottingham,
   index: 1
 - name: Nottingham Centre of Gravity, University of Nottingham
   index: 2
date: 3 July 2025
bibliography: paper.bib
---

# Summary

SpheriCo.jl is an open-source Julia package for simulating gravitational collapse in spherical symmetry, in both classical and semiclassical gravity. In the classical setup, it solves the gravitational collapse of a massless scalar field in general relativity. This problem has been extensively studied and led to the discovery of gravitational critical phenomena [@Choptuik:1992jv], with implications for fundamental questions such as the formation of naked singularities [@Gundlach:2007gc; @christodoulou1999instability]. In the semiclassical setting, the scalar field is promoted to a quantum operator, which is expanded into modes. The accuracy of the numerical approximation improves as the number of modes increases [@Berczi:2021hdh]. However, the equations satisfied by the different modes are stiff and can lead to code instabilities. SpheriCo.jl offers a solution to this problem and enables the study of semiclassical gravitational phenomena, such as computation of Hawking correlators [@Berczi:2024yhb].

# Statement of need

The primary motivation behind the development of SpheriCo.jl is the simulation of semiclassical gravitational collapse. When quantum effects of matter are taken into account, black holes can exhibit novel phenomena absent in the classical setting, such as Hawking radiation [@Hawking:1975vcx]. Numerical simulations can be a valuable tool in our efforts to understand these phenomena better and can act complementary to analytical calculations, which often assume a static classical geometry [@Levi:2015eea; @Lowe:1992ed]. To the best of the author's knowledge, the first numerical simulations of semiclassical gravitational collapse in spherical symmetry in four spacetime dimensions were presented in [@Berczi:2020nqy]. The code used in that work is available in [@Berczi_codes] (specifically, the file `semiclassical_collapse_Alcubierre.c`; for details on other codes, see [@Berczi:2024mav]). The main limitation of this earlier implementation arises from the form of the equations governing the quantum modes, which become unstable near the origin of the radial domain. As a result, simulations were restricted to short timescales and required particularly strong initial data (i.e. an apparent horizon forms very fast). SpheriCo.jl mitigates this problem by utilizing the summation-by-parts (SBP) operators developed by [@Gundlach:2010et]. This significantly improves numerical stability and enables longer simulations with a broader range of initial data.

In the semiclassical setup, matter is quantized while the spacetime geometry remains classical. As a result, a reliable and accurate solution to the classical gravitational collapse problem is a prerequisite for the functionality of SpheriCo.jl. Although this setup has been extensively studied, there are relatively few open-source codes that implement it; to the best of the author’s knowledge, only [@OllinSphere-BiB] and [@engrenage] are publicly available. The equations governing the classical collapse of a massless scalar field contain terms proportional to $1/r$, which can lead to numerical instabilities near $r=0$. While there is substantial literature addressing this challenge—such as the use of Evans’ method on a centered grid [@SuarezFernandez:2020wqv; @Evans_PhD], or the combination of staggered grids and artificial dissipation [@AlcubierreNumRel]—achieving a stable and convergent numerical implementation remains non-trivial. SpheriCo.jl can offer a relatively low-entry-level option for researchers interested in classical gravitational collapse, and/or a benchmark for the development of their own code. The package includes documentation, usage examples, and postprocessing tools to support this. Additionally, the use of the specific SBP operators to handle the $1/r$ terms is a unique feature compared to existing solutions.

# Key features

The use of the SBP operators of [@Gundlach:2010et] to resolve instabilities coming from terms of the form $\sim l/r$ is the main key feature of SpheriCo.jl. The different quantum modes are labelled by two positive integer numbers $k,l$, which increase as more modes are included in the simulation. The goal is to include as many of these modes as possible, which results in an increasing set of evolved variables, with increasingly stiffer equations that become unstable faster. Currently, SpheriCo.jl utilizes a second-order accurate version of these SBP operators, and can achieve the expected second order convergence (see [@Berczi:2024yhb] for relevant tests). In comparison, the code used in [@Berczi:2020nqy] uses higher order standard finite difference methods in an attempt to mitigate instabilities. These are much more expensive (require more grid points for their numerical approximation) and achieve  numerical stability for much shorter time and restricted initial data, while recovering lower order convergence rate.

In addition to numerical instabilities, another challenge of the semiclassical setup is the growing number of evolved variables with increasing modes, which can make the simulations slow. To increase the speed, the package uses the third-order accurate Adams-Bashford time integrator. The key difference in comparison to e.g. a commonly used Runge-Kutta integrator, is that the right-hand-side of the equations for the quantum modes needs to be calculated only once for each timestep, in comparison to three (for a third-order accurate Runge-Kutta). This feature, together with the second-order accurate SBPs and the native Julia parallelization on multiple threads, make non-trivial SpheriCo.jl simulations possible on workstations and laptops.

# Research projects to date using SpheriCo.jl

SpheriCo.jl has been used in [@Berczi:2024yhb] to explore semiclassical phenomena around dynamically forming apparent horizons (supercritical solutions), that are persived as dynamically forming black holes. More specifically, two-point correlation functions of Hawking pairs were calculated, with results that hint at a non-trivial correlation across the horizon of Hawking quanta. 

![The real part of the two-point correlation function for a subritical (no apparent horizon) and supercritical (apparent horizon forms) solutions. Taken from [@Berczi:2024yhb]](./quantum_correlators.png)

# Limitations and possible future improvements

In its current state, SpheriCo.jl is designed primarily to perform long, stable, and convergent simulations of semiclassical gravitational collapse. To achieve this, it combines second-order accurate SBP operators with a third-order accurate Adams-Bashforth time integrator, enabling faster and more computationally efficient simulations, though at the cost of accuracy. While the package supports significantly longer stable simulations and a broader range of initial data than previously possible, instabilities can still arise in certain regions of parameter space, particularly near criticality, the threshold between black hole formation and dispersion to flat space. In this regime, high accuracy near $r=0$ is crucial.

Mesh refinement is a valuable tool for overcoming this limitation, and SpheriCo.jl includes this feature. However, the current implementation is based on interpolation and projection to a smaller radial grid, which is stable only in the classical regime. A promising alternative would be to implement a formulation with a dynamical shift vector [@Rinne:2020asi], which could reduce computational cost and potentially enable near-critical explorations in the semiclassical setting.

Given the importance of accuracy for critical phenomena studies, implementing higher-order accurate SBP operators and time integrators is another possible improvement. Finally, controlling violations of the Hamiltonian and momentum constraints remains essential for obtaining reliable near-critical solutions. Even though the formulation used allows for damping of these constraints, in practise it does not perform well in this front. A possible reason for this failure is that it does not include the damping of reduction constraints, which has been shown to be essential [@Cors:2023ncc]. Enhancing the formulation with this feature is another possible future improvement.

# Acknowledgments

During the development of SpheriCo.jl, the author was supported by the Science and Technology Facilities Council (STFC) under Grant Nos. ST/V005596/1 and ST/X000672/1. The author also thanks Benjamin Berczi, Magdalena Eriksson, and Paul Saffin for their collaboration on [@Berczi:2024yhb], which both motivated and guided the development of this package.

# References