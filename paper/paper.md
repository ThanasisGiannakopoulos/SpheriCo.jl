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
    affiliation: "1; 2"
affiliations:
 - name: School of Physics and Astronomy, University of Nottingham
   index: 1
 - name: Nottingham Centre of Gravity, University of Nottingham
   index: 2
date: xx July 2025
bibliography: paper.bib
---

# Summary

<!-- SpheriCo.jl is an open-source Julia package for simulating gravitational collapse in spherical symmetry. The package is designed in a modular way, with the main module containing the infrastructure to solve partial differential equations with the method of lines (e.g. finite difference operators, time integrators etc.), and the submodules the elements related to the physical setup (e.g. the right-hand-side of the equations solved). There two submodules: the classical and the quantum. The first solves the gravitational collapse of a massless scalar field in general relativity. This setup has been extensively studying and led to the discovery of gravitational critical phenomena [@Choptuik:1992jv], which have implications on fundamental questions such as the formation of naked signularities [@Gundlach:2007gc;@christodoulou1999instability]. The second solves a semiclassical version of the problem, where the massless scalar field is promoted to a quantum operator. This operator is expanded in modes and its approximation improves with increasing number of modes [@Berczi:2021hdh]. The equations satisfied by the different modes are stiff and can lead to code instabilities. The package offers a solution to this problem and allows for explorations of semiclassical gravitational phenomena, such as Hawking correlators [@Berczi:2024yhb]. -->

SpheriCo.jl is an open-source Julia package for simulating gravitational collapse in spherical symmetry, in classical and semiclassical gravity. The classical setup solves the gravitational collapse of a massless scalar field in general relativity. This has been extensively studying and led to the discovery of gravitational critical phenomena [@Choptuik:1992jv], which have implications on fundamental questions such as the formation of naked signularities [@Gundlach:2007gc; @christodoulou1999instability]. In the semiclassical version of the problem, the massless scalar field is promoted to a quantum operator. This operator is expanded in modes and its numerical approximation improves with increasing number of modes [@Berczi:2021hdh]. However, the equations satisfied by the different modes are stiff and can lead to code instabilities. The package offers a solution to this problem and allows for explorations of semiclassical gravitational phenomena, such as Hawking correlators [@Berczi:2024yhb].

# Statement of need

The main motivation for the development of SpheriCo.jl is the semiclassical gravitational collapse setup. When quantum effects of matter are considered, black holes can produce exciting new phenomena that are not expected classically, such as Hawking radiation [@Hawking:1975vcx]. Numerical simulations can be a valuable tool in our efforts to understand these phenomena better and can act complementary to analytical calculations, which often assume a static classical geometry [@Levi:2015eea; @Lowe:1992ed]. To the best of the author's knowledge, the first numerical simulations of semiclassical gravitational collapse in spherical symmetry, in four spacetime dimensions, were presented in  [@Berczi:2020nqy]. The code used in this paper can be found in [@Berczi_codes] (that is the code ``semiclassical_collapse_Alcubierre.c'', for more details on the other codes see [@Berczi:2024mav]). The main limitation of this code is due the form of the equations satisfied by the quantum modes that lead to instabilities near the origin of the radial domain and allow only relatively short simulations with very strong initial data. SpheriCo.jl mitigates this problem by  utilizing the summation-by-parts (SBP) operators developed by [@Gundlach:2010et] and greatly enhances our ability to perform longer similations with more diverse initial data.

<!-- The main limitation of this code is due the form of equations satisfied by the quantum modes that lead to instabilities near the origin of the radial domain.  -->

<!-- More specifically, each quantum mode is labelled by two integers $k,l$. In spherical polar coordinates, the evolution equation satisfied by each mode is a wave equation that involves a term of the form $\sim l/r$. As more modes are included in the simulation, the number $l$ increases, making the equation more stiff and the accurate approximation of the quantum mode near $r=0$ harder.

To mitigate this problem, SpheriCo.jl utilizes the the second order accurate version of the summation-by-parts (SBP) operators developed by [@Gundlach:2010et]. This is the key factor that allows SpheriCo.jl to achieve stable simulations for much longer than previously possible, with more diverse initial data. This makes the simulations much more economic in comparison to the earlier code, which was trying to mitigate instabilities by using expensive higher order standard finite difference operators. The combination of SBP operators, thread parallelization for the quantum mode equations, and the use of third order Adams-Bashford time integrator []. -->

In the semiclassical setup matter is quantized but the geometry is treated classically. This implies that providing a satisfying solution to the classical gravitational collapse is a prerequisite for SpheriCo.jl. Even though this setup is well studied, there are not many open-source codes that solve it; with [@OllinSphere-BiB] and [@engrenage] being the only ones available to the best of the author's knowledge. The equations satisfied by the classical massless scalar field still include terms of the form $\sim 1/r$, which can cause instabilities near $r=0$. While there is extensive literature on this problem and various different approaches (see e.g. [@SuarezFernandez:2020wqv] for the use of Evans' method [@Evans_PhD] on a centered grid, and [@AlcubierreNumRel] for the combination of a staggered grid and artificial dissipation), achieving a succesful (stable and converging) numerical implementation is non-trivial. SpheriCo.jl can offer a relatively low-entry-level option for researchers interested in classical gravitational collapse, and/or a benchmark for the development of their own code. The documentation, examples, and postprocessing tools of the package can be helpful in this process. In addition, the use of the specific SBP operators to treat the $1/r$ terms of the classical equations is a unique feature in comparison to other available options.

<!-- 
The other direction for doing such simulations is to use full 3D codes and impose extra symmetry, which is expensive and harder to do for less experienced researchers. -->

# Key features

The use of the SBP operators of [@Gundlach:2010et] to resolve instabilities coming from terms of the form $\sim l/r$ is the main key feature of SpheriCo.jl. The different quantum modes are labelled by two positive integer numbers $k,l$, which increase as more modes are included in the simulation. The goal is to include as many of these modes as possible, which results in an increasing set of evolved variables, with increasingly stiffer equations that become unstable faster. Currently, SpheriCo.jl utilized a second-order accurate version of these SBP operators, and can achieve the expected second order convergence predicted by this numerical method (see [@Berczi:2024yhb] for relevant tests). In comparison, the code used in [@Berczi:2020nqy] uses higher order standard finite difference methods to obtain better numerical stability. These are much more expensive (require more grid points for their numerical approximation) and achieve  numerical stability for much shorter time and restricted initial data.

In addition to numerical instabilities, another challenge of the semiclassical setup is the growing number of evolved variables with increasing modes, which can make the simulations slow and expensive. To increase the speed of the simulation, SpheriCo.jl uses the third-order accurate Adams-Bashford time integrator. The key difference in comparison to e.g. a commonly used Runge-Kutta integrator, is that the right-hand-side of the equations for the quantum modes needs to be calculated only once for each timestep, in comparison to three (for a third-order accurate Runge-Kutta). This feature, together with the second-order accurate SBPs and the native Julia parallelization on multiple threads, make non-trivial SpheriCo.jl simulations possible on workstations and laptops.

# Research projects to date using SpheriCo.jl

SpheriCo.jl has been used in [@Berczi:2024yhb] to explore semiclassical phenomena around dynamically forming apparent horizons (supercritical solutions), that are persived as dynamically forming black holes. More specifically, two-point correlation functions of Hawking pairs were calculated, with results that hint at a non-trivial correlation across the horizon of Hawking quanta. 

![The real part of the two-point correlation function for a subritical (no apparent horizon) and supercritical (apparent horizon forms) solutions. Taken from [@Berczi:2024yhb]](./quantum_correlators.png)

# Limitations and possible future improvements

At its current state, SpheriCo.jl is designed with the main goal of performing long, stable, and converging simulations of the semiclassical gravitational setup. To do this, second-order accurate SBPs and third-order accurate Adams-Bashford time integrator are combined, to achieve cheaper and faster simulations, but at the expense of accuracy. Even though the package allows for much longer stable simulations than previously possible, and wider set of initial data, there are still regions of the parameter space where instabilities can develop. This limitation is particularly prominent close to criticallity (the region in parameter space between black hole formation and dispersion of initial data to flat space). In this case, high accuracy at the region near $r=0$ for long simulations is necessary.

Mesh refinement can be very useful in addressing this limitation and SpheriCo.jl has this feature. However, it is implemented only in a numerical fashion that includes interpolation and projection to a smaller radial grid, which is stable only for the classical setup. An alternative way that could improve the package would be to implement a formulation with a dynamical shift [@Rinne:2020asi]. This would be cheaper and may allow for near critical explorations in the semicalssical setup. Since high accuracy is very important in studying gravitational critical phenomena, implementing higher-order accurate SBPs and time integrators is another possible improvement. Finally, control of the violation of the Hamiltonian and momentum constraints is also important in obtaining good near critical solutions. Even though SpheriCo.jl implements a formulation that allows for damping of these constraints, it does not include the damping of reduction constraints (that come from introducing first order derivatives of a variable as evolved variables). This has also been shown to be essential in controlling the violation of the Hamiltonian and momentum constraints [@Cors:2023ncc], and is another possible future improvement.

# Acknowledgments

During the development of SpheriCo.jl the author has been supported by the Science and Technology Facilities Council (STFC) [Grant Nos. ST/V005596/1 and ST/X000672/1]. The author would like to thank Benjamin Berczi, Magdalena Eriksson, and Paul Saffin for the collaboration in [@Berczi:2024yhb], which motivated and drove the development of the package.

# References