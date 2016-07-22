project: JOREK
project_website: http://jorek.eu
summary: JOREK is a non-linear extended MHD code for toroidal X-point geometries
author: the JOREK team
project_dir: .
page_dir: pages
output_dir: ./doc
exclude_dir: ./doc
	     ./models
docmark: <
docmark_alt: *
predocmark: >
predocmark_alt: #
extra_filetypes: c //
		 sh #
graph: true
coloured_edges: true

The non-linear extended MHD code JOREK resolves realistic toroidal X-point geometries with a C1 continuous flux-surface aligned grid including main plasma, scrape-off layer and divertor region. It is based on robust fully implicit numerics, and includes divertor boundary conditions, 3D resistive wall effects, two-fluid effects and neoclassical flows.

The well established physics and numerics community around JOREK has strong connections to the relevant experiments, the ITER Organization and the respective ITPA Topical Groups.


## Key Physics Applications

* Edge Localized Modes (ELMs)
* Pellet triggering of Edge Localized Modes
* Error field penetration for resonant magnetic perturbations (RMPs)
* Mitigation / suppression of ELMs via RMP fields
* Massive gas injection triggered disruptions
* Stabilization of Tearing Modes via ECCD
* Vertical Displacement Events
* Resistive Wall Modes
* QH-Mode

## Physics Models

* Reduced and full MHD models
* Two-fluid and neoclassical effects
* Divertor model
* Pellet ablation model
* Neutrals model for Deuterium MGI — impurity MGI model under development
* Resistive wall extension — inclusion of Halo currents under development
* Electron Cyclotron Current Drive (ECCD) model
* Particle in cell model under development

## Numerics

* Flux-aligned 2D Bezier finite elements — generalization under development
* Toroidal Fourier expansion — generalization under development
* Fully implicit time stepping
* GMRES and Newton iterations
* Physics based preconditioning — Jacobian free preconditioning under development
* Taylor-Galerkin stabilization
* MPI + OpenMP parallelization
