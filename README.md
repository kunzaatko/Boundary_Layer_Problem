# A Tale of a Boundary-Layer Problem

This repository holds references and slides that I used for a talk at a university seminar called
[*Problems of Contemporary Mathematics*](https://people.fjfi.cvut.cz/vybirja2/Seminar/index.html). The
ideal length of the talk is about 2 hours. It is based on the SIAM article [*Surprises in a Classic Boundary-Layer Problem*](https://doi.org/10.1137/21M1436087) (Clark et al., 2023) and follows the structure of this article. The slides are made with `beamer` and include a number of motivational figures and some figures that appeared in the article and are reimplemented in `julia` also in this repo. A part of the talk is dedicated to illustrating how interactive figures can help understand the problem at hand and lead to pinpointing the correct approach to arrive at the solution. For a hint where to look for the sources mentioned above, take a look at the [`Contents`](#contents) section.

---
__Title:__ A Tale of a Boundary-Layer Problem

__Abstract:__ I will demonstrate the solution with all its steps (and missteps) to the boundary layer problem 
$$\varepsilon y'' = yy' - y,$$ 
such that 
$$y(0) = 1\quad \text{and} \quad y(1) = -1\quad \text{and where}\quad \varepsilon \ll 1.$$ 
The journey towards the solution will take us on an excurse through some basics of _perturbation theory_ and _asymptotic approximation matching_. I will demonstrate visualizations that will give you insight into the reasoning behind the steps taken towards the solution serving as an exhibition for the leverage in problem solving accessible to us through the computational power at our fingertips.

---

## Contents

- [`Slides_Interactive/`:](Slides_Interactive) Contains a [`Pluto`](https://github.com/fonsp/Pluto.jl) notebook with interactive figures
- [`Slides/figures/sources/`:](Slides/figures/sources/) Contains the source `julia` code used to generate the figures in the `beamer` presentation
- [`Slides/figures/`:](Slides/figures/) Contains the generated figures that appear in the `beamer` presentation
- [`Slides/src/`:](Slides/src/) Contains all the $\LaTeX$ source code that was used to compile the slides
- [`Slides/build/Slides_A-Tale-of-a-Boundary-Layer-Problem/`:](Slides/build/Slides_A-Tale-of-a-Boundary-Layer-Problem/) Contains the `PDF` slides

## Figures

Here is a showcase of some figures that were created for this project

| Comparison of the solutions (approximations) of the simple pendulum problems |Comparison of the solutions (approximations) to the projectile problem |
| :---: | :---: |
| ![Comparison of the solutions (approximations) of the simple pendulum problems](https://user-images.githubusercontent.com/56647779/236527361-fde1ffdd-d220-4914-b4da-9c116f9818e9.svg)|![Comparison of the solutions (approximations) to the projectile problem ](https://user-images.githubusercontent.com/56647779/236527612-40502f38-7cdc-4bb5-a8fc-4bac7229552c.svg)|

| Phase portrait vector field of the boundary value problem |Solutions of the Boundary layer problem in the phase portrait |
| :---: | :---: |
|  ![Phase portrait of the boundary value problem](https://user-images.githubusercontent.com/56647779/236527908-bc6da3b7-7807-40cf-b364-1657b0ba4ca8.svg) | ![Solutions of the Boundary layer problem in the phase portrait](https://user-images.githubusercontent.com/56647779/236527766-1e229407-41b3-4ba4-a357-82dee2015cf2.svg) |

For the solution to the boundary layer problem in the $x$ and $y$ plane and other figures, take a look into the slides!

## Usage

### Interactive Notebook
__Requirements:__
- `julia` (It is possible to run the notebook using )

>_Note: You have to activate the project and download the dependencies. It is as easy as running `update` in the `julia`
package REPL. You can find how to do this at the [`Pkg.jl`](https://pkgdocs.julialang.org/v1/) documentation page._

It is as easy as running 

```bash
$ julia --project  --eval "using Pluto; Pluto.run()"
```
in the [`Slides_Interactive`](Slides_Interactive/) directory and selecting `slides.jl` in the notebook UI. You can run it online in `binder` that is linked (top-right corner) from
static [`HTML` version](Slides_Interactive/slides.html) of the interactive slides.

The `GLMakie` figures that are more interactive can be run by re-evaluating the cells that define the figures.

### $\LaTeX$ Compilation
__Requirements:__ 
- `tectonic` (but any $\LaTeX$ compiler will do if you know what you are doing)
- `biber`

>_Note: Currently there is a problem with the version of `biblatex` included in the `tectonic` bundle, so if `biber`
fails to run, you should look over here: [`tectonic#35`](https://github.com/tectonic-typesetting/tectonic/issues/35)_

The slides were compiled by [`tectonic`](https://github.com/tectonic-typesetting/tectonic) using
```bash
$ tectonic -X build
```
ran in the [`Slides`](Slides/) directory. You can disable / enable presenter notes by modifying the `beamer` package
options in the [`preamble`](Slides/src/_preamble.tex)

| Notes | No notes|
| :---: | :---: |
| `\documentclass[x11names]{beamer}` | `\documentclass[x11names, notes]{beamer}` |

### Figures Generation
__Requirements:__
- `julia`
- `inkscape` (used for creating `PDF_TEX` from the `SVG` files)

They are generate using `CairoMakie` as `SVG` graphics and converted to `PDF_TEX` using `inkscape` by running the `julia` script

```bash
$ julia --project  --eval "include(\"generate_figures.jl\")"
```
inside the [`Slides/figures/sources`](Slides/figures/sources/) directory.

>_Note: You have to activate the project and download the dependencies. It is as easy as running `update` in the `julia`
package REPL. You can find how to do this at the [`Pkg.jl`](https://pkgdocs.julialang.org/v1/) documentation page._

## License

This work is triple-licensed under Apache 2.0 and GPL 2.0 (or any later version) and CC0 1.0 Universal .
You can choose between one of them if you use this work.

`SPDX-License-Identifier: Apache-2.0 OR GPL-2.0-or-later OR CC0-1.0-universal`
