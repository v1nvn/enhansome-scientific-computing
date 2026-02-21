# Awesome Scientific Computing [![Awesome](https://awesome.re/badge.svg)](https://awesome.re) with stars

[<img src="https://nschloe.github.io/awesome-scientific-computing/sunglasses.svg" align="right" width="30%">](#readme)

> Useful resources for scientific computing and numerical analysis.

Scientific computing and numerical analysis are research fields that aim to provide
methods for solving large-scale problems from various areas of science with the help of
computers. Typical problems are ordinary and partial differential equations (ODEs,
PDEs), their discretizations, and the solution of linear algebra problems arising from
them.

## Contents

* [Basic linear algebra](#basic-linear-algebra)
* [Multi-purpose toolkits](#multi-purpose-toolkits)
* [Finite Elements](#finite-elements)
* [Meshing](#meshing)
* [Data formats](#data-formats)
* [Sparse linear solvers](#sparse-linear-solvers)
* [Visualization](#visualization)
* [Other libraries and tools](#other-libraries-and-tools)
* [Community](#community)

## Basic linear algebra

* [OpenBLAS](https://www.openblas.net) - Optimized BLAS library based on GotoBLAS2.
  (C and Assembly, BSD, [GitHub](https://github.com/OpenMathLib/OpenBLAS) ⭐ 7,290 | 🐛 125 | 🌐 C | 📅 2026-02-20)
* [BLIS](https://github.com/flame/blis) ⭐ 2,610 | 🐛 121 | 🌐 C | 📅 2025-11-11 - High-performance BLAS-like dense linear algebra libraries.
  (C, BSD, GitHub)
* [BLAS](https://netlib.org/blas/) - Standard building blocks for performing basic vector and matrix operations.
  (Fortran, public domain, [GitHub](https://github.com/Reference-LAPACK/lapack/tree/master/BLAS) ⭐ 1,808 | 🐛 176 | 🌐 Fortran | 📅 2026-01-20)
* [LAPACK](https://netlib.org/lapack/) - Routines for solving systems of linear equations, linear least-squares, eigenvalue problems, etc.
  (Fortran, BSD, [GitHub](https://github.com/Reference-LAPACK/lapack) ⭐ 1,808 | 🐛 176 | 🌐 Fortran | 📅 2026-01-20)
* [Ginkgo](https://ginkgo-project.github.io/) - High-performance manycore linear algebra library, focus on sparse systems.
  (C++, BSD, [GitHub](https://github.com/ginkgo-project/ginkgo) ⭐ 559 | 🐛 182 | 🌐 C++ | 📅 2026-02-20)
* [Eigen](https://eigen.tuxfamily.org/index.php?title=Main_Page) - C++ template library for linear algebra.
  (C++, MPL 2, [GitLab](https://gitlab.com/libeigen/eigen))
* [blaze](https://bitbucket.org/blaze-lib/blaze) - High-performance C++ math library for dense and sparse arithmetic.
  (C++, BSD, Bitbucket)

## Multi-purpose toolkits

* [NumPy](https://numpy.org/) - Fundamental package needed for scientific computing with Python.
  (Python, BSD, [GitHub](https://github.com/numpy/numpy) ⭐ 31,479 | 🐛 2,311 | 🌐 Python | 📅 2026-02-20)
* [SciPy](https://scipy.org) - Python modules for statistics, optimization, integration, linear algebra, etc.
  (Python, mostly BSD, [GitHub](https://github.com/scipy/scipy/) ⭐ 14,476 | 🐛 1,771 | 🌐 Python | 📅 2026-02-20)
* [DifferentialEquations.jl](https://docs.sciml.ai/DiffEqDocs/stable/) - Toolbox for solving different types of differential equations numerically. (Julia, MIT, [GitHub](https://github.com/SciML/DifferentialEquations.jl) ⭐ 3,059 | 🐛 168 | 🌐 Julia | 📅 2026-02-08)
* [PETSc](https://petsc.org/release/) - Parallel solution of scientific applications modeled by PDEs.
  (C, 2-clause BSD, [GitLab](https://gitlab.com/petsc/petsc))
* [DUNE Numerics](https://www.dune-project.org) - Toolbox for solving PDEs with grid-based methods.
  (C++, GPL 2, [GitLab](https://gitlab.dune-project.org/core/))

## Finite Elements

* [MOOSE](https://mooseframework.inl.gov/) - Multiphysics Object Oriented Simulation Environment.
  (C++, LGPL 2.1, [GitHub](https://github.com/idaholab/moose) ⭐ 2,157 | 🐛 2,629 | 🌐 C++ | 📅 2026-02-21)
* [MFEM](https://mfem.org) - Free, lightweight, scalable C++ library for finite element methods.
  (C++, BSD-3-Clause, [GitHub](https://github.com/mfem/mfem) ⭐ 2,104 | 🐛 234 | 🌐 C++ | 📅 2026-02-21)
* [deal.II](https://dealii.org) - Software library supporting the creation of finite element codes.
  (C++, LGPL 2.1, [GitHub](https://github.com/dealii/dealii) ⭐ 1,624 | 🐛 646 | 🌐 C++ | 📅 2026-02-20)
* [SfePy](https://sfepy.org) - Simple Finite Elements in Python.
  (Python, BSD, [GitHub](https://github.com/sfepy/sfepy) ⭐ 824 | 🐛 76 | 🌐 Python | 📅 2026-02-11)
* [libMesh](https://libmesh.github.io) - Framework for the numerical simulation of PDEs using unstructured discretizations.
  (C++, LGPL 2.1, [GitHub](https://github.com/libMesh/libmesh) ⭐ 736 | 🐛 322 | 🌐 C | 📅 2026-02-20)
* [Firedrake](https://www.firedrakeproject.org) - Automated system for the solution of PDEs using the finite element method.
  (Python, LGPL 3, [GitHub](https://github.com/firedrakeproject/firedrake) ⭐ 632 | 🐛 366 | 🌐 Python | 📅 2026-02-20)
* [scikit-fem](https://github.com/kinnala/scikit-fem) ⭐ 608 | 🐛 11 | 🌐 Python | 📅 2026-01-30 - Simple finite element assemblers.
  (Python, BSD/GPL, GitHub)
* [Netgen/NGSolve](https://ngsolve.org) - High performance multiphysics finite element software.
  (C++, LGPL 2.1, [GitHub](https://github.com/NGSolve/netgen) ⭐ 366 | 🐛 110 | 🌐 C++ | 📅 2026-02-16)
* [libceed](https://libceed.readthedocs.io/en/latest/index.html) - Code for Efficient Extensible Discretizations.
  (C, 2-clause BSD, [GitHub](https://github.com/CEED/libCEED) ⭐ 245 | 🐛 54 | 🌐 C | 📅 2026-02-20)
* [FEniCS](https://fenicsproject.org) - Computing platform for solving PDEs in Python and C++.
  (C++/Python, LGPL 3, [GitHub](https://github.com/FEniCS)/[Bitbucket](https://bitbucket.org/fenics-project/))
* [FreeFEM](https://freefem.org) - High level multiphysics-multimesh finite element language.
  (C++, LGPL, [GitHub](https://github.com/FreeFem))

## Meshing

### Triangular and tetrahedral meshing

* [CGAL](https://www.cgal.org) - Algorithms for computational geometry.
  (C++, mixed LGPL/GPL, [GitHub](https://github.com/CGAL/cgal) ⭐ 5,757 | 🐛 662 | 🌐 C++ | 📅 2026-02-18)
* [trimesh](https://trimesh.org) - Loading and using triangular meshes with an emphasis on watertight surfaces.
  (Python, MIT, [GitHub](https://github.com/mikedh/trimesh) ⭐ 3,488 | 🐛 474 | 🌐 Python | 📅 2026-02-10)
* [pygmsh](https://github.com/nschloe/pygmsh) ⭐ 951 | 🐛 60 | 🌐 Python | 📅 2023-10-04 - Python interface for Gmsh.
  (Python, GPL 3, GitHub)
* [TetWild](https://arxiv.org/abs/1908.03581) - Generate tetrahedral meshes for triangular surface meshes.
  (C++, GPL 3, [GitHub](https://github.com/Yixin-Hu/TetWild) ⭐ 697 | 🐛 31 | 🌐 C++ | 📅 2023-04-27)
* [pygalmesh](https://github.com/meshpro/pygalmesh) ⭐ 665 | 🐛 33 | 🌐 C++ | 📅 2024-07-22 - Python interface for CGAL's 3D meshing capabilities.
  (Python, GPL 3, GitHub)
* [MeshPy](https://mathema.tician.de/software/meshpy/) - Quality triangular and tetrahedral mesh generation.
  (Python, MIT, [GitHub](https://github.com/inducer/meshpy) ⭐ 572 | 🐛 19 | 🌐 C++ | 📅 2026-01-12)
* [fTetWild](https://arxiv.org/abs/1908.03581) - Same as TetWild, but faster.
  (C++, MPL 2, [GitHub](https://github.com/wildmeshing/fTetWild) ⭐ 543 | 🐛 31 | 🌐 C++ | 📅 2025-11-13)
* [TriWild](https://cims.nyu.edu/gcl/papers/2019-TriWild.pdf) - Robust triangulation with curve constraints.
  (C++, MPL 2, [GitHub](https://github.com/wildmeshing/TriWild) ⭐ 260 | 🐛 6 | 🌐 C++ | 📅 2025-12-22)
* [dmsh](https://github.com/meshpro/dmsh) ⭐ 222 | 🐛 5 | 📅 2023-03-13 - Simple generator for unstructured triangular meshes, inspired by distmesh.
  (Python, proprietary, GitHub)
* [SeismicMesh](https://github.com/krober10nd/SeismicMesh) ⚠️ Archived - Parallel 2D/3D triangle/tetrahedral mesh generation with sliver removal.
  (Python and C++, GPL 3, GitHub)
* [Gmsh](https://gmsh.info) - Three-dimensional finite element mesh generator with pre- and post-processing facilities.
  (C++, GPL, [GitLab](https://gitlab.onelab.info/gmsh/gmsh))
* [TetGen](https://www.wias-berlin.de/software/index.jsp?id=TetGen) - Quality tetrahedral mesh generator and 3D Delaunay triangulator.
  (C++, AGPLv3)
* [Triangle](https://www.cs.cmu.edu/~quake/triangle.html) - Two-dimensional quality mesh generator and Delaunay triangulator.
  (C, *nonfree software*)
* [distmesh](https://persson.berkeley.edu/distmesh/) - Simple generator for unstructured triangular and tetrahedral meshes.
  (MATLAB, GPL 3)

### Quadrilateral and hexahedral meshing

* [QuadriFlow](https://stanford.edu/~jingweih/papers/quadriflow/) - Scalable and robust quadrangulation from triangulation.
  (C++, BSD, [GitHub](https://github.com/hjwdzh/QuadriFlow) ⭐ 805 | 🐛 8 | 🌐 C++ | 📅 2019-12-07)

### Mesh tools

* [meshio](https://github.com/nschloe/meshio) ⭐ 2,261 | 🐛 244 | 🌐 Python | 📅 2024-07-23 - I/O for various mesh formats, file conversion.
  (Python, MIT, GitHub)
* [pmp-library](https://www.pmp-library.org/) - Polygon mesh processing library.
  (C++, MIT with Employer Disclaimer, [GitHub](https://github.com/pmp-library/pmp-library/) ⭐ 1,450 | 🐛 17 | 🌐 C++ | 📅 2026-02-19)
* [optimesh](https://github.com/meshpro/optimesh) ⭐ 626 | 🐛 3 | 📅 2024-01-21 - Triangular mesh smoothing.
  (Python, proprietary, GitHub)
* [Mmg](https://www.mmgtools.org/) - Robust, open-source & multidisciplinary software for remeshing.
  (C, LGPL 3, [GitHub](https://github.com/MmgTools/mmg) ⭐ 451 | 🐛 41 | 🌐 C | 📅 2026-02-04)
* [meshplex](https://github.com/meshpro/meshplex) ⭐ 109 | 🐛 8 | 📅 2025-07-23 - Fast tools for simplex meshes.
  (Python, proprietary, GitHub)
* [MOAB](https://sigma.mcs.anl.gov/moab-library/) - Representing and evaluating mesh data.
  (C++, mostly LGPL 3, [Bitbucket](https://bitbucket.org/fathomteam/moab/))

## Data formats

* [Zarr](https://zarr.readthedocs.io/en/stable/) - Format for the storage of chunked, compressed, N-dimensional arrays.
  (Python, MIT, [GitHub](https://github.com/zarr-developers/zarr-python) ⭐ 1,919 | 🐛 549 | 🌐 Python | 📅 2026-02-20)
* [HDF5](https://www.hdfgroup.org/solutions/hdf5/) - Data model, library, and file format for storing and managing data.
  (C/Fortran, BSD, [GitHub](https://github.com/HDFGroup/hdf5) ⭐ 895 | 🐛 307 | 🌐 C | 📅 2026-02-21)
* [NetCDF](https://www.unidata.ucar.edu/software/netcdf) - Software libraries and data formats for array-oriented scientific data.
  (C/C++/Fortran/Java/Python, [custom open-source
  license](https://www.unidata.ucar.edu/software/netcdf/licensing),
  [GitHub](https://github.com/Unidata/netcdf-c/) ⭐ 582 | 🐛 305 | 🌐 C | 📅 2026-02-18)
* [XDMF](https://xdmf.org/) - eXtensible Data Model and Format for data from High Performance Computing codes.
  (C++, [GitLab](https://gitlab.kitware.com/xdmf/xdmf))

## Sparse linear solvers

* [hypre](https://computing.llnl.gov/projects/hypre-scalable-linear-solvers-multigrid-methods) - Library of high-performance preconditioners and solvers.
  (C, Apache 2.0/MIT, [GitHub](https://github.com/hypre-space/hypre) ⭐ 820 | 🐛 162 | 🌐 C | 📅 2026-02-18)
* [PyAMG](https://pyamg.readthedocs.io/en/latest/) - Algebraic Multigrid Solvers in Python.
  (Python, MIT, [GitHub](https://github.com/pyamg/pyamg) ⭐ 638 | 🐛 36 | 🌐 Python | 📅 2026-01-02)
* [SuperLU](https://portal.nersc.gov/project/sparse/superlu/) - Direct solution of large, sparse, nonsymmetric systems of linear equations.
  (C, mostly BSD, [GitHub](https://github.com/xiaoyeli/superlu) ⭐ 326 | 🐛 27 | 🌐 C | 📅 2026-02-06)

## Visualization

* [F3D](https://f3d.app/) - Cross-platform, fast, and minimalist 3D viewer with scientific visualization tools.
  (C++, BSD, [GitHub](https://github.com/f3d-app/f3d) ⭐ 4,136 | 🐛 270 | 🌐 C++ | 📅 2026-02-21)
* [PyVista](https://docs.pyvista.org/) - 3D plotting and mesh analysis through a streamlined interface for VTK.
  (Python, MIT, [GitHub](https://github.com/pyvista/pyvista) ⭐ 3,530 | 🐛 740 | 🌐 Python | 📅 2026-02-21)
* [vedo](https://vedo.embl.es) - Library for scientific analysis and visualization of 3D objects based on VTK.
  (Python, MIT, [GitHub](https://github.com/marcomusy/vedo) ⭐ 2,237 | 🐛 164 | 🌐 Python | 📅 2026-02-18)
* [Polyscope](https://polyscope.run/) - Viewer and user interface for 3D geometry processing.
  (C++, MIT, [GitHub](https://github.com/nmwsharp/polyscope) ⭐ 2,141 | 🐛 118 | 🌐 C++ | 📅 2026-02-18)
* [Mayavi](https://docs.enthought.com/mayavi/mayavi/) - 3D scientific data visualization and plotting in Python.
  (Python, BSD, [GitHub](https://github.com/enthought/mayavi) ⭐ 1,392 | 🐛 483 | 🌐 Python | 📅 2025-10-01)
* [yt](https://yt-project.org/) - Toolkit for analysis and visualization of volumetric data.
  (Python, BSD, [GitHub](https://github.com/yt-project/yt) ⭐ 540 | 🐛 473 | 🌐 Python | 📅 2026-02-20)
* [TTK](https://topology-tool-kit.github.io/) - Topological data analysis and visualization.
  (C++/Python, BSD, [GitHub](https://github.com/topology-tool-kit/ttk) ⭐ 466 | 🐛 28 | 🌐 C++ | 📅 2026-01-30)
* [morphologica](https://github.com/ABRG-Models/morphologica) ⚠️ Archived - Header-only, modern OpenGL code to visualize numerical simulations at runtime. (C++, Apache 2.0, GitHub)
* [ParaView](https://www.paraview.org) - Multi-platform data analysis and visualization application based on VTK.
  (C++, BSD, [GitLab](https://gitlab.kitware.com/paraview/paraview))
* [VTK](https://vtk.org/) - Process images and create 3D computer graphics.
  (C++, BSD, [GitLab](https://gitlab.kitware.com/vtk/vtk))

## Other libraries and tools

* [cvxpy](https://www.cvxpy.org/) - Modeling language for convex optimization problems.
  (Python, Apache 2.0, [GitHub](https://github.com/cvxpy/cvxpy) ⭐ 6,112 | 🐛 265 | 🌐 C++ | 📅 2026-02-19)
* [FFTW](http://www.fftw.org) - Discrete Fourier transforms in one or more dimensions, of arbitrary input size, real and complex.
  (C, GPL2, [GitHub](https://github.com/FFTW/fftw3) ⭐ 3,023 | 🐛 177 | 🌐 C | 📅 2026-02-16)
* [PyWavelets](https://pywavelets.readthedocs.io/en/latest/) - Wavelet transforms in Python.
  (Python, MIT, [GitHub](https://github.com/PyWavelets/pywt) ⭐ 2,345 | 🐛 81 | 🌐 Python | 📅 2026-02-01)
* [OpenFOAM](https://www.openfoam.com) - Free, open source CFD (computational fluid dynamics) software.
  (C++, GPL 3, [GitHub](https://github.com/OpenFOAM/OpenFOAM-dev) ⭐ 1,974 | 🐛 7 | 🌐 C++ | 📅 2026-02-20)
* [pyGAM](https://pygam.readthedocs.io/en/latest/) - Generalized Additive Models in Python.
  (Python, Apache 2.0, [GitHub](https://github.com/dswah/pyGAM) ⭐ 957 | 🐛 112 | 🌐 Python | 📅 2026-01-22)
* [Qhull](http://www.qhull.org) - Convex hull, Delaunay triangulation, Voronoi diagram, halfspace intersection about a point, etc.
  (C/C++, [custom open source license](http://www.qhull.org/COPYING.txt),
  [GitHub](https://github.com/qhull/qhull/) ⭐ 808 | 🐛 29 | 🌐 C | 📅 2025-09-07)
* [quadpy](https://github.com/sigma-py/quadpy) ⭐ 785 | 🐛 31 | 📅 2023-04-04 - Numerical integration (quadrature, cubature) in Python.
  (Python, proprietary, GitHub)
* [Chebfun](https://www.chebfun.org/) - Computing with functions to about 15-digit accuracy.
  (MATLAB, BSD, [GitHub](https://github.com/chebfun/chebfun) ⭐ 660 | 🐛 171 | 🌐 MATLAB | 📅 2025-09-27)
* [Dedalus](https://dedalus-project.org/) - Solve partial differential equations with spectral methods.
  (Python, GPL 3, [GitHub](https://github.com/DedalusProject/dedalus) ⭐ 658 | 🐛 61 | 🌐 Python | 📅 2026-01-21)
* [FiPy](https://www.ctcms.nist.gov/fipy/) - Finite-volume PDE solver.
  (Python, [custom open-source
  license](https://www.nist.gov/open/copyright-fair-use-and-licensing-statements-srd-data-software-and-technical-series-publications),
  [GitHub](https://github.com/usnistgov/fipy) ⭐ 600 | 🐛 168 | 🌐 Python | 📅 2026-02-20)
* [PyGMO](https://esa.github.io/pygmo/) - Massively parallel optimization.
  (Python/C++, MPL 2, [GitHub](https://github.com/esa/pygmo2) ⭐ 520 | 🐛 49 | 🌐 C++ | 📅 2024-08-10)
* [pyMOR](https://pymor.org/) - Model Order Reduction with Python.
  (Python, 2-clause BSD, [GitHub](https://github.com/pymor/pymor/) ⭐ 340 | 🐛 123 | 🌐 Python | 📅 2026-02-19)
* [shenfun](https://shenfun.readthedocs.io/en/latest/) - High-performance Python library for the spectral Galerkin method.
  (Python, BSD-2, [GitHub](https://github.com/spectralDNS/shenfun) ⭐ 226 | 🐛 29 | 🌐 Python | 📅 2025-12-16)
* [orthopy](https://github.com/sigma-py/orthopy) ⭐ 189 | 🐛 8 | 📅 2024-02-15 - Compute orthogonal polynomials efficiently.
  (Python, proprietary, GitHub)
* [NFFT](https://www-user.tu-chemnitz.de/~potts/nfft/) - Nonequispaced fast Fourier transform.
  (C/MATLAB, GPL 2, [GitHub](https://github.com/NFFT/nfft) ⭐ 184 | 🐛 22 | 🌐 C | 📅 2026-02-01)
* [HPDDM](https://github.com/hpddm/hpddm) ⭐ 154 | 🐛 1 | 🌐 C++ | 📅 2026-02-19 - High-performance unified framework for domain decomposition methods.
  (C++, LGPL 3, GitHub)
* [PyDMD](https://github.com/mathLab/PyDMD) ⭐ 117 | 🐛 0 | 🌐 Python | 📅 2025-03-06 - Dynamic Mode Decomposition (DMD) in Python.
  (Python, MIT, GitHub)
* [accupy](https://github.com/sigma-py/accupy) ⭐ 107 | 🐛 2 | 🌐 Python | 📅 2021-09-06 - Accurate sums and dot products for Python.
  (Python, GPL 3, GitHub)
* [GSL](https://www.gnu.org/software/gsl/) - Random number generators, special functions, and least-squares fitting etc.
  (C/C++, GPL 3, [Savannah](https://savannah.gnu.org/projects/gsl))
* [SLEPc](https://slepc.upv.es) - Scalable Library for Eigenvalue Problem Computations.
  (C, 2-clause BSD, [GitLab](https://gitlab.com/slepc/slepc))
* [preCICE](https://precice.org/) - Coupling library for partitioned multi-physics simulations (FSI, CHT, and more).
  (C++, LGPL 3, [GitHub](https://github.com/precice/))

## Community

* [SciComp StackExchange](https://scicomp.stackexchange.com/) - Computational Science on the StackExchange network.
* [Wolfgang Bangerth's video class](https://www.math.colostate.edu/~bangerth/videos.html) - MATH 676: Finite element methods in scientific computing.
* [Nick Higham's blog](https://nhigham.com/) - Mostly on MATLAB, general computing advice.
* [Nick Trefethen's Video Lectures](https://people.maths.ox.ac.uk/trefethen/videos.html) - 36 video lectures on approximation theory/practice and scientific computing.
* [John D. Cook's blog](https://www.johndcook.com/blog/) - Feats of scientific computing.
* [Jack Dongarra's software list](https://netlib.org/utk/people/JackDongarra/la-sw.html) - List of freely available software for the solution of linear algebra problems.
* [NA Digest](https://netlib.org/na-digest-html/) - Collection of articles on topics related to numerical analysis and those who practice it.
* [Gabriel Peyré on Bluesky](https://bsky.app/profile/gabrielpeyre.bsky.social) - One post a day on computational mathematics.
* [Discord: Numerical Software](https://discord.com/invite/hnTJ5MRX2Y) - Discord messaging server on numerical software.
