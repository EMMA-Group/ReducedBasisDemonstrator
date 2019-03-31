Reduced Basis Demonstrator
Example of the deformation gradient based Reduced Basis method for hyperelasticity

This software package is related to the research article
Authors: Oliver Kunc and Felix Fritzen
Title:   Finite strain homogenization using a reduced basis and efficient sampling
Journal: Mathematical and Computational Applications
         Special Issue "Machine Learning, Low-Rank Approximations and Reduced Order Modeling in Computational Mechanics"
Year:    2019
Volume:  24
Number:  1
DOI   ...
URL   dx.doi.org/...

Start the program by running ReducedBasisDemonstrator
See the paper for more information or contact the authors.

---

Explanatory and general notes

The Reduced Basis (RB) method approximates the deformation gradient F.
The RB matrix B contains the RB elements as columns. All matrices are
stored in row-major format (in contrast to MATLAB/Octave default, hence
some transpositions near the reshape commands). The columns of the RB
matrix B contain deformation gradient fluctuations (9 components) at each
quadrature point (total number: N_qp) for each basis element (total
number: N), thus the size(B) = [ 9*N_qp , N ]
All main steps contain in-line references to the corresponding formulas
in the open access paper.
The example data contains multiple sets of approximatively uniformly
distributed directions in 5 dimensions, the RB of an example micro-
structure, and the quadrature weights. The microstructure is a the unit
cube with an off-center cubical inclusion. Due to the nature of the
method, no mesh and not even point coordinates are necessary. For the
sake of shortness, this information is not contained in this repository.
Since the original mesh is very simple, all quadrature weights are equal
in this example. The code, however, is general and does not exploit this.

For bugs, comments, or suggestions please contact
Oliver Kunc: kunc@mechbau.uni-stuttgart.de
