# Neighbor-Net implementation audit

SplitsPy 0.0.10's circular ordering was repeatable on the main 17-individual matrix. Its bundled constrained-weight solver was not: five repeated calls produced 45/55 positive splits and a maximum reconstructed-distance change of 0.006760057871787628. This was caught before manuscript integration. Its `distances.ls_fit` also initializes a residual accumulator to zero and multiplies it, reporting a spurious perfect fit; it was never used in reported analyses.

Final implementation uses the original authors' Neighbor-Net circular ordering, explicitly enumerates all 136 circular splits for 17 taxa, and fits nonnegative weights with SciPy `optimize.nnls`. KKT optimality, repeatability and additive-quartet tests are checked. Fitted-distance residuals are computed independently. The package is not patched. The phylogenetic outline only draws final weighted splits; the complete split system is retained as NEXUS/TSV.

An automatic phangorn Neighbor-Net is also retained. Different circular ordering can yield small alternative splits, so automatic ordering agreement and same-order weight verification are distinguished. Tree-based locality results are independently reproduced by Biopython and ape.
