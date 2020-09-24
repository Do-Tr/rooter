# Rooter
Rooter is a Bayesian molecular clock rooting method described in (Huelsenbeck <i>et al</i>., 2002). Rooter calculates the posterior probability distribution of the root position under the molecular clock.

Huelsenbeck, J. P., Bollback, J. P., and Levine, A. M. 2002. Inferring  the  root  of  a  phylogenetic  tree. Systematic biology, 51(1): 32–43.
# Prerequisites
To compile this program, you will need to download the Eigen library, which you can find at http://eigen.tuxfamily.org/. Make certain to change the appropriate line in the Makefile, pointing your compilation to the Eigen install on your computer.
# Installation

```S
cd rooter
make
make install
```
# Command-line usage
Rooter can be used from the command line with arguments specifying input data and parameters. Only PHYLIP alignment format is supported. Trees can be read as newick. To run Rooter use:

```S
rooter --
```

