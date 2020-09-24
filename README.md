# Rooter
Rooter is a Bayesian molecular clock rooting method described in (Huelsenbeck <i>et al</i>., 2002). Rooter calculates the posterior probability distribution of the root position under the molecular clock.

Huelsenbeck, J. P., Bollback, J. P., and Levine, A. M. 2002. Inferring  the  root  of  a  phylogenetic  tree. Systematic biology, 51(1): 32–43.
# Prerequisites
To compile this program, you will need to download the Eigen library, which you can find at http://eigen.tuxfamily.org/. Make certain to change the appropriate line in the Makefile (modify the INC definition), pointing your compilation to the Eigen install on your computer.
# Installation

```S
git clone https://github.com/NielsenBerkeleyLab/rooter.git
cd rooter
make
make install
```
# Command-line usage
Rooter can be used from the command line with arguments specifying input data and parameters. Only PHYLIP alignment format is supported. Only newick format is supported for tree files. For tip-dates <dates.tsv>, use a tab-separated format where the first column corresponds to the tip-name in the tree and the second column is the date with the format YYYY-MM-DD. To run Rooter use:

```S
rooter -i <input.phylip> -t <input.nwk> -c <dates.tsv> -l <number of MCMC iterations> -s <MCMC sampling frequency> -p <MCMC printing frequency> -o <output_prefix>
```
# Outputs
Rooter prints a log file <code><output_prefix>.log</code>with a summary of each root position of the credible set along with its posterior probability and cumulative probability of the root branch. It also prints out the trees <code><output_prefix>.t</code>for each bipartition in Nexus format.
