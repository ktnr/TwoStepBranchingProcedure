# Nested tree-search for two-dimensional orthogonal packing problems [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.6673875.svg)](https://doi.org/10.5281/zenodo.6673875)

This repository is an attempt to reproduce the nested branch-and-bound algorithm by Clautiaux, Carlier and Moukrim (2007) to solve two-dimensional orthogonal packing problems (2D-OPP). 

The 2D-OPP asks, given a bin and non-rotatable items, can the items be placed in the bin without overlapping? The 2D-OPP is an important subproblem in many combinatorial optimization problems and, despite its simple formulation, is notoriously difficult to solve. There are instances with just 20-30 items that cannot be proven to be feasible or infeasible with current methods.

# Algorithmic outline

The outer branch-and-bound algorithm is called two-step branching procedure (TSBP) and is a relaxation of the 2D-OPP, where the goal is to find feasible x-coordinate placements. The relaxation simplifies the problem such that items need not be placed contiguously along the y-axis. For each feasible set of x-coordinate placements, the x-coordiantes are fixed and the inner branch-and-bound algorithm is called to search for corresponding feasible y-coordinate placements, i.e., a solution to the original 2D-OPP. The inner method is called left-most active only (LMAO) algorithm. A search tree is built that, at each node, spawns a new child node for each remaining unplaced item in addition to a deactivation node. Each normal node corresponds to a placement point and an item, whereas the deactivation node corresponds only to a placement point. 

The first placement point is the origin of the bin: (0, 0). Then, items are placed successively according to the branching order at the left-most downward *active* placement point. If a placement is feasible, the item is placed and the new child nodes for the remaining items are generated by determining the new active position. If a placement is infeasible, the current node is pruned. 
The deactivation node deactivates the current active placement point and determines the new active placement point in order to spawn new child nodes. The deactivation procedure interacts with the rest of the algorithm such that (almost) no equivalent placements are enumerated, which drastically reduces the size of the search tree.

The left-most active only algorithm can also be called as an independent solution algorithm but usually performs worse.

# Improvements

The original paper describes two symmetry reductions in the outer branch-and-bound tree: pseudo symmetry and block equivalence. These procedures are rather complex and are not implemented. Instead, a different, simple but effective symmetry reduction technique is implemented, which covers the same symmetry reductions as pseudo-symmetry and many of block equivalence. In fact, it covers all pseudo-symmetries and all pair-wise block equivalence reductions.

The branch-and-bound algorithm implicitly enumerates all normal patterns, cp. Herz (1972) resp. Beasely (1985). However, the placements still contain rotational symmetries. A simple way to exclude those symmetries is by applying a domain reduction according to Soh et al. (2010) to a single item, i.e., restricting the feasible placement points of a single item to the lower left corner of the bin. 
This is valid because any feasible solution can be transformed into a solution where any single item is within its reduced domain. The proof is analogous to the one of preprocessing step 1 in C??t?? and Iori (2018). 
Thereby, domain reduction covers all pair-wise block equivalence reductions. Merely block equivalences with three or more blocks, where blocks that do not contain the domain reduced item can be swapped, are not covered.
Domain reduction is also applicable to the inner branch-and-bound.

Additional improvements over the original algorithm include
- parallelization at the root node
- use of modern memory allocation

# Usage

Run the CLI-app with the following arguments, cp. [`.vscode/launch.json`](https://github.com/ktnr/TwoStepBranchingProcedure/blob/2fa7152b532571ff98a6fe40e790e19998e4d716/.vscode/launch.json):
```
Options:
  -h,--help                   Print this help message and exit
  -i,--inputdir TEXT REQUIRED The directory where the input file -f resides
  -f,--file TEXT REQUIRED     The input file name
  -o,--outdir TEXT REQUIRED   The output directory
  -p,--param TEXT             The .json parameter full file path
  -t,--timeSuffix BOOLEAN     If time should be appended to the output path (1=true, 0=false)

For example:
>tsbp.exe -i ${workspaceFolder}/data/input/CJCM08/ -f E00N10.json -o ${workspaceFolder}/data/output/
$tsbp -i ${workspaceFolder}/data/input/CJCM08/ -f E00N10.json -o ${workspaceFolder}/data/output/
```

Caution, the memory usage for the program increases rapidly. For hard instances which require a few billion nodes to solve, memory consumption can reach multiple hundreds of GB. To avoid excessive memory usage, we abort at 95% memory utilization.

https://github.com/ktnr/TwoStepBranchingProcedure/blob/2fa7152b532571ff98a6fe40e790e19998e4d716/tsbp/src/BranchAndBound.cpp#L2510-L2515

Reducing the memory utilization by releasing already explored parts of the search tree (https://github.com/ktnr/TwoStepBranchingProcedure/issues/11) will allow to explore more nodes but will have limited effect on the runtime.

# Evaluation

Hard instances generally exhibit low area waste ($\epsilon$) between 2% and 7%, see the original paper. Particularly difficult are instances where all items are smaller than half the container dimensions, see [ktnr/BinPacking2D/BPP-Subproblems](https://github.com/ktnr/BinPacking2D/tree/master/data/input/OPP/BPP-Subproblems).

This repository cannot replicate the performance (number of explored nodes for infeasible instances) of the original paper, not even the performance of the simplest version without block equivalence and pseudo symmetry reductions. So, there might still be a bug in the code, inefficiencies, or invalid assumptions. If you happen to have knowledge about implementation details of the original paper, please do get in touch.

The constraint programming appraoch of Clautiaux et al. (2008) is often cited to achieve the best performance to date. However, we couldn't reproduce their performance, even with a comparable model using modern soft- and hardware. For details please refer to https://github.com/google/or-tools/discussions/3177.

# Future research and pointers

The 2D-OPP is an important subproblem to many combinatorial optimization problems. Solving 2D-OPPs more efficiently may further progress for a variety of other problems. For example, the 2D-OPP is a subproblem to 2D bin packing problems. Some comparatively small 2D bin packing instances that have been introduced more than 30 years ago remain open to this day. Some of them because the arising 2D-OPPs are too hard. For suggestions on how to improve the solution times see the [issues page](https://github.com/ktnr/TwoStepBranchingProcedure/issues) and in particular https://github.com/ktnr/TwoStepBranchingProcedure/issues/5 for a list of references that might be compatible with this algorithm or might even be better suited to solve hard 2D-OPPs.

An interactive, educational version of the basic nested branch-and-bound algorithm is available [here](https://github.com/mike-ghes/CS9544-2dOrthogonalPacking). 

# References

- Clautiaux, F., Carlier, J., & Moukrim, A. (2007). A new exact method for the two-dimensional orthogonal packing problem. European Journal of Operational Research, 183(3), 1196-1211.
- Martello, S., & Vigo, D. (1998). Exact solution of the two-dimensional finite bin packing problem. Management science, 44(3), 388-399.
- Soh, T., Inoue, K., Tamura, N., Banbara, M., & Nabeshima, H. (2010). A SAT-based method for solving the two-dimensional strip packing problem. Fundamenta Informaticae, 102(3-4), 467-487.
- Herz, J. C. (1972). Recursive computational procedure for two-dimensional stock cutting. IBM Journal of Research and Development, 16(5), 462-469.
- Beasley, J. E. (1985). Algorithms for unconstrained two-dimensional guillotine cutting. Journal of the Operational Research Society, 36(4), 297-306.
- C??t??, J. F., & Iori, M. (2018). The meet-in-the-middle principle for cutting and packing problems. INFORMS Journal on Computing, 30(4), 646-661.
- Clautiaux, F., Jouglet, A., Carlier, J., & Moukrim, A. (2008). A new constraint programming approach for the orthogonal packing problem. Computers & Operations Research, 35(3), 944-959.
