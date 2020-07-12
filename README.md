# V-Polyhedral Disjunctive Cuts
#### Aleksandr M. Kazachkov
#### Based on joint work with Egon Balas

This project contains the code for one implementation of the V-polyhedral disjunctive cut paradigm.

## Installation

1. Clone the code using `git clone git@github.com:akazachk/vpc.git`.
2. Install `Cbc` by running `scripts/install_coin.sh`. If there are any problems, please start an issue (in this project, [coinbrew](https://github.com/coin-or/coinbrew), or [Cbc](https://github.com/coin-or/Cbc)).
3. [Optional] Install Gurobi or CPLEX.
4. Choose the appropriate options in the makefile under `Variables user should set`: `PROJ_DIR`, for the location of the repository, and `COIN_OR`, for where Cbc is installed. If you wish to use Gurobi, set `USE_GUROBI=1` under `Options for solvers`, and set `GUROBI_DIR` appropriately. Similarly, set `USE_CPLEX=1` and `CPLEX_DIR` if you wish to use CPLEX.
5. There are two compilation modes: `debug` and `release`. These can be compiled with `make [debug or release]`, which create the executable `vpc` in a new subdirectory `Debug` or `Release` of the main folder.
6. You can test the code with `make test_debug` or `make test_release`, depending on which version you compiled. Alternatively, run `test/run_test.sh`, which assumes the `Debug` version has been compiled.

## Dependencies

The code relies on [Cbc](https://github.com/coin-or/Cbc). It was extensively tested with version 2.9 (up to revision 2376 in the subversion history), while performance in the trunk and 2.10 versions has been more unstable. I have not yet tracked down the precise reason; see the related open GitHub issue #12. In particular, as of 2020/03/27, we need to compile `Cbc` with `SAVE_NODE_INFO` defined to enable access to the `parentNode` code in the `CbcModel` files.

Some of the scripts use `bash`, but can probably be adapted to other shells.

The user needs to export the shell variable `VPC_DIR`, pointing to the repository location, or define it in `scripts/install_coin.sh`, `test/run_test.sh`, and other scripts.

You may need to use a compatible version of `clang` or `g++`, and make sure that the same compiler is used when running `install_coin.sh` and the one used to generate `libgurobi_c++.a`, which can be rebuilt in the `Gurobi` directory under `src/build`.

You may need to install several dependencies, such as `libbz2-dev`, or comment out from `makefile` any linking to the relevant libraries you are missing.

## Execution

In more detail, you run the code with:

```
[Debug or Release]/vpc -f filename -d num_disj_terms
```

The flag `-f` is used to specify the filename. Then `-d` takes an integer (the number of disjunctive terms to generate, though we also use this for the number of split disjunctions). Other options are described when calling the code with `--help` or `-h`.

## Advanced Execution

To run experiments with a set of instances, use `scripts/run_experiments.sh`, which will put output in the `results` folder.
You will either need to export the shell variable `VPC_DIR`, pointing to the root directory of the repository, or enter it when prompted by the script.
The script takes three arguments:
1. the full path to an instance (in .mps or .lp format, possibly compressed with `gzip` or `bz2`) or a set of instances given in a file with extension `.instances` or `.batch`, with a specific format detailed below
2. the path to a directory where you wish to save results
3. the type of experiments to run; options are: "preprocess" (presolve instances), "bb0" (run each instance 7 times without cuts to gather baseline statistics), or "bb" (run each instance with VPCs from 2, 4, 8, 16, 32, and 64 term disjunctions)

If the first argument to `run_experiments.sh` is not an instance, it needs to be a file with extension `.instances` or `.batch`.
The former is for a set of instances that should be run one-by-one, and the latter is for sets of instances that should be run in parallel.
For `.instances` files, the lines must be _relative_ paths to instances, relative to the directory `data/instances`, such as `original/miplib2/bm23` (the extension can be left out from each line).
For `.batch` files, a new batch is indicated by a line that ends with a forward slash, e.g., `batchname/`.
Under each batch, instances should be indicated with relative paths just as in `.instances` files.
For examples, see `scripts/test.instances` and `scripts/test.batch`.

If you run experiments in batch mode, or generally save results in `/path/to/results/*/vpc-{type}.csv`, where `*` is a set of folders, then you can call
```
scripts/merge.sh /path/to/results type
```
to merge (and sort) the results to `/path/to/results/vpc-{type}.csv`. Note that `preprocess` is a special type, for which the script will merge results in `cleaning_log.csv` instead of `vpc-{type}.csv`.

## Details

There are many parameters that can be set from the command line. Run with `-h` or `--help` option to see these parameters. There are several more that are not currently able to be set from the command line, as they are assumed to generally be "constants"; a description of these can be found in `VPCParameters.hpp`.

The key parts of the code are the classes `Disjunction` (an abstract class), `CglVPC`, and `PRLP`. The user must specify a disjunction to `CglVPC` and a set of parameters (given as the struct `VPCParameters`). The disjunction should inherit from `Disjunction` and provide a set of disjunctive terms (the vector `terms`) and their number (`num_terms`). Each term is a `DisjunctiveTerm` defined in `Disjunction.hpp` and contains the following members: 
1. `CoinWarmStart* basis`, the optimal basis for that term
2. `obj`, the objective for that term (used in `CglVPC` for ensuring the term is correctly reconstructed)
3. `changed_var`, `changed_bound`, `changed_value`, `ineqs` vectors, which state which bounds need to be changed or inequalities to be added to define the disjunction.

In addition, a `Disjunction` has a few optional members: 
1. `name`, a name for the disjunction
2. the set of bounds changed or inequalities added at the root node that apply to all disjunctive terms
3. `integer_sol` and `integer_obj`, an integer-feasible solution and its objective value
4. `timer`, a way to do timing via `VPCTimeStats`

Before running `generateCuts` from `CglVPC`, ensure that the cut limit is properly set (especially when cuts are done in rounds or multiple disjunctions are used per round, the user must decide how to set the cut limit before each call to `generateCuts`).

## Future

There are many things left to be implemented in the future:
1. Cuts that do not cut off the LP optimum (this capability will be controlled by the parameter `PRLP_beta`).
2. Objective functions based on intermediate nodes of the partial branch-and-bound tree, to cut off points other than the LP optimum.
3. Dynamic disjunctions, in which cutting and branching are alternated.
4. Strengthening cuts.

## To install `Cbc`
You may need `gfortran`, `pkg-conf`, `LAPACK`, `BLAS`.

## Contact Information
Aleksandr M. Kazachkov,
akazachk AT discreteopt DOT com
