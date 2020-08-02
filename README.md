# V-Polyhedral Disjunctive Cuts
#### [Aleksandr M. Kazachkov](https://akazachk.github.io)
#### Based on joint work with Egon Balas

This project contains the code for one implementation of the V-polyhedral disjunctive cut paradigm.

## Installation

1. Clone the code using `git clone git@github.com:akazachk/vpc.git`.
2. Install [Cbc](https://github.com/coin-or/Cbc) by running [`scripts/install_coin.sh`](scripts/install_coin.sh). If there are any problems, please start an issue (in this project, [coinbrew](https://github.com/coin-or/coinbrew), or [Cbc](https://github.com/coin-or/Cbc)). This step can be customized (e.g., if Cbc is already installed), but it is crucial that the latest Cbc version is compiled with the macro `SAVE_NODE_INFO` defined.
3. [Optional] Install Gurobi or CPLEX.
4. Choose the appropriate options in the [`makefile`](makefile) under `Variables user should set`: `PROJ_DIR`, for the location of the repository, and `COIN_OR`, for where Cbc is installed. If you wish to use Gurobi, set `USE_GUROBI=1` under `Options for solvers`, and set `GUROBI_DIR` appropriately. Similarly, set `USE_CPLEX=1` and `CPLEX_DIR` if you wish to use CPLEX.
5. There are two compilation modes: `debug` and `release`. These can be compiled with `make [debug or release]`, which creates the executable `vpc` in a new subdirectory `Debug` or `Release` of the main folder.
6. You can test the code with `make test_debug` or `make test_release`, depending on which version you compiled. Alternatively, run [`test/run_test.sh`](test/run_test.sh), which assumes the `Debug` version has been compiled.

## Dependencies

The VPC code relies on [Cbc](https://github.com/coin-or/Cbc). It was previously extensively tested with Cbc version 2.9 (up to revision 2376 in the subversion history), though this backwards compatibility is no longer actively maintained. Performance in the main/trunk and 2.10 versions has been more unstable. Specifically, Cbc version 2.10 will probably *not* work. See the related GitHub issue [#12](https://github.com/akazachk/vpc/issues/12). In particular, as of 2020/03/27, we need to compile the development branch of Cbc with `SAVE_NODE_INFO` defined (using `ADD_CXXFLAGS="-DSAVE_NODE_INFO"`) to enable access to the `parentNode` code in `CbcModel`. In addition, if the Cbc commit is before [0f6ffed](https://github.com/coin-or/Cbc/commit/0f6ffed4c26daaf75edac2f87b70f3cc40cb12fd), we need to comment out the line `currentNode_ = NULL` in `CbcModel.cpp` around [line CbcModel.cpp:15392](https://github.com/coin-or/Cbc/blob/53f34cfea21360091608b02a041a962b2be7d6bc/src/CbcModel.cpp#L15390-L15391).

For Cbc, you may need `gfortran`, `pkg-conf`, `LAPACK`, and `BLAS`. It may also be necessary to use `--with-cplex=false` as an option in the `coinbrew` commands, if [Osi](https://github.com/coin-or/Osi)'s configure script detects the CPLEX library through `CPXgetstat` but the CPLEX include directory is not found (see https://github.com/coin-or/coinbrew/issues/49).

Some of the scripts use `bash`, but can probably be adapted to other shells.

There shoud be an environment variable `VPC_DIR` pointing to the local repository location, or this variable can be defined in each of the scripts: [`scripts/install_coin.sh`](scripts/install_coin.sh), [`test/run_test.sh`](test/run_test.sh), and others.

You may need to use a compatible version of `clang` or `g++`, and make sure that the same compiler is used when running [`install_coin.sh`](scripts/install_coin.sh) and the one used to generate `libgurobi_c++.a`, which can be rebuilt in the Gurobi directory (located at, say, `${GUROBI_HOME}`), under `${GUROBI_HOME}/src/build`.

There are some optional libraries, such as `libbz2-dev`, that are linked to in [`makefile`](makefile). If these are missing, remove the corresponding linking in the makefile.

## Execution

In more detail, to run the code, a basic call will look like:

```
[Debug or Release]/vpc -f filename -d num_disj_terms
```

The flag `-f` is used to specify the filename. Then `-d` takes an integer (the number of disjunctive terms to generate, though we also use this for the number of split disjunctions). Other options are described when calling the code with `--help` or `-h`.

## Advanced Execution

To run experiments with a set of instances, use [`scripts/run_experiments.sh`](scripts/run_experiments.sh), which will by default put output in the `results` folder.
You will either need to export the shell variable `VPC_DIR`, pointing to the root directory of the repository, or enter it when prompted by the script.
The script takes four arguments:
1. the full path to an instance (in .mps or .lp format, possibly compressed with `gzip` or `bz2`) or a set of instances given in a file with extension `.instances` or `.batch`, with a specific format detailed below
2. the path to a directory where you wish to save results
3. the type of experiments to run; options are: "preprocess" (presolve instances), "bb0" (run each instance 7 times without cuts to gather baseline statistics), or "bb" (run each instance with VPCs from 2, 4, 8, 16, 32, and 64 term disjunctions)
4. [optional] extra arguments you wish to pass, in quotation marks, which will overwrite any defaults, and specifically it will not loop over various b&b tree depths, so the `-d` option needs to be set manually in this case

For example, you can call 
```
scripts/run_experiments.sh data/instances/original/miplib2/bm23.mps.gz results bb "-d 32 -t 120"
```
to run instance `bm23.mps`, send output to `results`, use up to a 32-term disjunction, and enforce a timelimit of 120 seconds (as opposed to the default of an hour).

If the first argument to [`run_experiments.sh`](scripts/run_experiments.sh) is not an instance, it needs to be a file with extension `.instances` or `.batch`.
The former is for a set of instances that should be run one-by-one, and the latter is for sets of instances that should be run in parallel.
For `.instances` files, the lines must either be full paths and end with .mps/.lp, such as `${VPC_DIR}/data/instances/original/miplib2/bm23.mps`, or be _relative_ paths to instances *without* the extension, relative to the directory `data/instances`, such as `original/miplib2/bm23` (the extension *must* be left out from each line in this case).
For `.batch` files, a new batch is indicated by a line that ends with a forward slash, e.g., `batchname/`.
Under each batch, instances should be indicated with relative paths just as in `.instances` files.
For examples, see [`scripts/test.instances`](scripts/test.instances) and [`scripts/test.batch`](scripts/test.batch).

If you run experiments in batch mode, or generally save results in `/path/to/results/*/vpc-{type}.csv`, where `*` is a set of folders, then you can call
```
scripts/merge.sh /path/to/results type
```
to merge (and sort) the results to `/path/to/results/vpc-{type}.csv`.
Note that `preprocess` is a special type, for which the script will merge results in `cleaning_log.csv` instead of `vpc-{type}.csv`.
The [`merge.sh`](scripts/merge.sh) script also takes a third optional argument for where to save the output.

## Details

There are many parameters that can be set from the command line. Run with `-h` or `--help` option to see these parameters. There are several more that are not currently able to be set from the command line, as they are assumed to generally be "constants"; a description of these can be found in [`VPCParameters.hpp`](include/VPCParameters.hpp).

The key parts of the code are the classes [`Disjunction`](include/Disjunction.hpp) (an abstract class), [`CglVPC`](include/CglVPC.hpp), and [`PRLP`](include/PRLP.hpp). The user must specify a disjunction to `CglVPC` and a set of parameters (given as the struct `VPCParameters`). The disjunction should inherit from `Disjunction` and provide a set of disjunctive terms (the vector `terms`) and their number (`num_terms`); see for example the [`SplitDisjunction`](include/SplitDisjunction.hpp) class. Each term is a `DisjunctiveTerm` defined in [`Disjunction.hpp`](include/Disjunction.hpp) and contains the following members: 
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

## Contact Information
[Aleksandr M. Kazachkov](https://akazachk.github.io),
[akazachk AT discreteopt DOT com](https://akazachk.github.io/email)
