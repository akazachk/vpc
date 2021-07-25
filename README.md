# V-Polyhedral Disjunctive Cuts
#### [Aleksandr M. Kazachkov](https://akazachk.github.io)
#### Based on joint work with Egon Balas

This project contains the code for one implementation of the *V-polyhedral disjunctive cut paradigm*. Some more details can be found in Chapter 4 of my [dissertation](https://akazachk.github.io/pubs/dissertation.pdf), [talk slides from the INFORMS Annual Meeting in 2018](https://akazachk.github.io/pubs/INFORMS2018.pdf), [poster at NemFest 2017](https://akazachk.github.io/pubs/Nemfest2017.pdf), or video at a [2021 GERAD seminar](https://youtu.be/TyWemNwPGq0?t=4).

## Documentation

If [`doxygen`](https://www.doxygen.nl/index.html) is installed, type `make doxygen` for documentation on classes and functions in the code. This documentation is still incomplete.

## Installation

1. Check that the [dependencies](#dependencies) are present in your system. You may need to install [`homebrew`](https://brew.sh) on a Mac in order to download the dependencies, and installation on Windows has never been tested.
2. Clone the code by typing (in your terminal) `git clone git@github.com:akazachk/vpc.git` if you have a GitHub account and your SSH key has been added to it, or `git clone git://github.com/akazachk/vpc.git` otherwise. This will create a directory `vpc` and download the project code there.
3. Save this location via `export VPC_DIR=$PWD/vpc` (modified as needed if you change the directory you are located in, or the name of the `vpc` folder). Check that this worked via the output of `echo $VPC_DIR`. It will help (but is not necessary) to add the full path to the `vpc` directory via a definition of `VPC_DIR` in your `.bashrc` or `.zprofile` preference files in your home directory.
4. Install [Cbc](https://github.com/coin-or/Cbc) by running [`setup/install_coin.sh`](setup/install_coin.sh). If there are any problems, please start an issue (in this project, [coinbrew](https://github.com/coin-or/coinbrew), or [Cbc](https://github.com/coin-or/Cbc)). This step can be customized (e.g., if Cbc is already installed), but it is crucial that the latest Cbc version is compiled with the macro `SAVE_NODE_INFO` defined.
5. [Optional] Install Gurobi or CPLEX.
6. Choose the appropriate options in the [`Makefile`](Makefile) under `Variables user should set`: `PROJ_DIR`, for the location of the repository, and `COIN_OR`, for where Cbc is installed. If you wish to use Gurobi, set `USE_GUROBI=1` under `Options for solvers`, and set `GUROBI_DIR` appropriately. Similarly, set `USE_CPLEX=1` and `CPLEX_DIR` if you wish to use CPLEX. If you are not using Gurobi or CPLEX, make sure to set `USE_GUROBI=0` and `USE_CPLEX=0`, as appropriate.
7. There are two compilation modes: `debug` and `release`. These can be compiled with `make [debug or release]`, which creates the executable `vpc` in a new subdirectory `Debug` or `Release` of the main folder.
8. You can test the code with `make test_debug` or `make test_release`, depending on which version you compiled. Alternatively, run [`test/run_test.sh`](test/run_test.sh), which assumes the `Debug` version has been compiled.

## Dependencies

1. *Operating system*: This code has been tested on Linux and Mac operating systems. It is assumed that the installer has some familiarity and comfort with running commands in the terminal. Root access may be needed to get the proper dependencies installed. If you do not have root access, contact your system adminstrator for support.

2. *Compiler*: You will need a compiler (typically `g++` on Linux and `clang` on Mac). Some of the scripts use `bash`, but can probably be adapted to other shells. Your `bash` version should be at least version 4. Recently, MacOS comes with `zsh` instead of `bash`, and most of the time this should be fine, but if needed, use `homebrew` to install `bash` (via `homebrew install bash`).

3. *Branch-and-bound code*: The VPC code relies on [Cbc](https://github.com/coin-or/Cbc), which is installed via the [`setup/install_coin.sh`](setup/install_coin.sh) script. The script requires you to have `wget`.

  The project was previously extensively tested with Cbc version 2.9 (up to revision 2376 in the subversion history), though this backwards compatibility is no longer actively maintained. Performance in the main/trunk and 2.10 versions has been more unstable. Specifically, Cbc version 2.10 will probably *not* work. See the related issue [#12](https://github.com/akazachk/vpc/issues/12). In particular, as of 2020/03/27, we need to compile the development branch of Cbc with `SAVE_NODE_INFO` defined (using `ADD_CXXFLAGS="-DSAVE_NODE_INFO"`) to enable access to the `parentNode` code in `CbcModel`. In addition, if the Cbc commit is before [0f6ffed](https://github.com/coin-or/Cbc/commit/0f6ffed4c26daaf75edac2f87b70f3cc40cb12fd), we need to comment out the line `currentNode_ = NULL` in `CbcModel.cpp` around [line CbcModel.cpp:15392](https://github.com/coin-or/Cbc/blob/53f34cfea21360091608b02a041a962b2be7d6bc/src/CbcModel.cpp#L15390-L15391).

4. *Cbc dependencies*: For Cbc, you may need `gfortran`, `pkg-conf`, `LAPACK`, and `BLAS`. It may also be necessary to use `--with-cplex=false` as an option in the `coinbrew` commands, if [Osi](https://github.com/coin-or/Osi)'s configure script detects the CPLEX library through `CPXgetstat` but the CPLEX include directory is not found (see https://github.com/coin-or/coinbrew/issues/49).

5. *Environment variables*: There shoud be an environment variable `VPC_DIR` pointing to the local repository location, or this variable can be defined in each of the scripts: [`setup/install_coin.sh`](setup/install_coin.sh), [`test/run_test.sh`](test/run_test.sh), and others.

6. You may need to check compatibility with your version of `clang` or `g++` (for example, for the [inline variable features](https://github.com/akazachk/vpc/commit/0974799ec01e4a9135a58d48c03c5afa09756419#diff-b16b01e4c6bc7937a16cdc866415fcb3cc1c882a79dda4543dc626056a698f00) that this code uses from from C++17, at least version 7 of `g++` is needed).
If using Gurobi, make sure that the same compiler is used when running [`install_coin.sh`](setup/install_coin.sh) and the one used to generate `libgurobi_c++.a`, which can be rebuilt in the Gurobi directory (located at, say, `${GUROBI_HOME}`), under `${GUROBI_HOME}/src/build`.
The `Makefile` also assumes that `git` is at least version 2, to use the `-C` option to get the version of `Cbc` and `Clp`.

6. *Optional*: There are some optional libraries, such as `libbz2-dev`, that are linked to in [`Makefile`](Makefile). If these are missing, _remove the corresponding linking in the Makefile_.

## Execution

In more detail, to run the code, a basic call will look like:

```
${VPC_DIR}/[Debug or Release]/vpc -f filename -d num_disj_terms
```

The flag `-f` is used to specify the filename. Then `-d` takes an integer (the number of disjunctive terms to generate, though we also use this for the number of split disjunctions). Other options are described when calling the code with `--help` or `-h`.

## Advanced Execution

To run experiments with a set of instances, you can either use (1) [`scripts/run_experiments.sh`](scripts/run_experiments.sh), which will by default put output in the `results` folder, or (2) set up a batch of runs using [`scripts/prepare_batch.sh`](scripts/prepare_batch.sh), which will save a list of commands to a file `job_list_$MODE.txt`, where `$MODE` is the last argument into the `prepare_batch.sh` script.
You will either need to export the shell variable `VPC_DIR`, pointing to the root directory of the repository, or enter it when prompted by the script.

The `run_experiments.sh` script takes four arguments:
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
The [`merge.sh`](scripts/merge.sh) script also takes a third optional argument for where to save the output.

## Details

There are many parameters that can be set from the command line. Run with `-h` or `--help` option to see these parameters. There are several more that are not currently able to be set from the command line, as they are assumed to generally be "constants"; a description of these can be found in [`VPCParameters.hpp`](include/VPCParameters.hpp).

The key parts of the code are the classes [`Disjunction`](include/common/Disjunction.hpp) (an abstract class), [`CglVPC`](include/CglVPC.hpp), and [`PRLP`](include/PRLP.hpp). The user must specify a disjunction to `CglVPC` and a set of parameters (given as the struct `VPCParameters`). The disjunction should inherit from `Disjunction` and provide a set of disjunctive terms (the vector `terms`) and their number (`num_terms`); see for example the [`SplitDisjunction`](include/SplitDisjunction.hpp) class. Each term is a `DisjunctiveTerm` defined in [`Disjunction.hpp`](include/common/Disjunction.hpp) and contains the following members: 
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
4. Strengthening cuts -- see the related [MIP 2021 talk](https://www.youtube.com/watch?v=axqaOED4CXQ).

## Contact Information
[Aleksandr M. Kazachkov](https://akazachk.github.io),
[akazachk AT discreteopt DOT com](https://akazachk.github.io/email)
