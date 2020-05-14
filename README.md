# V-Polyhedral Disjunctive Cuts
#### Aleksandr M. Kazachkov
#### Based on joint work with Egon Balas

This project contains the code for one implementation of the V-polyhedral disjunctive cut paradigm.

## Dependencies

The code relies on [Cbc](https://github.com/coin-or/Cbc). It was extensively tested with version 2.9 (up to revision 2376 in the subversion history), while performance in the trunk and 2.10 versions has been more unstable. I have not yet tracked down the precise reason; see the related open GitHub issue. In particular, as of 2020/03/27, I need to use custom `CbcModel` files to make the code work with the latest `Cbc` versions. I provide these files in the `lib` directory.

The user needs to export the shell variable `PROJ_DIR`, pointing to the repository location, or define it in `makefile`, `scripts/install_coin.sh`, `test/run_test.sh`, and possibly other locations.

You may need to use a compatible version of `g++`, and make sure that the same compiler is used when running `install_coin.sh` and the one used to generate `libgurobi_c++.a`, which can be rebuilt in the `Gurobi` directory under `src/build`.

## Compilation

There are two compilation modes: `debug` and `release`. These can be compiled with `make [debug or release]`, which create the executable `vpc` in a new subdirectory `Debug` or `Release` of the main folder.

## Execution

You can test the code by running `scripts/run_test.sh`. This assumes the `Debug` version has been compiled.

In more detail, you run the code with:

```
[Debug or Release]/vpc -f filename -d num_disj_terms/num_disjunctions
```

The flag `-f` is used to specify the filename. Then `-d` takes an integer (the number of disjunctive terms to generate, though we also use this for the number of split disjunctions). Other options are described when calling the code with `--help` or `-h`.

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

## Contact Information
Aleksandr M. Kazachkov
aleksandr DOT kazachkov AT polymtl DOT ca
