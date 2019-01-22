# V-Polyhedral Disjunctive Cuts
#### Aleksandr M. Kazachkov
#### Based on joint work with Egon Balas

This project contains the code for one implementation of the V-polyhedral disjunctive cut paradigm

## Dependencies

The code is tested with `Cbc 2.9`, which can be obtained via the script `install_coin.sh` in the `scripts` directory.

## Compilation

There are two compilation modes: `debug` and `release`. These can be compiled with `make [debug or release]`, which create the executable `vpc` in a new subdirectory `Debug` or `Release` of the main folder.

## Execution

You can test the code by running `scripts/run_test.sh`. This assumes the `Debug` version has been compiled.

## Contact Information
Aleksandr M. Kazachkov
aleksandr DOT kazachkov AT polymtl DOT ca
