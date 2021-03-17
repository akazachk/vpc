# Instructions for setting up conda environment

To export to a different system use:
        conda env export > environment.yml
and to recreate the environment:
        conda env create -f environment.yml

You also need to add the following to your bashrc, and make sure to link against the `CONDA_LIB` directory when compiling your code.
    export CONDA_LIB=${HOME}/.conda/envs/vpc/lib
    export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${CONDA_LIB}
