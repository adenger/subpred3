#!/bin/bash

# Make sure that conda env is activated, and screen session is started.

# for nb in *.ipynb ; do jupyter nbconvert --to notebook --inplace --execute $nb ; done
for nb in $(find -name '*.ipynb' | sort); do jupyter nbconvert --to notebook --inplace --execute $nb ; done