#!/bin/bash

# Make sure that conda env is activated, and screen session is started.

files=(meta_*evaluation*)

for file in "${files[@]}"; do
    # echo "$file"
    jupyter nbconvert --to notebook --inplace --execute "$file"
done
