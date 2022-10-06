#!/bin/bash

# Make sure that conda env is activated, and screen session is started.


# human=(human_*evaluation*)
# yeast=(yeast_*evaluation*)
# files=("${human[@]}" "${human[@]}")
files=(yeast_*evaluation*)

for file in "${files[@]}"; do
    # echo "$file"
    jupyter nbconvert --to notebook --inplace --execute "$file"
done
