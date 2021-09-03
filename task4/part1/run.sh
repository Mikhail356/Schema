#!/bin/bash

g++ -Ofast -ffast-math -march=native *.cpp

for g in 1 2 3 4
do
    for v in 1 2 3 4
    do
        ./a.out 0.1 1 1000 1e-4 $v $g 1e-3 > include_in_tex$v$g.tex
    done
done


