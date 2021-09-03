#!/bin/bash

g++ -Ofast -march=native -ffast-math *.cpp
./a.out 0.1 1 1000 1e-3 1 1 123
python3 create_pic.py

