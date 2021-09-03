#!/bin/bash

M="1000"
H="0.01"
E="1e-3"

g++ -Ofast -ffast-math -march=native *.cpp

./a.out 2 0.1 1 $M $H $E
python3 create_pic_V.py
python3 create_pic_G.py 
mv plot_V.png plot_0_V_0.png
mv plot_G.png plot_0_G_0.png

