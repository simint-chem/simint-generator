#!/bin/bash

mkdir -p dat

mkdir -p dat/100.10    dat/100.10.rand
./gen_num.py -r -np 10 -ns 100 -o dat/100.10.rand/1.dat
./gen_num.py -r -np 10 -ns 100 -o dat/100.10.rand/2.dat
./gen_num.py -r -np 10 -ns 100 -o dat/100.10.rand/3.dat
./gen_num.py -r -np 10 -ns 100 -o dat/100.10.rand/4.dat
./gen_num.py -np 10 -ns 100 -o dat/100.10/1.dat
./gen_num.py -np 10 -ns 100 -o dat/100.10/2.dat
./gen_num.py -np 10 -ns 100 -o dat/100.10/3.dat
./gen_num.py -np 10 -ns 100 -o dat/100.10/4.dat

mkdir -p dat/2.2    dat/2.2.rand
./gen_num.py -r -np 2 -ns 2 -o dat/2.2.rand/1.dat
./gen_num.py -r -np 2 -ns 2 -o dat/2.2.rand/2.dat
./gen_num.py -r -np 2 -ns 2 -o dat/2.2.rand/3.dat
./gen_num.py -r -np 2 -ns 2 -o dat/2.2.rand/4.dat
./gen_num.py -np 2 -ns 2 -o dat/2.2/1.dat
./gen_num.py -np 2 -ns 2 -o dat/2.2/2.dat
./gen_num.py -np 2 -ns 2 -o dat/2.2/3.dat
./gen_num.py -np 2 -ns 2 -o dat/2.2/4.dat


mkdir -p dat/1.1.rand
./gen_num.py -r -np 1 -ns 1 -o dat/1.1.rand/1.dat
./gen_num.py -r -np 1 -ns 1 -o dat/1.1.rand/2.dat
./gen_num.py -r -np 1 -ns 1 -o dat/1.1.rand/3.dat
./gen_num.py -r -np 1 -ns 1 -o dat/1.1.rand/4.dat
