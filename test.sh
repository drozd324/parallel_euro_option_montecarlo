#!/bin/bash

echo "main"
time ./main

echo " "

echo "parallel_main_omp"
time ./parallel_main_omp

echo " "

echo "parallel_main_mpi"
time mpirun -n 4 ./parallel_main_mpi
