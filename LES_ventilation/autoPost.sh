#!/bin/bash

# A
cd ./A1_5.45/
mpirun -np 12 ~/CharLES/nextgen_scalar/bin/charles.exe -i 00_post_avg.in

cd ../A1_5.90/
mpirun -np 12 ~/CharLES/nextgen_scalar/bin/charles.exe -i 00_post_avg.in

# C
cd ../C1_5/
mpirun -np 12 ~/CharLES/nextgen_scalar/bin/charles.exe -i 00_post_avg.in

cd ../C1_5.45/
mpirun -np 12 ~/CharLES/nextgen_scalar/bin/charles.exe -i 00_post_avg.in

cd ../C1_5.90/
mpirun -np 12 ~/CharLES/nextgen_scalar/bin/charles.exe -i 00_post_avg.in

# D
cd ../D1_5/
mpirun -np 12 ~/CharLES/nextgen_scalar/bin/charles.exe -i 00_post_avg.in
cd ../D1_5.45/
mpirun -np 12 ~/CharLES/nextgen_scalar/bin/charles.exe -i 00_post_avg.in
cd ../D1_5.90/
mpirun -np 12 ~/CharLES/nextgen_scalar/bin/charles.exe -i 00_post_avg.in


# E
cd ../E1_5/
mpirun -np 12 ~/CharLES/nextgen_scalar/bin/charles.exe -i 00_post_avg.in
cd ../E1_5.45/
mpirun -np 12 ~/CharLES/nextgen_scalar/bin/charles.exe -i 00_post_avg.in
cd ../E1_5.90/
mpirun -np 12 ~/CharLES/nextgen_scalar/bin/charles.exe -i 00_post_avg.in
cd ../E1_10/
mpirun -np 12 ~/CharLES/nextgen_scalar/bin/charles.exe -i 00_post_avg.in
cd ../E1_20/
mpirun -np 12 ~/CharLES/nextgen_scalar/bin/charles.exe -i 00_post_avg.in

# Z
cd ../Z1/
mpirun -np 12 ~/CharLES/nextgen_scalar/bin/charles.exe -i 00_post_avg.in
cd ../Z1.45/
mpirun -np 12 ~/CharLES/nextgen_scalar/bin/charles.exe -i 00_post_avg.in






