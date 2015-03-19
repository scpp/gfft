#!/bin/bash

declare -a tharr=("1" "2" "3" "4")
declare -a parr=("1" "2" "3" "4" "5" "6")


echo "============ Powers of 2 ================"

#for n in "${parr[@]}"
for n in `seq 1 12`
do
    cd ./build
    cmake -H.. -B. -DMYAUTO=1 -DMYNUM="(1<<${n})" .; time make
    src/gfft_inp
    read -p "Press [Enter] key to continue..."
done
