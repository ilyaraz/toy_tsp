#!/bin/bash

for ((i = 10; i <= 200; i += 10))
do
    echo "n = $i"
    touch log$i.txt
    for ((j = 0; j < 50; ++j))
    do
        python ../generators/random_euclidean.py $i | ../tsp 2>>log$i.txt
    done
done
