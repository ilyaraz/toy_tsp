#!/bin/bash

for ((i = 10; i <= 200; i += 10))
do
    echo "n = $i"
    touch ../logs/log$i.txt
    for ((j = 0; j < 50; ++j))
    do
        python ../generators/random_euclidean.py $i | ../hk 2>>../logs/log$i.txt
        echo -n "."
    done
    echo " done"
done
