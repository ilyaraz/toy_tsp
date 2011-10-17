#!/bin/bash

for ((i = 50; i <= 500; i += 50))
do
    echo "n = $i"
    touch ../logs/log$i.txt
    for ((j = 0; j < 10; ++j))
    do
        python ../generators/random_euclidean.py $i | ../hk 2>>../logs/log$i.txt
        echo -n "."
    done
    echo " done"
done
