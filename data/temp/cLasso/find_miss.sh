#!/bin/bash

for i in $(seq 1 1 8000)
do
    if ! [ -f "./$i.RData" ]; then
        echo $i
    fi
done

