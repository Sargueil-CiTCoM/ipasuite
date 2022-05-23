#!/bin/bash

#set -xe
snakemake $* -c8 --use-conda --keep-going

for log in logs/*.log ; do
    if [ -s "$log" ] ; then
        echo $log
        echo -ne '\t'
        cat $log
    fi
done
