#!/bin/bash

CPUS=1
RAM=8
TEMP="."

while [ "$1" != "" ]; do
    case $1 in
    -h | --help )               usage; exit 1
                                ;;
    --cpus )                    shift; CPUS=$1
                                ;;
    --ram )                     shift; RAM=$1
                                ;;
    -T )                        shift; TEMP=$1
                                ;;
    * )                         BEDFILE=$1
                                ;;
    esac
    shift
done

sort -T $TEMP -S "$RAM"G --parallel=$CPUS -k1,1V -k2,2n -k3,3n -k6,6 -k11,11 -k12,12 -k4,4 -k14,14 $BEDFILE
