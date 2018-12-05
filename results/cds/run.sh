#! /bin/bash

../../sources/comparing_distance_scale/main.py ../../data/13.txt ../../data/14.txt
gnuplot ../../sources/scripts/cds.gnu
gnuplot ../../sources/scripts/cds_without_ej.gnu


