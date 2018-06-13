#!/bin/bash

base=`echo $1 |awk -F/ '{print $2}'`
filename=`echo $base |awk -F_ '{print $1}'`

./quality_control.sh quality_control $filename ES_cell_data '_1.fq' 'raw'
./simulate.sh simulate $filename ES_cell_data
./quality_control.sh qualitycontrol $filename Simulation/data/simulated '_1.fq' "simulated"
./quantify.sh Kallisto $filename
./quantify.sh eXpress $filename
./quantify.sh Salmon $filename
./quantify.sh RSEM $filename
./quantify.sh Sailfish $filename
./quantify_real_data.sh Kallisto $filename
