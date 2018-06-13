#!/bin/bash

if [ ! "$(ls -A Simulation/indices/Salmon_SMEM)" ]; then
  ./Simulation/Salmon-0.8.2_linux_x86_64/bin/salmon index -t Simulation/ref/reference.transcripts.fa -i Simulation/indices/Salmon_SMEM/transcripts_index_SMEM --type fmd -p 2
fi

if [ ! "$(ls -A Simulation/indices/Salmon_quasi)" ]; then
  ./Simulation/Salmon-0.8.2_linux_x86_64/bin/salmon index -t Simulation/ref/reference.transcripts.fa -i Simulation/indices/Salmon_quasi/transcripts_index_quasi --type quasi -k 31 -p 2
fi

if [ ! "$(ls -A Simulation/indices/Kallisto)" ]; then
  ./Simulation/kallisto_linux-v0.43.1/kallisto index -i Simulation/indices/Kallisto/transcripts.idx Simulation/ref/reference.transcripts.fa
fi

if [ ! "$(ls -A Simulation/indices/Sailfish)" ]; then
  export LD_LIBRARY_PATH=`pwd`/SailfishBeta-0.10.0_CentOS5/lib:$LD_LIBRARY_PATH
  export PATH=`pwd`/SailfishBeta-0.10.0_CentOS5/bin:$PATH
  LC_ALL=C ./Simulation/SailfishBeta-0.10.0_CentOS5/bin/sailfish index -p 2 -t Simulation/ref/reference.transcripts.fa -o Simulation/indices/Sailfish/ -k 31
fi

if [ ! "$(ls -A Simulation/indices/STAR)" ]; then
  Simulation/STAR/bin/Linux_x86_64/STAR --runThreadN 2 --runMode genomeGenerate --genomeDir Simulation/indices/STAR --genomeFastaFiles $2  --sjdbGTFfile $1
fi
