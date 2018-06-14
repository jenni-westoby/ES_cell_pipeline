#!/bin/bash

#rename command line arguments
path_to_ref_fasta=$1
path_to_ref_gtf=$2

./setup.sh setup
./RSEM_ref.sh make_ref $path_to_ref_gtf $path_to_ref_fasta
./make_indexes.sh $path_to_ref_gtf $path_to_ref_fasta
gunzip ES_cell_data/*


for i in ES_cell_data/*_1.fastq;
do
  num_jobs=`bjobs | wc -l`
  max_jobs=30

  #This prevents the number of queued jobs greatly exceeding 30.
  while [[ $num_jobs -gt $max_jobs ]];
  do
    sleep 100
    num_jobs=`bjobs | wc -l`
  done

  bsub -n8 -R"span[hosts=1]" -c 99999 -G team_hemberg -q normal -o $TEAM/temp.logs/output.$i -e $TEAM/temp.logs/error.$i -R"select[mem>100000] rusage[mem=100000]" -M100000 ./cell_level_analysis.sh $i
done

#make clean results matrices
./make_matrix.sh make_matrix RSEM
./make_matrix.sh make_matrix eXpress
./make_matrix.sh make_matrix Kallisto
./make_matrix.sh make_matrix Sailfish
./make_matrix.sh make_matrix Salmon
./make_matrix.sh make_matrix ground_truth
python generate.py Kallisto_real `pwd` Simulation/Kallisto_real_results
./clean_data.sh

#move data - some of this should move into setup.sh
cp Simulation/results_matrices/clean* raw_results/data/

#format data to make figures
cd raw_results
Rscript Figure2.R
Rscript Figure4.R
Rscript Figure5a.R
Rscript Figure5b.R
Rscript Figure6.R
Rscript SupplementaryFigure10.R
Rscript SupplementaryFigure11.R
Rscript SupplementaryFigure12.R

#make figure pdfs
cd ../figures/scripts
Rscript Figure2.R
Rscript Figure5.R
Rscript Figure6.R
Rscript SupplementaryFigure10.R
Rscript SupplementaryFigure11.R
Rscript SupplementaryFigure12.R
Rscript SupplementaryFigure16.R
Rscript SupplementaryFigure17.R
