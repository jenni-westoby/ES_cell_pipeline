#!/bin/bash

./setup.sh setup
./RSEM_ref.sh make_ref path_to_ref_gtf path_to_ref_fasta
./make_indexes.sh

for i in ES_cell_data/*_1.fastq;
do
  ./cell_level_analysis.sh $i
done

#make clean results matrices
./make_matrix.sh make_matrix RSEM
./make_matrix.sh make_matrix eXpress
./make_matrix.sh make_matrix Kallisto
./make_matrix.sh make_matrix Sailfish
./make_matrix.sh make_matrix Salmon_align
./make_matrix.sh make_matrix Salmon_SMEM
./make_matrix.sh make_matrix Salmon_quasi
./make_matrix.sh make_matrix ground_truth
./make_matrix.sh make_matrix Kallisto_real
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
