# README to describe the test data:

## Files:
___
* 0_percent_affected_phecode_matrix.txt
* 5_percent_affected_phecode_matrix.txt
* 20_percent_affected_phecode_matrix.txt
* 50_percent_affected_phecode_matrix.txt
* 80_percent_affected_phecode_matrix.txt
* cftr_and_vw_test_phecode_matrix.txt
* gene_info.txt
* shared_segment_test_data.txt

## Description:
___
* **_percent_affected_phecode_matrix.txt**
    * These files have three columns: grids, 286.11, and 499. The first column has the random grid ids. The second column has the affect status for vonwillibrands while the third column has the status for cftr. The number of affected individuals depends on what percentage the file is Ex: 5_percent means 5/100 people are affected with phenotype. There are 100 random individuals who were chosen from the larger set of MEGA individuals

* **cftr_and_vw_test_phecode_matrix.txt**
    * This file is the same as the _percent_affected_phecode_matrix.txt files but has the original affected statuses from the biobanks

* **gene_info.txt**
    * has information about the gene in four columns. First column is name. SEcond column is the chromosome. Third column is the gene start position. And the fourth column is the gene end position. These positions are taken from GnomAD build37

* **shared_segment_test_data.txt**
    * file that has the shared ibd pairwise segment data. The columns are: pair id 1, individual 1 phase, pair id 2, individual 2 phase, chromosome, segment start, segment end, and segment length. The pairs are matched to those 100 indivduals in the _percent_affected_phecode_matrix.txt.
