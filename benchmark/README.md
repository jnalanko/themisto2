# Benchmarks (work in progess)

Three datasets:

* 1k random genomes from GTDB representative microbial genomes
* 3682 E. coli genomes
* Very small test dataset with 4 E. coli genomes

TODO: put the data somewhere, Zenodo?

There is one subdirectory per testcase. The sequences go into the subsubdirectory seqs.

There is a script `setup.bash` which will build unitigs with ggcat, and create the
required file-of-files. You should run the script so that the working directory
is at the dataset subdirectory.

There is a script `run_all.bash` that contains the sbwt and themisto2 commands.
This also needs to be run so that the working directory is in the dataset
subdirectory.

