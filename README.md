qDNAmod
======
qDNAmod is a toolkit for quantitative detection of DNA modification heterogeneity from SMRT sequencing data. 

System
======
Linux

Dependencies
============
1. R from http://www.r-project.org/
(Add directory containing binary file 'Rscript' to your enviroment variable PATH)

2. HDF5 from http://www.hdfgroup.org/HDF5/release/obtainsrc.html. 
(Specify enviroment varible HDF5_HOME as the directory of installed HDF5)

Installation
============
Unzip the downloaded file, enter the directory and type ./install.sh 

Add (the unzipped directory of qDNAmod)/bin to your PATH

Set evironment variable QDNAMODHOME=(the unzipped directory of qDNAmod)

(g++ is required for installation)

Make sure your computer can access internet during the installation.

License
=======
GPL2


Input
=====
The Input of qDNAmod is .cmp.h5 files (aligned reads).


Example
=======
http://bioinfo.au.tsinghua.edu.cn/member/zfeng/example_release.tar.gz

Maintainer
==========
Zhixing Feng

Contact
=======
zxfeng.thu@gmail.com


Usage
=====
There are 3 steps: pileup reads, learn prior and detect modification proportion.

Step 1: pileup reads. Firstly, we need to pileup reads position by position.

	qDNAmod_pileup  [-r <reagent>] [-m <mapQVthreshold>] <cmpH5file> <outdir>

	where:
   		-r <reagent>,  --reagent <reagent>
     		chemistry used for SMRT sequencing, candidates are: "C2", default is "C2"

   		-m <mapQVthreshold>,  --mapQVthreshold <mapQVthreshold>
     		minimal mapQV, default is 255

   		<cmpH5file>
     		(required)  aligned data in cmpH5 format

   		<outdir>
     		(required)  output directory

Step 2: learn prior. To reduce uncertainty of estimation, we need to roughly estimate IPD distribution of modified bases.

	qDNAmod_prior  <native pileup dir> <WGA pileup dir> <outdir>

	where:

   		<native pileup dir>
     		(required)  native pileup data directory

   		<WGA pileup dir>
     		(required)  WGA pileup data directory

   		<outdir>
     		(required)  output directory

Step 3: detect modification proportion. Now we fit a Bayesian mixture model using prior estimated in step 2.

	qDNAmod_detect <native pileup dir> <WGA pileup dir> <prior dir> <outdir>

	where:
   		<native pileup dir>
     		(required)  native pileup data directory

   		<WGA pileup dir>
     		(required)  WGA pileup data directory

   		<prior dir>
     		(required)  directory containing estimated prior

   		<outdir>
     		(required)  output directory

Output
======
In the <outdir> of qDNAmod_detect, there are .txt files with the names "detect_<chromosome_name>.txt" corresponding to results for <chromosome_name>. 
In each .txt file, the columns are:

column 1 (locus): genome locus.

column 2 (strand): strand, 0 means forward, 1 means backward.

column 3 (prop): estimated modification proportion.

column 4 (N_1): expectation of number of kinetic variant bases.  

column 5 (N_0): expectation of number of normal bases. 

column 6 (avg_n): average number of times each base being sequenced.

column 7 (cvg_wga): coverage of WGA sample. 






