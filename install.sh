cd src
make clean
make
cd ..

R CMD INSTALL library/Rcpp
R CMD INSTALL library/locfdr
R CMD INSTALL library/seqPatch
R CMD INSTALL library/pbh5
R CMD INSTALL library/h5r
Rscript install_Rpackages.R

