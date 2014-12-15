# compile source code
cd src
make clean
make
cd ..

# compile HDF5
CURHOME=$PWD

cd HDF5/hdf5
./configure --prefix=$CURHOME/library/hdf5 --enable-cxx
make 
make check
make install
make check-install

# setup environment variables
cd $CURHOME

echo "export QDNAMODHOME=$CURHOME" > ./setenv.sh
echo "export LD_LIBRARY_PATH=$CURHOME/library/hdf5/lib:\$LD_LIBRARY_PATH" >> ./setenv.sh
echo "export PATH=$CURHOME/bin:\$PATH" >> ./setenv.sh

chmod u+x ./setenv.sh
source ./setenv.sh

# install R packages

R CMD INSTALL --configure-args="--with-hdf5=$CURHOME/library/hdf5" library/h5r
R CMD INSTALL library/Rcpp
R CMD INSTALL library/locfdr
R CMD INSTALL library/seqPatch
Rscript install_Rpackages.R




