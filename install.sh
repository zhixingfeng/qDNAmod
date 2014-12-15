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
echo "export QDNAMODHOME=$CURHOME" > ./setenv.sh
echo "export LD_LIBRARY_PATH=$CURHOME/library/hdf5/lib:\$LD_LIBRARY_PATH" >> ./setenv.sh
chmod u+x ./setenv.sh
source ./setenv.sh

# install R packages
cd $CURHOME

R CMD INSTALL --configure-args="--with-hdf5=$CURHOME/library/hdf5" library/h5r
if [ $?!=0 ];then
	echo '[ERROR] h5r installation failed.'
	exit 1
fi
R CMD INSTALL library/Rcpp
if [ $?!=0 ];then
        echo '[ERROR] Rcpp installation failed.'
        exit 1
fi
R CMD INSTALL library/locfdr
if [ $?!=0 ];then
        echo '[ERROR] locfdr installation failed.'
        exit 1
fi
R CMD INSTALL library/seqPatch
if [ $?!=0 ];then
        echo '[ERROR] seqPatch installation failed.'
        exit 1
fi
R CMD INSTALL library/pbh5
if [ $?!=0 ];then
        echo '[ERROR] pbh5 installation failed.'
        exit 1
fi
Rscript install_Rpackages.R
if [ $?!=0 ];then
        echo '[ERROR] online R packages installation failed. Check your network.'
        exit 1
fi




