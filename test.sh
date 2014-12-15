CURHOME=$PWD
R CMD INSTALL --configure-args="--with-hdf5=$CURHOME/library/hdf5" library/h5r
if [ $?!=0 ];then
        echo '[ERROR] h5r installation failed.'
	exit 1
fi



