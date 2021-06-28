# Compile the test mechanism under FOAM_USER_LIBBIN for consistency.
testMechLibdDir=$FOAM_USER_LIBBIN/unittests/pyjacTestMech/build

rm -rf $testMechLibdDir

cmake -B$testMechLibdDir -H. -DCMAKE_C_COMPILER=cc 

pushd $testMechLibdDir
    make
popd

# generate a symbolic link in FOAM_USER_LIBBIN to avoid modification of LD_LIBRARY_PATH
ln -sf $testMechLibdDir/libc_pyjac_test.so $FOAM_USER_LIBBIN/libc_pyjac_test.so
