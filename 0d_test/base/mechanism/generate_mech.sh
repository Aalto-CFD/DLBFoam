#!/bin/bash
./Allclean
rm -rf pyJac
mkdir pyJac
cp -r cmake/* pyJac

cd pyJac
python -m pyjac --lang c --last_species N2 --input ../chemkin/chem.cti
./runCmake.sh

cd ../
cd create_OF_mechanism/
python3 fit_thermoProps.py
./makeMechFileForOF.sh

cp -r ../pyJac/build/libc_pyjac.so foam_mech_files
mv foam_mech_files ../../constant

