#!/bin/sh
NC='\e[0m' # No Color
GREEN='\e[1;32m'
YELLOW='\e[1;33m'
DARKGRAY='\e[1;30m'

echo -e "${YELLOW}Validate ISAT-DLBFoam coupling on a chemFoam test:${DARKGRAY}"

# Source tutorial run functions
. $WM_PROJECT_DIR/bin/tools/RunFunctions

# Redirect all output to /dev/null
exec > /dev/null 2>&1

rm -rf gri
cp -r $WM_PROJECT_DIR/test/chemistry/gri .

cd gri

    application=$(getApplication)

    runApplication chemkinToFoam \
                chemkin/chem.inp chemkin/therm.dat chemkin/transportProperties \
                constant/reactions constant/speciesThermo

    # Set ISAT
    sed -i '/chemistry       on;/a#includeEtc     "caseDicts/solvers/chemistry/TDAC/chemistryPropertiesFlame.cfg"\n#remove         reduction' constant/chemistryProperties
    foamDictionary -expand -entry tabulation/tolerance -set 1e-5 constant/chemistryProperties

    # Test standard implementation
    runApplication -s standard chemFoam
    mv chemFoam.out chemFoam_standard.out
    (cd validation && ./Allrun $*)
    cp validation/OF_vs_CHEMKINII.eps OF_vs_CHEMKINII_standard.eps

    # Test DLB
    foamDictionary -entry libs -set '("libchemistryModel_DLB.so" )' system/controlDict
    foamDictionary -entry loadbalancing -set {} constant/chemistryProperties
    foamDictionary -entry loadbalancing/active -add false constant/chemistryProperties
    foamDictionary -entry refmapping -set {} constant/chemistryProperties
    foamDictionary -entry refmapping/active -add false constant/chemistryProperties

    foamDictionary -entry chemistryType/method -set loadBalanced constant/chemistryProperties
    runApplication -s DLB chemFoam
    mv chemFoam.out chemFoam_DLB.out
    (cd validation && ./Allrun $*)
    cp validation/OF_vs_CHEMKINII.eps OF_vs_CHEMKINII_DLB_ISAT.eps

    # Re-enable output for the unit test
    exec > /dev/tty 2>&1

    # Unit test to compare outputs
    if diff chemFoam_standard.out chemFoam_DLB.out > /dev/null ; then
         echo -e "${GREEN}PASSED.${NC}"
    else
        echo -e "${RED}FAILED. Check the ISAT-DLBFoam implementation.${NC}"
    fi

cd ../

rm -rf gri
