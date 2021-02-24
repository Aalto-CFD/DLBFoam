RED='\e[1;31m'
NC='\e[0m' # No Color
GREEN='\e[1;32m'
YELLOW='\e[1;33m'
DARKGRAY='\e[1;30m'

echo
echo -e "${YELLOW}Validate pyJac chemistry module against reference data on perfectly stirred reactor:${DARKGRAY}"
echo

./PSRTest -case testCase

if [ $? -eq 0 ]; then
    echo
    echo -e "${GREEN}Validation checks PASSED.${NC}"
else
    echo
    echo -e "${RED}Validation FAILED.${NC}"
    echo -e "${RED}See output above for further information.${NC}"
fi



echo
echo -e "${YELLOW}Test error control on mechanism consistency in pyjacLoadBalancedChemistryModel:${DARKGRAY}"
cp testCase/constant/thermophysicalProperties testCase/constant/thermophysicalProperties.orig
cp testCase/constant/thermophysicalProperties.illdefined testCase/constant/thermophysicalProperties
./PSRTest -case testCase > /dev/null 2>&1
if [ $? -eq 1 ]; then
    echo -e "${GREEN}PASSED.${NC}"
else
    echo -e "${RED}FAILED. pyjacLoadBalancedChemistryModel.C did not catch the species mismatch error.${NC}"
fi
cp testCase/constant/thermophysicalProperties.orig testCase/constant/thermophysicalProperties
rm testCase/constant/thermophysicalProperties.orig
