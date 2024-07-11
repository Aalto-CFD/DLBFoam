# ISAT-DLBFoam Coupling Validation Script

This script validates the ISAT-DLBFoam coupling on a chemFoam test case. It performs the following steps:

1. Copies the necessary files from $WM_PROJECT_DIR/test/chemistry/gri to a local directory 'gri'.
2. Converts chemistry files using chemkinToFoam.
3. Configures ISAT settings in constant/chemistryProperties.
4. Runs chemFoam with standard implementation and stores results.
5. Runs chemFoam with DLBFoam implementation and stores results.
6. Compares the outputs from standard and DLB implementations using a unit test.

## Dependencies

- OpenFOAM environment with necessary libraries and tools.

## Usage

Ensure that your OpenFOAM environment is properly set up and configured. Modify paths and settings in the script as necessary to match your specific setup.
