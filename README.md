# pyJacChemistryModel
![OpenFOAM 8](https://img.shields.io/badge/OpenFOAM-8-brightgreen)
![OpenFOAM 7](https://img.shields.io/badge/OpenFOAM-7-brightgreen)

An OpenFOAM chemistry model utilizing analytical Jacobian formulation, using [pyJac](https://github.com/SLACKHA/pyJac).

## pyJac Library Generation

install [pyJac](https://github.com/SLACKHA/pyJac).

We will generate the jacobian library for the GRI mechanism, using an OpenFOAM tutorial.

### In OpenFOAM-6 environment, do the following:

```
cd $WM_PROJECT_USER_DIR
mkdir pyJac
cp $FOAM_TUTORIALS/combustion/chemFoam/gri/chemkin/chem.inp pyJac
cp $FOAM_TUTORIALS/combustion/chemFoam/gri/chemkin/therm.dat pyJac
cd pyJac
python -m pyjac --lang c --last_species N2 --input chem.inp --thermo therm.dat
```

### Generate the library

```
python -m pyjac.libgen --source_dir ./out --lang c -out $FOAM_USER_LIBBIN
```

### Compile the pyJac chemistry model inside src/thermophysicalModels/chemistryModel

```
wmake libso
```

## Usage

Check the tutorials given in tutorials folder.
