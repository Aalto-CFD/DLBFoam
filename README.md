# pyJacChemistryModel

An OpenFOAM chemistry model utilizing analytical Jacobian formulation, using [pyJac](https://github.com/SLACKHA/pyJac).

## pyJac Library Generation

install [pyJac](https://github.com/SLACKHA/pyJac).

We will generate the jacobian library for the gri mechanism, using an OpenFOAM tutorial.

### In OpenFOAM-6 environment, do the following:

```
cd $WM_PROJECT_USER_DIR
mkdir pyJac_gri
cp $FOAM_TUTORIALS/combustion/chemFoam/gri/chemkin/chem.inp pyJac_gri
cp $FOAM_TUTORIALS/combustion/chemFoam/gri/chemkin/therm.dat pyJac_gri
cd pyJac_gri
python -m pyjac --lang c --last_species N2 --input chem.inp --thermo therm.dat
```

### Generate the library

```
python -m pyjac.libgen --lang c -out $FOAM_USER_LIBBIN
```

### Compile the pyJac chemistry model inside src/thermophysicalModels/chemistryModel

```
wmake libso
```

## Usage

Check the tutorials given in tutorials folder.