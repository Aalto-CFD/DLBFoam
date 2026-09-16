#!/bin/bash
set -e # Any subsequent(*) commands which fail will cause the shell script to exit immediately

# Set some default values:
mechanism=unset
pyjac_output=unset
compile=0
help=0

# Absolute path to this script
script_path=$(readlink -f "$0")
cmake_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# parent dir of output
input_rootdir=$(dirname "$(realpath -s $pyjac_output)")
lib_path=$input_rootdir/lib
foam_dir=$input_rootdir/foam


usage()
{
  echo "Usage: arguments  [ -m | --mechanism MECHANISM]
                          [ -i | --input path/to/pyjac/output ]
                          [ -c | --compile ]
                          [ -h | --help ]"
  exit 2
}


arg_check()
{
  if [ "$1" == "unset" ]; then
    echo
    echo "ERROR: Required arguments not given:"
    echo
    usage
    exit 2
  fi
}


check_cmake()
{
  cmake_com="cmake"
  if ! type $cmake_com &> /dev/null ; then
    echo "cmake command could not be found."
    exit 0
  fi
}


init_build_dir()
{
  # Copy original pyjac output for further modification and building
  mkdir -p  "$lib_path"
  rm -rf "$lib_path"/*
  cp -r "$pyjac_output" "$lib_path/src"

  # copy CMake directives to pyjac output
  cp "$cmake_dir"/CMakeLists.txt "$lib_path/"
  cp "$cmake_dir"/runCmake.sh "$lib_path/"
}


# DLBFoam requires function calls for number of species
modify_pyjac_src()
{
  if [ ! -f "$lib_path/src/chem_utils.h" ]; then
    echo "ERROR: $lib_path/src/chem_utils.h not found."
    exit 0
  fi
  sed -i '4 a int PYJAC_NSP();\nint PYJAC_FWD_RATES();' "$lib_path/src/chem_utils.h"
  sed -i '3 a int PYJAC_NSP()\n{\n    return NSP;\n}\nint PYJAC_FWD_RATES()\n{\n    return FWD_RATES;\n}' "$lib_path/src/chem_utils.c"
}


compile_pyjac_lib()
{
  pushd "$lib_path" > /dev/null
    ./runCmake.sh
  popd > /dev/null
  echo
  echo
}


#Grep the mechanism file of the pyjac in order to retain the right species order
rewrite_species()
{
  mechPath="$lib_path/src/mechanism.h"
  #Determine the species within a pyjac output file
  grep -A100000 'Species Indexes' "$mechPath"  > tmp.txt
  grep -B100000 '*/' tmp.txt > tmp2.txt
  sed '1d;$d' tmp2.txt > tmp3.txt
  awk '{ print $2}' tmp3.txt > tmp4.txt
  rm -f tmp.txt tmp2.* tmp3.*
  NS=$(cat tmp4.txt | wc -l)

  # rewrite the species
  echo "species" >> tmp
  echo "$NS" >> tmp
  echo "(" >> tmp

  rm -f "$foam_dir/species.foam"
  cat tmp tmp4.txt > "$foam_dir/species.foam"
  printf ");" >> "$foam_dir/species.foam"
  rm -f tmp*
}


write_reactions_note()
{
  sed -i '2 i\// - The reaction kinetics is read from pyJac source files.' "$foam_dir/reactions.foam"
}


howto_include()
{
  inc_f="$foam_dir/include_example.txt"
  rm -f "$inc_f"
  printf "%s\n" \
    "system/controlDict:" \
    "libs" \
    "(" \
      "    \"libchemistryModel_DLB.so\"" \
      "    \"libODE_DLB.so\"" \
      "    \"$lib_path/build/libc_pyjac.so\"" \
    ");" \
    ""\
    "constant/thermophysicalProperties:" \
    "    #include \"$foam_dir/species.foam\"" \
    "    #include \"$foam_dir/thermo.foam\"" \
    "" \
    "constant/chemistryProperties:" \
    "    #include \"$foam_dir/reactions.foam\"" \
    "" > "$inc_f"
  cat "$inc_f"
}



# ------------------ main ---------------------- #

PARSED_ARGUMENTS=$(getopt -a -n pyjac2foam -o m:i:ch --long mechanism:,input:,compile,help -- "$@")
VALID_ARGUMENTS=$?
if [ "$VALID_ARGUMENTS" != "0" ]; then
  usage
fi

#echo "PARSED_ARGUMENTS is $PARSED_ARGUMENTS"
eval set -- "$PARSED_ARGUMENTS"
while :
do
  case "$1" in
    -m | --mechanism)   mechanism="$2"      ; shift 2 ;;
    -i | --input)       pyjac_output="$2"   ; shift 2 ;;
    -c | --compile)     compile=1           ; shift   ;;
    -h | --help)        help=1              ; shift   ;;
    # -- means the end of the arguments; drop this, and break out of the while loop
    --) shift; break ;;
    # If invalid options were passed, then getopt should have reported an error,
    # which we checked as VALID_ARGUMENTS when getopt was called...
    *) echo "Unexpected option: $1 - this should not happen."
       usage ;;
  esac
done

if [ $help == 1 ]; then
    usage
    exit 2
fi

arg_check "$mechanism"
arg_check "$pyjac_output"

if [ ! -d "$pyjac_output" ]; then
   echo "ERROR: $pyjac_output path not found."
   exit 0
fi

if ! command -v ct2foam &>/dev/null; then
    echo "ERROR: ct2foam not found. Check your dependencies." >&2
    exit 1
fi


# -- Execute functionalities -- #

init_build_dir
modify_pyjac_src
if [ $compile == 1 ]; then
    check_cmake
    compile_pyjac_lib
fi

ct2foam -i "$mechanism" -o "$foam_dir" -T 1000.0

rewrite_species
write_reactions_note

printf '%s\nExample how to include files in OpenFOAM:\n'
howto_include

printf '%s\n\nExiting Succesfully.\n'
