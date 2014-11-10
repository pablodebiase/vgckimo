#!/bin/bash
program=$(basename $0 | sed 's/\(.*\)\..*/\1/' )
module='praxis math'
levmarlibs='-llapack -lblas -lm' #-lf2c'
if [ "$1" == "intel" ]; then
  #intel optimized
  intellibs='-limf -lifcore'
  echo Compiling $(echo $program | tr '[:lower:]' '[:upper:]') using INTEL compilers:
  objects=''
  for mod in $module; do
    echo Compiling object $mod ...
    ifort -Warn all -ip -O3 -unroll -xHost -c $mod'.f90'
    objects=$mod.o' '$objects
  done
#compile levmar
  cobj='lm Axb misc lmlec lmbc lmblec lmbleic'
  cobjects=''
  for co in $cobj; do
    icc  -O3 -ip -xHost -unroll -Wall -c -o ${co}.o ${co}.c
    cobjects=${co}.o' '$cobjects
  done
  echo Compiling object $program ...
  icpc -Wall -openmp -ip -O3 -xHost -c $program'.cpp' 
  echo Linking $program ...
  icpc -Wall -openmp $intellibs $objects $cobjects $program'.o' -o $program'-intel' $levmarlibs 
else
  # gnu optimized
  gnulibs='-lgfortran'
  echo Compiling $(echo $program | tr '[:lower:]' '[:upper:]') using GNU compilers:
  objects=''
  for mod in $module; do
    echo Compiling object $mod ...
    gfortran -Wall -O3 -march=native -ffree-line-length-none -c $mod'.f90'
    objects=$mod.o' '$objects
  done
#compile levmar
  cobj='lm Axb misc lmlec lmbc lmblec lmbleic'
  cobjects=''
  for co in $cobj; do
    gcc   -O3 -funroll-loops -Wall -march=native  -c -o ${co}.o ${co}.c
    cobjects=${co}.o' '$cobjects
  done
  echo Compiling object $program ...
  g++ -fopenmp -Wall -O3 -march=native -c $program'.cpp'
  echo Linking $program ...
  g++ -fopenmp -Wall $objects $cobjects $program'.o' -o $program'-gnu' $gnulibs $levmarlibs
fi
echo Removing temporary files ...
rm -f $program'.o' $objects *__* $cobjects
echo $(echo $program | tr '[:lower:]' '[:upper:]') compiled.
