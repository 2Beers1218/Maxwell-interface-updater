#!/bin/bash

# Define the list
#
real_width=(1.0)
real_intensity=(1.0)

# Trainings set with only real imput
for rw in "${real_width[@]}"
do
  for ri in "${real_intensity[@]}"
  do
    # Load the default options:
    cp default_options.prm options.prm

    # Set the parameter
    sed -i "s/set beam width real       = 100/set beam width real       = ${rw}/g" options.prm
    sed -i "s/set beam intensity real   = 1.0/set beam intensity real   = ${ri}/g" options.prm

    # Create the trainings data
    mpirun -np 2 ./build/KirasFM
  
    # Store the interface
    cat interface_0_1_1.csv >> ../training_data/interface_0_1_1_test.csv
    cat interface_0_1_4.csv >> ../training_data/interface_0_1_4_test.csv

    cat interface_1_0_1.csv >> ../training_data/interface_1_0_1_test.csv
    cat interface_1_0_4.csv >> ../training_data/interface_1_0_4_test.csv

    cat solution_0_1.csv >> ../training_data/solution_0_1_test.csv
    cat solution_1_1.csv >> ../training_data/solution_1_1_test.csv


  done
done

real_width=(1.0)
imag_width=(1.0)
real_intensity=(1.0)
imag_intensity=(1.0)
# Trainings set with only real imput
for rw in "${real_width[@]}"
do
  for rc in "${real_center[@]}"
  do
    for ri in "${real_intensity[@]}"
    do
      for iw in "${imag_width[@]}"
      do
        for ic in "${imag_center[@]}"
        do
          for ii in "${imag_intensity[@]}"
          do
            # Load the default options:
            cp default_options.prm options.prm

            # Set the parameter
            sed -i "s/set beam width real       = 100/set beam width real       = ${rw}/g" options.prm
            sed -i "s/set beam center real      = 0.5/set beam center real      = ${rc}/g" options.prm 
            sed -i "s/set beam intensity real   = 1.0/set beam intensity real   = ${ri}/g" options.prm
            sed -i "s/set beam width imag       = 100/set beam width imag       = ${iw}/g" options.prm
            sed -i "s/set beam center imag      = 0.5/set beam center imag      = ${ic}/g" options.prm 
            sed -i "s/set beam intensity imag   = 0.0/set beam intensity imag   = ${ii}/g" options.prm

            # Create the trainings data
            mpirun -np 2 ./build/KirasFM

            # Store the interface
            cat interface_0_1_1.csv >> ../training_data/interface_0_1_1_test.csv
            cat interface_0_1_4.csv >> ../training_data/interface_0_1_4_test.csv

            cat interface_1_0_1.csv >> ../training_data/interface_1_0_1_test.csv
            cat interface_1_0_4.csv >> ../training_data/interface_1_0_4_test.csv

            cat solution_0_1.csv >> ../training_data/solution_0_1_test.csv
            cat solution_1_1.csv >> ../training_data/solution_1_1_test.csv
          done
        done
      done
    done
  done
done

