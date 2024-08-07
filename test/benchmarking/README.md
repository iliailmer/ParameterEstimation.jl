# Benchmarking Results

This directory contains the files used to generate the benchmarking results in our paper.

## Directory Structure

The directory `test_files` contains subdirectories corresponding to the four software platforms:
- `AMIGO2`
- `IQM`
- `ParameterEstimation`
- `SciML`

## File Naming Convention

The test files are named in the following format:
[DIFFERENTIAL MODEL NAME]_[TEST NUMBER]_[UPPER BOUND FOR SEARCH INTERVAL]

## Running the Tests

### AMIGO2 and IQM

To run AMIGO2 and IQM, make sure to change the `AMIGO2_PATH` and `IQM_PATH` in the scripts to your own path to the packages.

### Running All Tests

The file `run.sh` in each folder runs all tests in the system. The script captures the outputs and stores them as `.out` files in the `outputs` subdirectory.
