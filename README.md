# nmr-data-analysis
A collection of Fortran programs for scientific data processing, including time-unit conversion for solid-state NMR data.



# NMR Data Analysis

A collection of Fortran programs for scientific data processing, including time-unit conversion for solid-state NMR data.

## 1. Time Conversion

### Overview

This Fortran program converts the time values in an NMR data file from seconds to milliseconds while preserving the real and imaginary data columns.

### Source Code

`Time converter (s to ms).f`

### Input

The program reads a file named `fid.dat` containing three columns:

1. Time (seconds)
2. Real part
3. Imaginary part

### Processing

The first column is multiplied by 1000 to convert the time values from seconds to milliseconds. The real and imaginary parts are retained.

### Output

The converted data are written to `fid_new.dat`.

The output contains the time in fixed-point notation with six decimal places and the real and imaginary parts in scientific notation.

### Requirements

- GNU Fortran (gfortran)

### Compilation

```bash
gfortran "Time converter (s to ms).f" -o time_converter
```

### Execution

Place `fid.dat` in the working directory and run:

```bash
./time_converter
```

The program generates `fid_new.dat` in the same directory.

### Purpose

This project demonstrates scientific programming and data processing for solid-state NMR applications.

