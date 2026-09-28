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



---

## 2. Fourier Coefficient Calculation for a C5 Symmetry Sequence

### Overview

This Fortran program calculates complex Fourier coefficients associated with the I+ and I- operators for a C5 symmetry pulse sequence in solid-state NMR.

### Source Code

`C5_Fourier.f`

### Sequence Parameters

- **Symmetry sequence:** C5
- **Number of pulse elements:** 10
- **Fourier index for I+:** n = -2
- **Fourier index for I-:** n = -2

### Method

The program defines the phases and angular intervals of the ten pulse elements and evaluates the corresponding complex integrals for the I+ and I- operators.

The individual contributions are combined to obtain the final complex Fourier coefficients.

### Output

The program prints:

- The sequence identifier (C5)
- The Fourier index for each operator
- The real and imaginary parts of the coefficient associated with I+
- The real and imaginary parts of the coefficient associated with I-

### Requirements

- A Fortran compiler, such as GNU Fortran (gfortran).

### Compilation

```bash
gfortran C5_Fourier.f -o C5_Fourier
```

### Execution

```bash
./C5_Fourier
```

### Purpose

This project demonstrates analytical and numerical calculations of Fourier coefficients for symmetry-based pulse sequences in solid-state NMR.

This project demonstrates scientific programming and data processing for solid-state NMR applications.

