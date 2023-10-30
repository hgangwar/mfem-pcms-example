# This project is part of the Fundamentals of Finite Element Course.

## How to run the project and visulaize the results

### Visualize and varify the mesh before running the code on glvis

```
.\glvis-Windows\glvis.exe -m cylfuelcell3d-v2.msh
```

### Run and visualize the results on glvis

1. Get and build (with parallel option enabled) [mfem](https://mfem.org/building/).
2. Go to the mfem examples directory

```
cd mfem/examples
```

3. Clone this repository

```
git clone https://github.com/Fuad-HH/FEMProject.git
```

4. Build example 11 (or your code).

```
make ex11p
```

5. Clone the repository

```
git clone https://github.com/Fuad-HH/FEMProject.git
```

6. Run the example

```
mpirun -np 4 ./ex11p -m ./FEMProject/Mesh/cylfuelcell3d-v2.msh -rs 0 -rp 0 -n 1
```

7. To visualize the result on glvis

```
.\glvis-Windows\glvis.exe -np 4 -m mesh -g mode_00
```

## Project Description

We are simulating a fuel cell using finite element method, solving neutron flux and temperature with two different mfem solvers and coupling them.

## Project Structure
