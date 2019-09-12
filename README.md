# PDE Numerical Solver
Numerical methods to solve PDE problems in space-time mesh.

## Dependences
  * Codes made in Fortran 90
  * gfortran Compiler
  * Plot solution using Matlab/Octave
  
## Mathematical Model
  In this code, we solve a bidimensional partial equation in space-time using a finite differences algorithm. 
  
  ![equation 1](https://latex.codecogs.com/gif.latex?%5Cnabla%20%5Ccdot%20%5Cleft%28%20k%20%5Cnabla%20u%20%5Cright%29%20&plus;%20%5Cnabla%28vu%29%20&plus;%20%5Csigma%20u%20%3D%20f%28x%2Cy%29%20&plus;%20%5Cdisplaystyle%20%5Cfrac%7B%5Cpartial%20u%7D%7B%5Cpartial%20t%7D)
 
## How to use
### Running code
  The code runs in Fortran 90, you will need a fortran compiler, such as gfortran.
  The problem conditions are changed in code, then you'll need to compile every change:
  
  ```
    gfortran 2d_transient.f90 -o <your_exe_name>
  ```
  
  Then, run:
  
  * On Windows
  ```
    your_exe_name.exe
  ```
  * On Linux 
  ```
    ./your_exe_name.out
  ```
  
 After this, the code will generate three .out files.
 * output_mash.out : 
 Contain the points of domain discretization. 
 * output_solution.out :
 Contain the solution for each point in each time step
 * plot_info.out :
 Contain important information about the mesh and time step that will be used to plot solutions.
 
 ### Plotting solution
You will need Matlab or Octave to run the .m code. 
     Once with Matlab/Octave opened, just run the code with execute button 
     and watch the solutin change in time.
