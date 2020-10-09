G14DIS
======

## An example compile and run command on Windows:
```
g++ common.cpp element.cpp linearSystems.cpp matrix.cpp matrix_full.cpp mesh.cpp refinement.cpp quadrature.cpp solution.cpp solution_linear.cpp solution_nonlinear.cpp ../tests/hp_sin.cpp -o ../build/hp_sin.exe; ../build/hp_sin.exe
```

## Modified compile and run for dg example working on 2020-10-06
```
g++ common.cpp element.cpp linearSystems.cpp matrix.cpp matrix_full.cpp mesh.cpp refinement.cpp quadrature.cpp solution.cpp solution_linear.cpp solution_nonlinear.cpp solution_dg_linear.cpp ../tests/h_dg_boundary.cpp -o ../build/h_dg_boundary.exe; ../build/h_dg_boundary.exe
```