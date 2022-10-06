Some experimental codes in FORTRAN using PETSc.
To use these code one should install PETSc and set environmental variables:
```
PETSC_DIR
PETSC_ARCH
LD_LIBRARY_PATH
```

Following configuration is used for building PETSc 3.18.0 is:
```
./configure --with-cc=mpicc --with-cxx=mpicxx --with-fc=mpif90 \
--download-fblaslapack --download-scalapack --download-mumps --download-superlu --download-ptscotch \
--with-metis-include=/user/abhole/home/lib/metis-5.1.0/include \
--with-metis-lib=/user/abhole/home/lib/metis-5.1.0/lib/libmetis.a -lmetis \
--with-parmetis-include=/user/abhole/home/lib/parmetis-4.0.3/include \
--with-parmetis-lib=/user/abhole/home/lib/parmetis-4.0.3/lib/libparmetis.a -lparmetis -lmetis \
--with-hdf5-include=/user/abhole/home/lib/hdf5-1.8.18/include \
--with-hdf5-lib=/user/abhole/home/lib/hdf5-1.8.18/lib64/libhdf5.a \
--with-valgrind=1 \
--with-scalar-type=real \
--with-precision=double
```

