.SUFFIXES:
FC      = mpif90
GIT_VERSION = $(shell git describe --abbrev=8 --dirty --always)
FFLAGS  = -O3 -fdefault-real-8 -I$(PETSC_DIR)/arch-linux-c-debug/include \
		-DVERSION=\"$(GIT_VERSION)\"
LDFLAGS = -lpetsc -L$(PETSC_DIR)/arch-linux-c-debug/lib -L/usr/lib64

SRC  = $(wildcard *.F90)
OBJ  = $(patsubst %.F90,%.o,$(SRC))
MOD  = double.mod variables.mod

TARGET  = exe

$(TARGET): $(OBJ) $(MOD)
	$(FC) -o $(TARGET) $(OBJ) $(LDFLAGS)

double.mod: double.F90
	$(FC) -c $(FFLAGS) double.F90

variables.mod: variables.F90 double.mod
	$(FC) -c $(FFLAGS) variables.F90

%.o: %.F90 $(MOD)
	$(FC) $(FFLAGS) -c -o $@ $<

clean:
	rm -f *.mod *.o $(OBJ) $(TARGET)
