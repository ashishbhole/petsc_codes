FC      = mpif90
GIT_VERSION = $(shell git describe --abbrev=8 --dirty --always)
FFLAGS  = -O3 -fdefault-real-8 -I$(PETSC_DIR)/arch-linux-c-debug/include \
		-DVERSION=\"$(GIT_VERSION)\"
LDFLAGS = -lpetsc -L$(PETSC_DIR)/arch-linux-c-debug/lib -L/usr/lib64

SRC  = $(wildcard *.F90)
OBJ  = $(patsubst %.F90,%.o,$(SRC))

TARGET  = exe

$(TARGET): $(OBJ) #$(MOD) $(OMOD)
	$(FC) -o $(TARGET) $(OBJ) $(LDFLAGS)

%.o: %.F90
	$(FC) $(FFLAGS) -c -o $@ $<

clean:
	rm -f *.mod *.o $(OBJ) $(TARGET)
