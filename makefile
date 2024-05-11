

PETSC_ARCH ?= complex-opt

SRCDIR = ./SourceCode

SRCS := $(wildcard $(SRCDIR)/*.c)

TARGET = RSVDt_clean

CC = mpicc

CFLAGS = -fPIC -Wall -Wwrite-strings -Wno-unknown-pragmas -Wno-lto-type-mismatch -fstack-protector -fvisibility=hidden -g -O
CPPFLAGS = -I${SLEPC_DIR}/include -I${SLEPC_DIR}/$(PETSC_ARCH)/include -I${PETSC_DIR}/include -I${PETSC_DIR}/$(PETSC_ARCH)/include -I$(SRCDIR)
LDFLAGS = -Wl,-export-dynamic -Wl,-rpath,${SLEPC_DIR}/$(PETSC_ARCH)/lib -L${SLEPC_DIR}/$(PETSC_ARCH)/lib -lslepc -Wl,-rpath,${PETSC_DIR}/$(PETSC_ARCH)/lib -L${PETSC_DIR}/$(PETSC_ARCH)/lib -lpetsc -lpthread -lscalapack -lflapack -lfblas -lm -lX11 -ldl -lmpi_usempif08 -lmpi_usempi_ignore_tkr -lmpi_mpifh -lmpi -lgfortran -lm -lgfortran -lm -lgcc_s -lquadmath -lpthread -lstdc++ -lquadmath -ldl

$(TARGET): $(SRCS)
	$(CC) $(CFLAGS) $(CPPFLAGS) $^ $(LDFLAGS) -o $@

clean:
	rm -f $(TARGET)
