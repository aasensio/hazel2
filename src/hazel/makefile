# makefile for Hazel

# Choose MPI or serial version (version=mpi/serial/python)
ifndef $(version)
	version = serial
endif

# Choose desired compiler (compiler=gfortran/ifort)
ifndef $(compiler)
	compiler = ifort
endif

# Choose desired compiler (supercomputer=local/teide)
ifndef $(supercomputer)
	supercomputer = local
endif


#******************
# COMPILATION OPTIONS
# Add your options below
#******************
# In my desktop computer
ifeq ($(version), mpi)
	
#------------------------
# My local supercomputer
#------------------------
	ifeq ($(supercomputer), local)		
		ifeq ($(compiler), gfortran)
			NETCDF_INCLUDE = -I/usr/lib64/gfortran/modules
			NETCDF_LIB = -L/usr/pkg/netCDF/netCDF-4.1.1/lib -lnetcdff
			COMPILER = mpif90
			COMPILER_LINKS = #-L/net/nas4/fis/aasensio/scratch/MPICH/lib 
			COMPILER_OPTS = -c -O3 -ffree-line-length-none #-I/net/nas4/fis/aasensio/scratch/MPICH/include
			PREPROCESS = -D$(version) -cpp
		else
			NETCDF_INCLUDE = -I/net/nas4/fis/aasensio/scratch/LIBS/NETCDF/include
			NETCDF_LIB = -L/net/nas4/fis/aasensio/scratch/LIBS/NETCDF/lib -lnetcdf
			COMPILER = /net/nas4/fis/aasensio/scratch/MPICH/bin/mpif90
			COMPILER_LINKS = -L/net/nas4/fis/aasensio/scratch/MPICH/lib 
			COMPILER_OPTS = -c -O3 -I/net/nas4/fis/aasensio/scratch/MPICH/include #-ipo
			PREPROCESS = -D$(version) -fpp
		endif
	endif

#------------------------
# Teide supercomputer
#------------------------
	ifeq ($(supercomputer), teide)
		NETCDF_INCLUDE = -I/opt/envhpc/netcdf/4.3.2/intelmpi/include
		NETCDF_LIB = -L/opt/envhpc/netcdf/4.3.2/intelmpi/lib -lnetcdf
		COMPILER = mpiifort		 
		COMPILER_OPTS = -c -O3 -ipo -xAVX
		PREPROCESS = -D$(version) -fpp		
	endif

	HAZEL_DEP = hazel.f90 allen.o synth.o rt_coef.o SEE.o maths.o vars.o mpi_routines.o inversion.o io.o
	HAZEL_MAIN_OBJ = hazel.o
	HAZEL_MAIN_SRC = hazel.f90
	OBJECTS = vars.o maths.o allen.o svd.o io.o SEE.o rt_coef.o synth.o marquardt.o hazel.o DIRect.o \
		DIRserial.o DIRsubrout.o singleton.o mpi_routines.o inversion.o	

endif

ifeq ($(version), serial)
	ifeq ($(compiler), gfortran)
		COMPILER = gfortran
		COMPILER_OPTS = -c -O3 -ffree-line-length-none 
		PREPROCESS = -D$(version) -cpp
	else		
		COMPILER = ifort
		COMPILER_OPTS = -c -O3 
		PREPROCESS = -D$(version) -fpp
	endif
	HAZEL_DEP = hazel.f90 allen.o synth.o rt_coef.o SEE.o maths.o vars.o mpi_routines_fake.o inversion.o io.o
	HAZEL_MAIN_OBJ = hazel.o
	HAZEL_MAIN_SRC = hazel.f90
	OBJECTS = vars.o maths.o allen.o svd.o io.o SEE.o rt_coef.o synth.o marquardt.o hazel.o DIRect.o \
		DIRserial.o DIRsubrout.o singleton.o mpi_routines_fake.o inversion.o
endif

ifeq ($(version), python)
	COMPILER = gfortran
    COMPILER_OPTS = -c -ffree-line-length-none -fPIC -O3 -fno-automatic -ffast-math -funroll-loops
    HAZEL_DEP = hazel_py.f90 allen.o synth.o rt_coef.o SEE.o maths.o vars.o io_py.o
    HAZEL_MAIN_OBJ = hazel_py.o
    HAZEL_MAIN_SRC = hazel_py.f90
    OBJECTS = vars.o maths.o allen.o svd.o io_py.o SEE.o rt_coef.o synth.o hazel_py.o singleton.o
    PREPROCESS = -D$(version) -cpp
endif


#******************
# Hazel
#******************
LIBS = $(NETCDF_LIB) -lcurl
FINAL_EXECUTABLE = hazel

print-% : ; $(info $* is $(flavor $*) variable set to [$($*)]) @true
	
hazel: $(OBJECTS)
ifeq ($(version), serial)
	$(COMPILER) $(COMPILER_LINKS) $(OBJECTS) $(LIBS) $(NETCDF_INCLUDE)  -o $(FINAL_EXECUTABLE)
	cp $(FINAL_EXECUTABLE) ../run
	cp $(FINAL_EXECUTABLE) ../synthGUI
	cp $(FINAL_EXECUTABLE) ../invGUI
endif

ifeq ($(version), mpi)
	$(COMPILER) $(COMPILER_LINKS) $(OBJECTS) $(NETCDF_INCLUDE) $(LIBS) -o $(FINAL_EXECUTABLE)
	cp $(FINAL_EXECUTABLE) ../runMPI	
endif
		
clean: 
	find . -maxdepth 2 -name "*.o" -delete ; find . -maxdepth 1 -name "*.mod" -delete ;
	find . -maxdepth 1 -name "*.f90~" -delete ; find . -maxdepth 2 -name "*.a" -delete ;
	find . -maxdepth 1 -name "hazel" -delete ; find . -maxdepth 1 -name "*.so" -delete; 
	find . -maxdepth 1 -name "pyhazel.*" -delete;
	find ../ -maxdepth 2 -name "hazel" -delete 
	
vars.o: vars.f90
	$(COMPILER) $(COMPILER_OPTS)  vars.f90

maths.o: maths.f90 vars.o singleton.o
	$(COMPILER) $(COMPILER_OPTS)  maths.f90

allen.o: allen.f90 maths.o vars.o
	$(COMPILER) $(COMPILER_OPTS)  allen.f90
	
svd.o: svd.f90
	$(COMPILER) $(COMPILER_OPTS)  svd.f90

io.o: io.f90 maths.o vars.o
	$(COMPILER) $(COMPILER_OPTS) $(NETCDF_INCLUDE) $(PREPROCESS) io.f90

io_py.o: io_py.f90 maths.o vars.o
	$(COMPILER) $(COMPILER_OPTS) io_py.f90

SEE.o: SEE.f90 allen.o maths.o vars.o
	$(COMPILER) $(COMPILER_OPTS)  SEE.f90

rt_coef.o: rt_coef.f90 vars.o maths.o
	$(COMPILER) $(COMPILER_OPTS) $(PREPROCESS) rt_coef.f90

synth.o: synth.f90 rt_coef.o SEE.o vars.o
	$(COMPILER) $(COMPILER_OPTS)  synth.f90

marquardt.o: marquardt.f90 maths.o svd.o synth.o rt_coef.o SEE.o vars.o
	$(COMPILER) $(COMPILER_OPTS) $(NETCDF_INCLUDE) marquardt.f90
	
DIRect.o: DIRect.f90
	$(COMPILER) $(COMPILER_OPTS)  DIRect.f90
	
DIRserial.o: DIRserial.f90
	$(COMPILER) $(COMPILER_OPTS)  DIRserial.f90
	
DIRsubrout.o: DIRsubrout.f90
	$(COMPILER) $(COMPILER_OPTS)  DIRsubrout.f90

singleton.o: singleton.f90
	$(COMPILER) $(COMPILER_OPTS)  singleton.f90

mpi_routines.o: mpi_routines.f90
	$(COMPILER) $(COMPILER_OPTS) mpi_routines.f90

mpi_routines_fake.o: mpi_routines_fake.f90
	$(COMPILER) $(COMPILER_OPTS) mpi_routines_fake.f90

inversion.o: inversion.f90 marquardt.o
	$(COMPILER) $(COMPILER_OPTS) $(NETCDF_INCLUDE) inversion.f90
	
$(HAZEL_MAIN_OBJ): $(HAZEL_DEP)
	$(COMPILER) $(COMPILER_OPTS) $(NETCDF_INCLUDE) $(PREPROCESS) $(HAZEL_MAIN_SRC)
