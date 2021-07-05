
OBJECTS = vars.o singleton.o maths.o allen.o svd.o io_py.o SEE.o rt_coef.o synth.o hazel_py.omarquardt.o hazel.o DIRect.o \
		DIRserial.o DIRsubrout.o singleton.o mpi_routines.o inversion.o	
	
hazel: $(OBJECTS)
	
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
