CXX = g++

EXEC = ml-classification.exe 

CG_learn =  nn-class.o

ifdef MAC_OBLAS  			# macbook with openblas installed using macports
	CXXFLAGS = -std=c++11 -O2 -fopenmp -Wall -DMAC_OBLAS
	INCLUDE = -I/opt/local/include
	LDLIBS = -L/opt/local/lib -lopenblas -lpthread -lm -ldl
else ifdef OBLAS             		# linux with openblas
	CXXFLAGS = -std=c++11 -O2 -fopenmp -Wall -DOBLAS
	INCLUDE = -I${BLAS_INC}
	LDLIBS = -L${BLAS_LIB} -lopenblas -lpthread -lm -ldl
else           				# linux with mkl (default)

	CXXFLAGS = -std=c++11 -O2 -fopenmp -Wall
	INCLUDE  = -DMKL_ILP64 -m64 -I${MKLROOT}/include
	LDLIBS   = -L${MKLROOT}/lib/intel64 -lmkl_scalapack_ilp64 -lmkl_cdft_core -lmkl_intel_ilp64 -lmkl_core -lmkl_sequential -lmkl_blacs_openmpi_ilp64 -lpthread -lm -ldl
	LDLIBS  += -L${MPI_LIB} -lmpi
endif



# MAC with intel cpu suggestion for mkl linking
#dynamic
#ONEAPI_ADD = ${ONEAPI_ROOT}/compiler/2023.0.0/mac/compiler/lib
#LDLIBS = -L${MKLROOT}/lib -lmkl_scalapack_ilp64 -lmkl_intel_ilp64 -lmkl_intel_thread -lmkl_core -lmkl_blacs_mpich_ilp64 -lpthread -lm -ldl -L${ONEAPI_ADD} -liomp5
#LDLIBS = -L/opt/local/lib -lopenblas -lpthread -lm -ldl

install : $(CG_learn)   
	$(CXX) $(CG_learn) $(CXXFLAGS) $(LDLIBS) -o ./$(EXEC)

clean:
	rm -f *.o   
	rm -f $(EXEC)


%.o: %.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDE) -c  $(<) $(outputs)



