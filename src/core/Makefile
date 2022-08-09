CPP_FLAGS=-O3 -std=c++17
C_FLAGS=-O3
BLAS_INC=-I/usr/local/opt/openblas/include
BLAS_LINK=-lcblas -L/usr/local/opt/openblas/lib

all: test_math_util test_ist_spectra
	
test_math_util: test_math_util.o math_util.o
	clang++ ${CPP_FLAGS} test_math_util.o math_util.o -o test_math_util

test_ist_spectra: math_util.o ist_spectra.o mmat_automatic.o mmat_cpp_safe_binding.o
	clang++ ${CPP_FLAGS} math_util.o ist_spectra.o mmat_automatic.c mmat_cpp_safe_binding.cpp test_ist_spectra.cpp ${BLAS_INC} ${BLAS_LINK} -o test_ist_spectra

ist_spectra.o: math_util.o interval_search.o mmat_cpp_safe_binding.o
	clang++ ${CPP_FLAGS} -c interval_search.o ist_spectra.cpp mmat_cpp_safe_binding.o math_util.o -o ist_spectra.o

mmat_automatic.o:
	clang ${C_FLAGS} -c mmat_automatic.c ${BLAS_INC} -o mmat_automatic.o

interval_search.o:
	clang++ ${CPP_FLAGS} -c interval_search.cpp -o interval_search.o

mmat_cpp_safe_binding.o: mmat_automatic.o
	clang++ ${CPP_FLAGS} -c mmat_cpp_safe_binding.cpp mmat_automatic.o ${BLAS_INC} -o mmat_cpp_safe_binding.o

test_math_util.o:
	clang++ ${CPP_FLAGS} -c test_math_util.cpp -o test_math_util.o

math_util.o:
	clang++ ${CPP_FLAGS} -c math_util.cpp

.PHONY : clean
clean :  
	rm *.o
