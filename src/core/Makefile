CPPFLAGS=-O3 -std=c++17
CFLAGS=-O3
BLAS_INC=-I/usr/local/opt/openblas/include
BLAS_LINK=-lcblas -L/usr/local/opt/openblas/lib

all: test_math_util test_ist_spectra
	
test_math_util: test_math_util.o math_util.o
	clang++ ${CPPFLAGS} test_math_util.o math_util.o -o test_math_util

test_ist_spectra: math_util.o ist_spectra.o mmat_automatic.o mmat_cpp_safe_binding.o test_ist_spectra.cpp
	clang++ ${CPPFLAGS} interval_search.o math_util.o mmat_automatic.o mmat_cpp_safe_binding.o ist_spectra.o test_ist_spectra.cpp ${BLAS_INC} ${BLAS_LINK} -o test_ist_spectra

ist_spectra.o: math_util.o interval_search.o mmat_cpp_safe_binding.o mmat_automatic.o ist_spectra.cpp
	clang++ ${CPPFLAGS} -c interval_search.o ist_spectra.cpp mmat_automatic.o mmat_cpp_safe_binding.o math_util.o -o ist_spectra.o

mmat_automatic.o: mmat_automatic.c
	clang ${CFLAGS} -c mmat_automatic.c ${BLAS_INC} -o mmat_automatic.o

interval_search.o: interval_search.cpp mmat_automatic.o
	clang++ ${CPPFLAGS} -c interval_search.cpp mmat_automatic.o -o interval_search.o

mmat_cpp_safe_binding.o: mmat_automatic.o mmat_cpp_safe_binding.cpp
	clang++ ${CPPFLAGS} -c mmat_cpp_safe_binding.cpp mmat_automatic.o ${BLAS_INC} -o mmat_cpp_safe_binding.o

test_math_util.o: test_math_util.cpp
	clang++ ${CPPFLAGS} -c test_math_util.cpp -o test_math_util.o

math_util.o: math_util.cpp
	clang++ ${CPPFLAGS} -c math_util.cpp

.PHONY : clean
clean :  
	rm *.o
