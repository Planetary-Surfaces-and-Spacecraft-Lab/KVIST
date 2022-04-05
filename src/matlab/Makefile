all: compile_mex.m
	matlab -nodesktop -nosplash -nodisplay -r 'compile_mex; exit'
	echo "Done compiling"
	if [ -d "../../lib" ]; then rm -rf ../../lib; fi
	mkdir ../../lib
	cp *.m *.mex* ../../lib

.PHONY : clean
clean :  
	rm -rf codegen
	rm *.mex*    
	rm -rf ../../lib
