all: compile_mex.m
	matlab -nodesktop -nosplash -nodisplay -r 'compile_mex; exit'
	echo "Done compiling"
	if [ -d "../../bin" ]; then rm -rf ../../bin; fi
	mkdir ../../bin
	cp *.m *.mex* ../../bin

.PHONY : clean
clean :  
	rm -rf codegen
	rm *.mex*    
	rm -rf ../../bin
