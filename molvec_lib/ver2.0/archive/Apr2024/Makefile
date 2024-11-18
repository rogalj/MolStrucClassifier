
CC=g++
CX=mpicxx
#CXX = g++
#CC  = gcc 
all: 
#	$(CX) -O3 -lm -Wall -o sa readinput.cpp symmetryfunc.cpp nn.cpp steinhardt.cpp main.cpp 
#$(CC) -O3 -lm -Wall -o sa readinput.cpp main.cpp symmetryfunc.cpp molvec.cpp nn.cpp 
	$(CC) -O3 -lm -Wall -o sa readinput.cpp main.cpp  molvec.cpp nn.cpp 
#	$(CX) -O3 -lm -Wall -o sa readinput.cpp main.cpp  molvec.cpp nn.cpp 
#	$(CX) -O3 -lm -Wall -o sa nn.cpp
#	$(CC) -O3 -lm -Wall -fsanitize=address -fno-omit-frame-pointer -o sa_test readinput.cpp main.cpp  molvec.cpp nn.cpp
#	$(CC) -O3 -lm -Wall -fsanitize=leak -fno-omit-frame-pointer -o sa_test readinput.cpp main.cpp  molvec.cpp nn.cpp

#
gprof:
	$(CC) -O3 -mtune=native -march=native -mfpmath=both -lm -Wall -pg -o sa_test readinput.cpp main.cpp  molvec.cpp nn.cpp

#lib:
#	MAC
#	$(CX) -O3 -lm -Wall -dynamiclib -o libsa.dylib readinput.cpp symmetryfunc.cpp nn.cpp main.cpp 
#	LINUX
#	$(CX) -O3 -lm -Wall -fPIC -c readinput.cpp symmetryfunc.cpp nn.cpp steinhardt.cpp main.cpp 
#	$(CX) -O3 -lm -Wall -shared -o libsa.so readinput.o symmetryfunc.o nn.o steinhardt.o main.o

libggmtlm:
#       LINUX
	 $(CX) -O3 -lm -Wall -fPIC -c readinput.cpp symmetryfunc.cpp nn.cpp steinhardt.cpp main.cpp 
	 #$(CX) -std=gnu++98 -O3 -lm -Wall -fPIC -c readinput.cpp symmetryfunc.cpp nn.cpp steinhardt.cpp main.cpp 
	 $(CX) -O3 -lm -Wall -shared -o libcv_nn_lm.so readinput.o symmetryfunc.o nn.o steinhardt.o main.o
	 #$(CX) -std=gnu++98 -O3 -lm -Wall -shared -o libsa-ggmt.so readinput.o symmetryfunc.o nn.o steinhardt.o main.o

#libggmt:
lib:	 
#       LINUX
#	$(CX) -O3 -Wall -fPIC -c readinput.cpp symmetryfunc.cpp nn.cpp molvec.cpp main.cpp
	$(CX) -O3 -Wall -fPIC -c readinput.cpp nn.cpp molvec.cpp main.cpp
	#$(CX) -std=gnu++98 -O3 -Wall -fPIC -c readinput.cpp symmetryfunc.cpp nn.cpp steinhardt.cpp main.cpp
	$(CX) -O3 -Wall -shared -o libcv_nn.so readinput.o molvec.o nn.o  main.o
	#$(CX) -std=gnu++98 -O3 -Wall -shared -o libsa-ggmt.so readinput.o symmetryfunc.o nn.o steinhardt.o main.o

lib_mac:
#	MAC
	$(CX) -O3 -lm -Wall -dynamiclib -o libcv_nn.dylib readinput.cpp symmetryfunc.cpp nn.cpp steinhardt.cpp main.cpp

staticlib:
	$(CX) -O3 -Wall -c readinput.cpp symmetryfunc.cpp nn.cpp steinhardt.cpp main.cpp
	ar rc orderggmt.a readinput.o symmetryfunc.o nn.o steinhardt.o main.o
	ranlib orderggmt.a

print:
	$(CX) -O3 -lm -Wall -o sa_print readinput.cpp symmetryfunc.cpp nn.cpp print.cpp histogram.cpp steinhardt.cpp main-printtraj.cpp

clean:
#	rm -f *.o sa sa_print libsa.dylib libsa.so
	rm -f *.o sa sa_print 
