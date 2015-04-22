NAME_OF_EXECUTABLE_FILE = motif

COMPILER = g++ 
WARNINGS = -Wchar-subscripts -Wparentheses -Wreturn-type -Wmissing-braces -Wundef -Wshadow
COMPILER_OPTS = -c -g -O0 $(WARNINGS)
 
LINKER = g++ 
LINKER_OPTS = -o $(NAME_OF_EXECUTABLE_FILE)
 
OBJS = benchmark.o motif.o main.o

$(NAME_OF_EXECUTABLE_FILE): $(OBJS) 
	$(LINKER) $(OBJS) $(LINKER_OPTS)

main.o: main.cpp motif.h benchmark.h
	$(COMPILER) $(COMPILER_OPTS) main.cpp

benchmark.o: benchmark.cpp benchmark.h
	$(COMPILER) $(COMPILER_OPTS) benchmark.cpp

motif.o: motif.cpp motif.h
	$(COMPILER) $(COMPILER_OPTS) motif.cpp 

clean:
	-rm -f *.o $(NAME_OF_EXECUTABLE_FILE)

tidy:
	-rm -rf ./doc