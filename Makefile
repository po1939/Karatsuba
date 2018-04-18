PROGS=driver
HEADERS=posint.hpp
CPPFLAGS=-O3 -Wall -Wno-sign-compare -Wno-unused-function
#CPPFLAGS=-Wall -Wextra -Wno-sign-compare -fprofile-arcs -ftest-coverage -g

# Default target
all: $(PROGS)

# Dependencies
$(PROGS): $(HEADERS:.hpp=.o)

# Rules to generate the final compiled program
$(PROGS): %: %.cpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -o $@ $^

# Generic rule for compiling C++ programs from source
# (Actually, make also defines this by default.)
%.o: %.cpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $<

.PHONY: clean all
clean:
	rm -f *.o $(PROGS)
