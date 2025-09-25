# Compiler and flags
CC = g++ 
CPPFLAGS = -g -std=c++11 -Wall -O3
LDFLAGS =


# Common source files (without mains)
COMMON_SRC = utilities.cpp TSP.cpp TSPSolution.cpp localSearch.cpp genetic_algorithm.cpp
COMMON_OBJ = $(COMMON_SRC:.cpp=.o)

# Executable-specific source files
MAIN_LS_SRC = main_ls.cpp
GENETIC_SRC = main_genetic.cpp

# Targets
EXECUTABLES = localSearch genetic

# Default target
all: $(EXECUTABLES)

# Rules for each executable
localSearch: $(COMMON_OBJ) $(MAIN_LS_SRC:.cpp=.o)
	$(CC) $(CPPFLAGS) $^ -o $@ 

genetic: $(COMMON_OBJ) $(GENETIC_SRC:.cpp=.o)
	$(CC) $(CPPFLAGS) $^ -o $@ 

%.o: %.cpp
	$(CC) $(CPPFLAGS) -c $< -o $@

# Clean target
clean:
	rm -rf $(COMMON_OBJ) $(MAIN_LS_SRC:.cpp=.o) $(GENETIC_SRC:.cpp=.o) $(EXECUTABLES)

.PHONY: all clean