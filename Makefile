CXX = g++
CXXFLAGS = -Wall -Wextra -std=c++17
LDFLAGS = -lm -lhts -lpthread -llzma -lz -lbz2
INCLUDES = -I/r1/people/bioinf/src/htslib-1.4/htslib
LIBS = -L/r1/people/bioinf/src/htslib-1.4/htslib

analyzeBAM = analyzeBAM
filterBAM = filterBAM
deamBAM = deamBAM
analyzeVCF = analyzeVCF

# Source files
SRCS1 = analyzeBAM.cpp
SRCS2 = filterBAM.cpp
SRCS3 = deamBAM.cpp
SRCS4 = analyzeVCF.cpp

# Object files
OBJS1 = $(SRCS1:.cpp=.o)
OBJS2 = $(SRCS2:.cpp=.o)
OBJS3 = $(SRCS3:.cpp=.o)
OBJS4 = $(SRCS4:.cpp=.o)

.PHONY: all clean

all: $(analyzeBAM) $(filterBAM) $(deamBAM) $(analyzeVCF)

$(analyzeBAM): $(OBJS1) 
	$(CXX) $(CXXFLAGS) $(INCLUDES) -o $(analyzeBAM) $(OBJS1) $(LDFLAGS) $(LIBS)

$(filterBAM): $(OBJS2) 
	$(CXX) $(CXXFLAGS) $(INCLUDES) -o $(filterBAM) $(OBJS2) $(LDFLAGS) $(LIBS)

$(deamBAM): $(OBJS3)
	$(CXX) $(CXXFLAGS) $(INCLUDES) -o $(deamBAM) $(OBJS3) $(LDFLAGS) $(LIBS)

$(analyzeVCF): $(OBJS4)
	$(CXX) $(CXXFLAGS) $(INCLUDES) -o $(analyzeVCF) $(OBJS4) $(LDFLAGS) $(LIBS)

analyzeBAM: analyzeBAM.cpp
filterBAM: filterBAM.cpp
deamBAM: deamBAM.cpp
analyzeVCF: analyzeVCF.cpp

.cpp.o:
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

clean:
	rm -f *.o
