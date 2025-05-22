CXX = g++
CXXFLAGS = -Wall -Wextra -std=c++17 -O3
LDFLAGS = -lm -lhts -lm -lpthread -llzma -lz -lbz2
INCLUDES = -I/r1/people/bioinf/src/htslib-1.4/htslib
LIBS = -L/r1/people/bioinf/src/htslib-1.4/htslib

MAIN1 = analyzeBAM
MAIN2 = filterBAM
MAIN3 = deamBAM

# Source files
SRCS1 = analyzeBAM.cpp
SRCS2 = filterBAM.cpp
SRCS3 = deamBAM.cpp

# Object files
OBJS1 = $(SRCS1:.cpp=.o)
OBJS2 = $(SRCS2:.cpp=.o)
OBJS3 = $(SRCS3:.cpp=.o)

.PHONY: all clean

all: $(MAIN1) $(MAIN2) $(MAIN3)

$(MAIN1): $(OBJS1) 
	$(CXX) $(CXXFLAGS) $(INCLUDES) -o $(MAIN1) $(OBJS1) $(LDFLAGS) $(LIBS)

$(MAIN2): $(OBJS2) 
	$(CXX) $(CXXFLAGS) $(INCLUDES) -o $(MAIN2) $(OBJS2) $(LDFLAGS) $(LIBS)

$(MAIN3): $(OBJS3)
	$(CXX) $(CXXFLAGS) $(INCLUDES) -o $(MAIN3) $(OBJS3) $(LDFLAGS) $(LIBS)

analyzeBAM: analyzeBAM.cpp
filterBAM: filterBAM.cpp
deamBAM: deamBAM.cpp

.cpp.o:
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $<  -o $@

clean:
	rm -f *.o *~ $(MAIN1) $(MAIN2) $(MAIN3)
