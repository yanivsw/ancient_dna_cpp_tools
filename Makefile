CXX = g++
CXXFLAGS = -Wall -Wextra -std=c++17
LDFLAGS = -lm -lhts -lpthread -llzma -lz -lbz2
INCLUDES = -I/r1/people/bioinf/src/htslib-1.4/htslib -Iinclude
LIBS = -L/r1/people/bioinf/src/htslib-1.4/htslib

# Directories
SRCDIR = src
INCDIR = include
OBJDIR = obj

# Executables
analyzeBAM = analyzeBAM
filterBAM = filterBAM
deamBAM = deamBAM
analyzeVCF = analyzeVCF

# Source files
SRCS1 = $(SRCDIR)/analyzeBAM.cpp
SRCS2 = $(SRCDIR)/filterBAM.cpp
SRCS3 = $(SRCDIR)/deamBAM.cpp
SRCS4 = $(SRCDIR)/analyzeVCF.cpp

# Object files
OBJS1 = $(OBJDIR)/analyzeBAM.o
OBJS2 = $(OBJDIR)/filterBAM.o
OBJS3 = $(OBJDIR)/deamBAM.o
OBJS4 = $(OBJDIR)/analyzeVCF.o

.PHONY: all clean

all: $(analyzeBAM) $(filterBAM) $(deamBAM) $(analyzeVCF)

# Create object directory if it doesn't exist
$(OBJDIR):
	mkdir -p $(OBJDIR)

$(analyzeBAM): $(OBJDIR) $(OBJS1) 
	$(CXX) $(CXXFLAGS) -o $(analyzeBAM) $(OBJS1) $(LIBS) $(LDFLAGS)

$(filterBAM): $(OBJDIR) $(OBJS2) 
	$(CXX) $(CXXFLAGS) -o $(filterBAM) $(OBJS2) $(LIBS) $(LDFLAGS)

$(deamBAM): $(OBJDIR) $(OBJS3)
	$(CXX) $(CXXFLAGS) -o $(deamBAM) $(OBJS3) $(LIBS) $(LDFLAGS)

$(analyzeVCF): $(OBJDIR) $(OBJS4)
	$(CXX) $(CXXFLAGS) -o $(analyzeVCF) $(OBJS4) $(LIBS) $(LDFLAGS)

# Pattern rule for object files
$(OBJDIR)/%.o: $(SRCDIR)/%.cpp | $(OBJDIR)
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

clean:
	rm -f $(OBJDIR)/*.o $(analyzeBAM) $(filterBAM) $(deamBAM) $(analyzeVCF)
	rmdir $(OBJDIR) 2>/dev/null || true
