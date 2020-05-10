BIN		:= bin
SRC		:= src
INC		:= include
LIB		:= lib

LIBRARIES	:= -L ~/lib -lsofa_c
INCLUDES	:= -I ~/include

SOURCEDIRS	:= $(shell find $(SRC) -type d)

SOURCES		:= $(wildcard $(patsubst %,%/*.cpp, $(SOURCEDIRS)))
EXECUTABLE	:= $(SOURCES:.cpp=.exe)
EXECUTABLE	:= $(subst $(SRC), $(BIN), $(EXECUTABLE))

CXX			:= g++
CXXFLAGS	:= -Wall -g -O2 `root-config --cflags --libs ` -lMinuit -I$(INC) $(CLIBS) $(INCLUDES) $(LIBRARIES)

all: $(EXECUTABLE)

.PHONY: clean
clean:
	-$(RM) $(EXECUTABLE)

run: all
	./$(EXECUTABLE)

$(EXECUTABLE): $(BIN)/%.exe: $(SRC)/%.cpp
	-$(CXX) $^  $(SRC)/$(*D)/*.cc -o $@ $(CXXFLAGS)