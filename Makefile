CXX			:= g++
CXXFLAGS	:= -Wall -g -O2 `root-config --cflags --libs ` -lMinuit

BIN		:= bin
SRC		:= src
INC		:= include
LIB		:= lib

LIBRARIES	:= -L ~/lib -lsofa_c
INCLUDES	:= -I ~/include

SOURCEDIRS	:= $(shell find $(SRC) -type d)

SOURCES		:= $(wildcard $(patsubst %,%/*.cpp, $(SOURCEDIRS)))
EXECUTABLE	:= $(SOURCES:.cpp=.exe)
EXECUTABLE	:= $(EXECUTABLE:.exe=.exe1)
EXECUTABLE	:= $(subst $(BIN)%,$(SRC)%,$(EXECUTABLE))

all: $(EXECUTABLE)

.PHONY: clean
clean:
	-$(RM) $(EXECUTABLE)

run: all
	./$(EXECUTABLE)

$(EXECUTABLE): $(BIN)/%.exe: $(SRC)/%.cpp $(SRC)/$(dir %)/*.cc
	$(CXX) $^ -o $@ $(CXXFLAGS) -I$(INC) -I$(INC)/$(dir %) $(CLIBS) $(INCLUDES) $(LIBRARIES)