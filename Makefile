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
CXXFLAGS	:= -Wall -g -O2 `root-config --cflags --libs ` -lMinuit -I$(INC) -L$(LIB) $(INCLUDES) $(LIBRARIES)

all: $(EXECUTABLE)

.PHONY: clean
clean:
	-$(RM) $(EXECUTABLE)
	-$(RM) $(OBJECTS)

run: all
	./$(EXECUTABLE)

bin/effmat.exe: src/effmat.cpp src/*.cc
	-$(CXX) $^ -o $@ $(CXXFLAGS)
bin/addweight_window.exe: src/addweight_window.cpp src/*.cc
	-$(CXX) $^ -o $@ $(CXXFLAGS)
bin/addweight_allsky.exe: src/addweight_allsky.cpp src/*.cc
	-$(CXX) $^ -o $@ $(CXXFLAGS)
bin/map.exe: src/map.cpp src/*.cc
	-$(CXX) $^ -o $@ $(CXXFLAGS)
bin/sim_par_select.exe: src/sim_par_select.cpp src/*.cc
	-$(CXX) $^ -o $@ $(CXXFLAGS)
bin/trackstat.exe: src/trackstat.cpp src/*.cc
	-$(CXX) $^ -o $@ $(CXXFLAGS)
bin/clonetree/par_select.exe: src/clonetree/par_select.cpp
	-$(CXX) $^ -o $@ $(CXXFLAGS)
bin/clonetree/window_select.exe: src/clonetree/window_select.cpp src/*.cc
	-$(CXX) $^ -o $@ $(CXXFLAGS)
bin/operator/par_select.exe: src/operator/par_select.cpp src/operator/*.cc
	-$(CXX) $^ -o $@ $(CXXFLAGS)
bin/rebase/par_select.exe: src/rebase/par_select.cpp src/rebase/*.cc
	-$(CXX) $^ -o $@ $(CXXFLAGS)
bin/rebase/window_select.exe: src/rebase/window_select.cpp src/rebase/*.cc src/*.cc
	-$(CXX) $^ -o $@ $(CXXFLAGS)
bin/main.exe: src/main.cpp
	-$(CXX) $^ -o $@ $(CXXFLAGS)