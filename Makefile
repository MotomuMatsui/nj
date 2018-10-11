CXX := g++

VPATH := src
INC   := -Ilib
LIB   := -Llib -lm

CXXFLAGS := -O3
CXXFLAGS += -std=c++11
CXXFLAGS += -march=native
CXXFLAGS += -fno-exceptions
CXXFLAGS += -Wall

OBJECTS := format.o
OBJECTS += messages.o
OBJECTS += nj.o
OBJECTS += ep.o
OBJECTS += ep_function.o
OBJECTS += main.o

.PHONY: all
all: nj clean

nj: $(OBJECTS)
	$(CXX) $(CXXFLAGS) -o $@ $^

$(OBJECTS): %.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $<

.PHONY: clean
clean:
	$(RM) $(OBJECTS)
