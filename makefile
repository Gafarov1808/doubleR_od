TARGET = determ
SOURCES = main.cpp
OBJECTS = $(SOURCES:.cpp=.o)

PYTHON = python3.11
PYTHON_VERSION = 3.11
MODULE = doubleR_module`python3.11-config --extension-suffix`
BINDINGS = bindings.cpp
BINDINGS_OBJECTS = bindings.o
MODULE_OBJECTS = $(BINDINGS_OBJECTS)

CXX = g++
CXXFLAGS = -Wall -std=c++17 -O3 -fPIC

PYTHON_INCLUDES = $(shell $(PYTHON)-config --includes)
PYBIND11_LDFLAGS = $(shell $(PYTHON)-config --ldflags)
PYBIND11_INCLUDES = $(shell $(PYTHON) -c "import pybind11; print(pybind11.get_include())" 2>/dev/null || echo "")

INCLUDES = $(PYTHON_INCLUDES) -I$(PYBIND11_INCLUDES)
PYTHON_LDFLAGS = $(shell $(PYTHON)-config --ldflags)

all: $(TARGET) module

$(TARGET): $(OBJECTS)
	$(CXX) $(OBJECTS) -o $(TARGET)

module: $(MODULE)

$(MODULE): $(BINDINGS_OBJECTS)
	$(CXX) -shared bindings.o -o $(MODULE) $(PYTHON_LDFLAGS) 

main.o: main.cpp doubleR.h
	$(CXX) $(CXXFLAGS) -c main.cpp -o main.o

bindings.o: bindings.cpp doubleR.h
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c bindings.cpp -o bindings.o

clean:
	rm -f main.o bindings.o $(TARGET) $(MODULE)

distclean: clean
	rm -f *.o *.so *.pyd *.dll

install: module
	pip install .

test: module
	python3 -c "import doubleR_module; print('Модуль успешно загружен!')"

.PHONY: all clean distclean module install test