TARGET = determ
SOURCES = main.cpp init_det.cpp
OBJECTS = $(SOURCES:.cpp=.o)

PYTHON = python3.12
PYTHON_VERSION = 3.12
MODULE = ODdoubleR$(shell python3.12-config --extension-suffix)
BINDINGS = bindings.cpp
BINDINGS_OBJECTS = bindings.o
MODULE_OBJECTS = $(BINDINGS_OBJECTS)

CXX = g++
CXXFLAGS = -Wall -std=c++17 -O3 -fPIC -I/usr/include/python3.12 $(shell python3 -m pybind11 --includes)

PYTHON_INCLUDES = $(shell $(PYTHON)-config --includes)
PYBIND11_LDFLAGS = $(shell $(PYTHON)-config --ldflags)
PYBIND11_INCLUDES = $(shell $(PYTHON) -c "import pybind11; print(pybind11.get_include())" 2>/dev/null || echo "")

INCLUDES = $(PYTHON_INCLUDES) -I$(PYBIND11_INCLUDES)
PYTHON_LDFLAGS = $(shell $(PYTHON)-config --ldflags)

all: $(TARGET) module

$(TARGET): $(OBJECTS)
	$(CXX) $(OBJECTS) -o $(TARGET)

module: $(MODULE)

$(MODULE): $(BINDINGS_OBJECTS) $(OBJECTS)
	$(CXX) -shared $^ -o $@ $(PYTHON_LDFLAGS)

main.o: main.cpp doubleR.h
	$(CXX) $(CXXFLAGS) -c main.cpp -o main.o

bindings.o: bindings.cpp doubleR.h
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c bindings.cpp -o bindings.o

init_det.o: init_det.cpp doubleR.h lambert.h
	$(CXX) $(CXXFLAGS) -c init_det.cpp -o init_det.o

clean:
	rm -f main.o bindings.o init_det.o $(TARGET) $(MODULE)

distclean: clean
	rm -f *.o *.so *.pyd *.dll

install: module
	pip install .

test: module
	python3 -c "import ODdoubleR; print('Модуль успешно загружен!')"

.PHONY: all clean distclean module install test