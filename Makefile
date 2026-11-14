CXX ?= g++
CPPFLAGS ?=
CXXFLAGS ?= -O3 -std=c++17 -Wall -Wextra -Wpedantic
LDFLAGS ?=
LDLIBS ?=

USE_ZLIB ?= 0
ifeq ($(USE_ZLIB),1)
CPPFLAGS += -DZLIB
LDLIBS += -lz
endif

ifeq ($(OS),Windows_NT)
TARGET := intergenic.exe
POWERSHELL ?= powershell
else
TARGET := intergenic
POWERSHELL ?= pwsh
endif
SOURCES := src/Intergenic.cpp src/Classes.cpp src/Utils.cpp
HEADERS := include/Classes.hpp include/Utils.hpp include/cxxopts.hpp

.PHONY: all clean test

all: $(TARGET)

$(TARGET): $(SOURCES) $(HEADERS)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -Iinclude $(SOURCES) -o $@ $(LDFLAGS) $(LDLIBS)

test: $(TARGET)
	$(POWERSHELL) -NoProfile -ExecutionPolicy Bypass -File tests/smoke.ps1 -Binary ./$(TARGET)
	$(POWERSHELL) -NoProfile -ExecutionPolicy Bypass -File tests/comprehensive.ps1 -Binary ./$(TARGET)

ifeq ($(OS),Windows_NT)
clean:
	powershell -NoProfile -Command "Remove-Item -LiteralPath '$(TARGET)' -Force -ErrorAction SilentlyContinue; exit 0"
else
clean:
	$(RM) $(TARGET)
endif
