CXX ?= c++
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
else
TARGET := intergenic
endif
SOURCES := src/Intergenic.cpp src/Classes.cpp src/Motif.cpp src/Utils.cpp
HEADERS := include/Classes.hpp include/Motif.hpp include/Utils.hpp include/cxxopts.hpp

.PHONY: all clean

all: $(TARGET)

$(TARGET): $(SOURCES) $(HEADERS)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -Iinclude $(SOURCES) -o $@ $(LDFLAGS) $(LDLIBS)

ifeq ($(OS),Windows_NT)
clean:
	powershell -NoProfile -Command "Remove-Item -LiteralPath '$(TARGET)' -Force -ErrorAction SilentlyContinue; exit 0"
else
clean:
	$(RM) $(TARGET)
endif
