CXX      ?= g++
CXXFLAGS ?= -std=c++17 -O2 -Wall -Wextra -Iinclude

BUILD := build
BINS  := $(BUILD)/EventPredictor $(BUILD)/PDKMCGenerator $(BUILD)/LunarPDKGenerator

.PHONY: all clean

all: $(BINS)

$(BUILD)/%: src/%.cpp $(wildcard include/*.h) | $(BUILD)
	$(CXX) $(CXXFLAGS) -o $@ $<

$(BUILD):
	mkdir -p $(BUILD)

clean:
	rm -f $(BINS)
