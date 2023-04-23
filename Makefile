CFLAGS = -Isrc/include

all: 
		g++ $(CFLAGS) -O3 -fomit-frame-pointer -fprefetch-loop-arrays -DNDEBUG src/textParser.cpp src/KMACS.cpp src/runFile.cpp src/suffixArray.c -o kmacs
