CC=g++
CFLAGS=-c -Wall -I/usr/include/boost -I./ -O3 -D_LARGEFILE64_SOURCE
SOURCES=Shingles.cpp AdjacencyMatrix.cpp vnmextract.cpp Trie.cpp
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=vnmextract

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(OBJECTS) -o $@

vnmextract.o: AdjacencyMatrix.h Trie.h Shingles.h vnmextract.cpp conf.h

AdjacencyMatrix.o: AdjacencyMatrix.cpp AdjacencyMatrix.h conf.h Utils.h

Trie.o: Trie.cpp Trie.h

Shingles.o: Shingles.cpp Shingles.h Utils.h conf.h 

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@

clean:
	rm -fr *.o vnmextract *~
