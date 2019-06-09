CC	= g++
CPPFLAGS= -Wno-deprecated -O3 -c -D_FILE_OFFSET_BITS=64 
LDFLAGS	= -O3 
SOURCES	= texttobin.cpp
OBJECTS	= $(SOURCES:.cpp=.o)
EXECUTABLE=texttobin

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE) : $(OBJECTS)
	$(CC) $(LDFLAGS) $(OBJECTS) -o $@

.cpp.o : 
	$(CC) $(CPPFLAGS) $< -o $@

clean:
	rm -f *.o
