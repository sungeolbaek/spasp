CC	= g++
CPPFLAGS= -O3 -c -D_FILE_OFFSET_BITS=64
LDFLAGS	= -O3 
SOURCES	= spasp_label.cpp
OBJECTS	= $(SOURCES:.cpp=.o)
EXECUTABLE=spasp_label

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE) : $(OBJECTS)
	$(CC) $(LDFLAGS) $(OBJECTS) -o $@

.cpp.o : 
	$(CC) $(CPPFLAGS) $< -o $@

clean:
	rm -f *.o
