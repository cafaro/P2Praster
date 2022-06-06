CC=clang++
CFLAGS=-Os -std=c++14 -I/opt/local/include
LDFLAGS=-lm -L/opt/local/lib -ligraph -lboost_program_options-mt
SOURCES=error.cpp statistic.cpp raster.cpp main.cpp
EXECUTABLE=raster

all:
	$(CC) $(CFLAGS) -o $(EXECUTABLE) $(SOURCES) $(LDFLAGS)

clean:
	rm -f *~ $(EXECUTABLE)
