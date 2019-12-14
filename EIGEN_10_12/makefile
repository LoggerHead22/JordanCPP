CC := g++
CFLAGS := -std=c++14 -Wall -pg  -ffast-math -march=native  -mfpmath=sse  

OBJS := gauss.o

a.out: main.o $(OBJS)
	$(CC) $(CFLAGS)-o $@ $^

clean:
	del *.o *.exe *.out

%.o: %.cpp *.h
	$(CC) $(CFLAGS) -c $<
