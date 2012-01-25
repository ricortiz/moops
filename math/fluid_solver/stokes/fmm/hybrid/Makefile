CC = nvcc
;CFLAGS = -O3 -lm -Xcompiler -fopenmp
CFLAGS = -O3 -ccbin icc -Xcompiler -openmp
OBJS = tree.o expansions.o
multipole : $(OBJS) driver.h driver.c
	$(CC) $(CFLAGS) -L. -lkernel driver.c -o multipole $(OBJS)
%.o: %.c driver.h
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	rm -rf *o multipole
