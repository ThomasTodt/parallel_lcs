# Makefile

# Compiler and flags
CC = gcc
CFLAGS = -O3 -fopenmp -Wall

# Targets
TARGETS = lcs lcs2

# Default rule
all: $(TARGETS)

# Compile lcs from lcs.c
lcs: lcs.c
	$(CC) $(CFLAGS) -o lcs lcs.c

# Compile lcs2 from lcs2.c
lcs2: lcs2.c
	$(CC) $(CFLAGS) -o lcs2 lcs2.c

# Clean rule to remove build artifacts
clean:
	rm -f $(TARGETS)
