filename = main

ifeq ($(OS), Windows_NT)
	EXE = $(filename).exe
else
	EXE = ./$(filename)
endif

CC = g++
CFLAGS = -Wshadow -Wno-unused-result -Wno-sign-compare -Wno-char-subscripts -std=c++20 -Wall -Wextra -g -DPIZZA

main:
	$(CC) $(CFLAGS) -o $(EXE) $(filename).cpp && $(EXE) $(filter-out $@,$(MAKECMDGOALS))
%:
	@:


in:
	$(CC) $(CFLAGS) -o $(EXE) $(filename).cpp && $(EXE) < input.txt
