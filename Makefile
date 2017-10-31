F_COMP := gfortran
F_FILES := $(wildcard src/*.f90)
OBJ_FILES := $(addprefix obj/,$(notdir $(F_FILES:.f90=.o)))



all: main.out

main.out: $(OBJ_FILES)
	$(F_COMP) -o $@ $^

obj/%.o: src/%.f90
	$(F_COMP) -c -o $@ $<

clean:
	rm -rf $(OBJ_FILES)	

run:
	./main.out