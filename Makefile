# Test programs
test/test_main.o: test/test_main.cpp test/catch/catch.hpp
	g++ test/test_main.cpp -c -o test/test_main.o

test/test_utils: test/test_main.o test/test_utils.cpp obj/utils.o
	g++ test/test_main.o test/test_utils.cpp obj/utils.o -o test/test_utils

test/test_hydrogenic: test/test_hydrogenic.cpp obj/hydrogenic.o obj/utils.o src/constants.hpp
	g++ test/test_hydrogenic.cpp obj/hydrogenic.o obj/utils.o -o test/test_hydrogenic

test/test_integrate: test/test_integrate.cpp obj/integrate.o obj/utils.o obj/boundary.o
	g++ test/test_integrate.cpp obj/integrate.o obj/utils.o obj/boundary.o -o test/test_integrate

test/test_atom: test/test_atom.cpp obj/atom.o obj/utils.o obj/integrate.o obj/hydrogenic.o src/constants.hpp
	g++ test/test_atom.cpp obj/atom.o obj/utils.o obj/integrate.o obj/hydrogenic.o -o test/test_atom

# OBJ files
obj/hydrogenic.o: src/hydrogenic.cpp src/hydrogenic.hpp src/constants.hpp 
	g++ src/hydrogenic.cpp -c -o obj/hydrogenic.o

obj/utils.o: src/utils.cpp src/utils.hpp
	g++ src/utils.cpp -c -o obj/utils.o

obj/integrate.o: src/integrate.cpp src/integrate.hpp
	g++ src/integrate.cpp -c -o obj/integrate.o

obj/boundary.o: src/boundary.cpp src/boundary.hpp
	g++ src/boundary.cpp -c -o obj/boundary.o

obj/atom.o: src/atom.cpp src/atom.hpp
	g++ src/atom.cpp -c -o obj/atom.o