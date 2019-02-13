# Test programs
test/test_main.o: test/test_main.cpp test/catch/catch.hpp
	g++ test/test_main.cpp -c -o test/test_main.o

test/test_utils: test/test_main.o test/test_utils.cpp obj/utils.o
	g++ test/test_main.o test/test_utils.cpp obj/utils.o -o test/test_utils

test/test_hydrogenic: test/test_hydrogenic.cpp obj/hydrogenic.o obj/utils.o src/constants.hpp
	g++ test/test_hydrogenic.cpp obj/hydrogenic.o obj/utils.o -o test/test_hydrogenic

test/test_integrate: test/test_integrate.cpp obj/integrate.o obj/utils.o 
	g++ test/test_integrate.cpp obj/integrate.o obj/utils.o -o test/test_integrate

# OBJ files
obj/hydrogenic.o: src/hydrogenic.cpp src/hydrogenic.hpp src/constants.hpp 
	g++ src/hydrogenic.cpp -c -o obj/hydrogenic.o

obj/utils.o: src/utils.cpp src/utils.hpp
	g++ src/utils.cpp -c -o obj/utils.o

obj/integrate.o: src/integrate.cpp src/integrate.hpp
	g++ src/integrate.cpp -c -o obj/integrate.o