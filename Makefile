# Test programs
test/test_hydrogenic: test/test_hydrogenic.cpp obj/hydrogenic.o obj/utils.o src/constants.hpp
	g++ test/test_hydrogenic.cpp obj/hydrogenic.o obj/utils.o -o test/test_hydrogenic

test/test_utils: test/test_utils.cpp obj/utils.o
	g++ test/test_utils.cpp obj/utils.o -o test/test_utils

# OBJ files
obj/hydrogenic.o: src/hydrogenic.cpp src/hydrogenic.hpp src/constants.hpp 
	g++ src/hydrogenic.cpp -c -o obj/hydrogenic.o

obj/utils.o: src/utils.cpp src/utils.hpp
	g++ src/utils.cpp -c -o obj/utils.o