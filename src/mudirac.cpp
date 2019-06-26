#include "mudirac.hpp"

int main(int argc, char *argv[])
{
    MuDiracInputFile config;

    if (argc > 1) {
        config.parseFile(argv[1]);
    }
}