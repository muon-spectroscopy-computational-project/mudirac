#include "mudirac.hpp"

int main(int argc, char *argv[])
{
    MuDiracInputFile config;

    if (argc > 1)
    {
        config.parseFile(argv[1]);
    }

    // Prepare the DiracAtom
    vector<string> lines = config.getStringValues("xrd_lines");

    for (int i = 0; i < lines.size(); ++i) {
        cout << lines[i] << '\n';
    }
}