#include "mudirac.hpp"

int main(int argc, char *argv[])
{
    string seed = "mudirac";
    MuDiracInputFile config;

    if (argc > 1)
    {
        config.parseFile(argv[1]);
        seed = splitString(argv[1], ".")[0];
    }

    // Set up logging
    AixLog::Severity log_verbosity;
    switch (config.getIntValue("verbosity"))
    {
    case 1:
        log_verbosity = AixLog::Severity::info;
        break;
    case 2:
        log_verbosity = AixLog::Severity::debug;
        break;
    case 3:
        log_verbosity = AixLog::Severity::trace;
        break;
    default:
        log_verbosity = AixLog::Severity::info;
        break;
    }
    AixLog::Log::init({make_shared<AixLog::SinkFile>(log_verbosity, AixLog::Type::normal, seed + ".log"),
                       make_shared<AixLog::SinkFile>(AixLog::Severity::warning, AixLog::Type::special, seed + ".err")});

    LOG(INFO) << "MuDirac, a muonic atomic solver\n";
    LOG(INFO) << "by Simone Sturniolo\n";
    LOG(INFO) << "Released under the MIT License (2019)\n";
    LOG(INFO) << " \n";

    DiracAtom da = config.makeAtom();

    vector<string> lines = config.getStringValues("xrd_lines");

    for (int i = 0; i < lines.size(); ++i)
    {
        vector<string> states = splitString(lines[i], "-");
        cout << states[0] << '-' << states[1] << '\n';
    }
}