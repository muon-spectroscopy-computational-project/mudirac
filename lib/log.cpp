/**
 * MuDirac - A muonic atom Dirac equation solver
 * by Simone Sturniolo (2019)
 * 
 * log.cpp
 * 
 * Functions for logging
 * 
 * @author Simone Sturniolo
 * @version 0.1 04/02/2019
 */

#include "log.hpp"

Log::Log(string filename)
{
    logfile.open(filename);
    clbuf = std::clog.rdbuf();
    fbuf = logfile.rdbuf();
    std::clog.rdbuf(fbuf);
}

Log::~Log()
{
    std::clog.rdbuf(clbuf);
    logfile.close();
    std::cout << "Log destroyed!\n";
}