/**
 * MuDirac - A muonic atom Dirac equation solver
 * by Simone Sturniolo (2019)
 * 
 * log.hpp
 * 
 * Functions for logging - header file
 * 
 * @author Simone Sturniolo
 * @version 0.1 04/02/2019
 */

#include <iostream>
#include <fstream>

using namespace std;

class Log
{
public:
    Log(string filename);
    ~Log();

private:
    std::ofstream logfile;
    std::streambuf *fbuf, *clbuf;
};