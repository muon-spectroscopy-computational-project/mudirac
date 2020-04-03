/**
 * MuDirac - A muonic atom Dirac equation solver
 * by Simone Sturniolo (2019-2020)
 * 
 * input.hpp
 * 
 * Functions to parse input files - header file
 * 
 * @author Simone Sturniolo
 * @version 1.0 20/03/2020
 */

#include <string>
#include <vector>
#include <map>
#include <fstream>

#include "utils.hpp"

using namespace std;

#ifndef MUDIRAC_INPUT
#define MUDIRAC_INPUT

/**
 * @class InputNode
 * 
 * @brief A container and parser for input file values
 * @note  A class that contains and parses single or multiple
 * values from an input file
 * 
 */
template <typename T>
class InputNode
{
private:
    bool single;
    bool case_sensitive;
    vector<T> default_value;
    vector<T> value;

    T parse(string s) const;

public:
    InputNode();
    InputNode(const InputNode<T> &t);
    InputNode(T default_value, bool case_sensitive = true);
    InputNode(vector<T> default_value, bool case_sensitive = true);

    void clear();
    int getSize() const;

    void setValue(T v);
    void setValues(vector<T> v);
    void parseValue(string s);
    T getValue(int i = 0) const;
    vector<T> getValues() const;
};

/**
 * @class InputFile
 * 
 * @brief  A container class for entire input files
 * @note   A container class that parses entire input files. Uses
 * InputNode elements as a schema to validate the files it's parsing.
 * 
 */
class InputFile
{
private:
    map<string, InputNode<string>> string_values;
    map<string, InputNode<bool>> bool_values;
    map<string, InputNode<double>> double_values;
    map<string, InputNode<int>> int_values;

public:
    InputFile(void);

    void copySchema(InputFile schema);
    void parseFile(string path);

    vector<string> getStringKeys();
    vector<string> getBoolKeys();
    vector<string> getDoubleKeys();
    vector<string> getIntKeys();
    string getStringValue(string k, int i = 0);
    bool getBoolValue(string k, int i = 0);
    double getDoubleValue(string k, int i = 0);
    int getIntValue(string k, int i = 0);
    vector<string> getStringValues(string k);
    vector<bool> getBoolValues(string k);
    vector<double> getDoubleValues(string k);
    vector<int> getIntValues(string k);

    void defineStringNode(string name, InputNode<string> node = InputNode<string>(""));
    void defineBoolNode(string name, InputNode<bool> node = InputNode<bool>(false));
    void defineDoubleNode(string name, InputNode<double> node = InputNode<double>(0));
    void defineIntNode(string name, InputNode<int> node = InputNode<int>(0));
};

#endif