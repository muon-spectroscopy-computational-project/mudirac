/**
 * MuDirac - A muonic atom Dirac equation solver
 * by Simone Sturniolo (2019)
 * 
 * input.cpp
 * 
 * Functions to parse input files
 * 
 * @author Simone Sturniolo
 * @version 0.1 04/02/2019
 */

#include "input.hpp"
#include <iostream>

/**
 * @brief  Instantiates an InputNode
 * @note   Instantiates an empty InputNode
 * @retval 
 */
template <typename T>
InputNode<T>::InputNode()
{
    this->single = false;
    this->default_value = vector<T>(0);
    this->clear();
}

/**
 * @brief  Instantiates an InputNode
 * @note   Instantiates an InputNode using a single default value
 * 
 * @param  default_value
 * @retval 
 */
template <typename T>
InputNode<T>::InputNode(T default_value)
{
    this->single = true;
    this->default_value = vector<T>(1, default_value);
    this->clear();
}

/**
 * @brief  Instantiates an InputNode
 * @note   Instantiates an InputNode using multiple default values.
 * The resulting InputNode will parse comma/space separated vectors.
 * 
 * @param  default_value
 * @retval 
 */
template <typename T>
InputNode<T>::InputNode(vector<T> default_value)
{
    this->single = false;
    this->default_value = vector<T>(default_value);
    this->clear();
}

/**
 * @brief  Copy constructor
 * @note   Copy constructor
 * 
 * @param  &t:  Object to copy
 * @retval 
 */
template <typename T>
InputNode<T>::InputNode(const InputNode<T> &t)
{
    single = t.single;
    default_value = t.default_value;
    value = t.value;
}

// Parsers, type by type
template <>
int InputNode<int>::parse(string s) const
{
    return stoi(s);
}

template <>
float InputNode<float>::parse(string s) const
{
    return stof(s);
}

template <>
double InputNode<double>::parse(string s) const
{
    return stod(s);
}

template <>
string InputNode<string>::parse(string s) const
{
    return s;
}

/**
 * @brief  Clears the contents of an InputNode
 * @note   Returns the contents of an InputNode to the default value.
 * 
 * @retval None
 */
template <typename T>
void InputNode<T>::clear()
{
    value = vector<T>(default_value);
}

/**
 * @brief  Returns the size of an InputNode
 * @note   Returns the number of values held by the InputNode.
 * 
 * @retval 
 */
template <typename T>
int InputNode<T>::getSize() const
{
    return single ? 1 : value.size();
}

/**
 * @brief  Returns a specific value from the InputNode
 * @note   Returns a specific value of index i from the InputNode.
 * 
 * @param  i:   Index of the requested value
 * @retval      Value
 */
template <typename T>
T InputNode<T>::getValue(int i) const
{
    // Avoiding pointless mistakes
    i = single ? 0 : i;
    return value[i];
}

/**
 * @brief  Returns all values from the InputNode
 * @note   Returns all values from the InputNode as a vector.
 * 
 * @retval Vector of values
 */
template <typename T>
vector<T> InputNode<T>::getValues() const
{
    return vector<T>(value);
}

/**
 * @brief  Set the InputNode value
 * @note   Set the InputNode value manually with a single value.
 * Will throw an error if the InputNode is vector-valued.
 * 
 * @param  v:       New value
 * @retval None
 */
template <typename T>
void InputNode<T>::setValue(T v)
{
    if (single)
    {
        value[0] = v;
    }
    else
    {
        throw "Can't set single value on vector-valued InputNode";
    }
}

/**
 * @brief  Set the InputNode values
 * @note   Set the InputNode values manually with a vector.
 * Will throw an error if the InputNode is single-valued.
 * 
 * @param  v:       New value
 * @retval None
 */
template <typename T>
void InputNode<T>::setValues(vector<T> v)
{
    if (single && v.size() > 1)
    {
        throw "Can't set vector of values for single-valued InputNode";
    }
    value = v;
}

/**
 * @brief  Parses a string into values
 * @note   Parses a string into one or multiple values for the InputNode.
 * If the InputNode accepts multiple values, this will split the string
 * by commas, spaces and/or tabs, merging separators, and parse the
 * individual fragments.
 * 
 * @param  s:       String to parse
 * @retval None
 */
template <typename T>
void InputNode<T>::parseValue(string s)
{
    if (single)
    {
        value[0] = parse(s);
    }
    else
    {
        // Split string values
        vector<string> svals = splitString(s, ", \t", true);
        value.clear();
        for (int i = 0; i < svals.size(); ++i)
        {
            value.push_back(parse(svals[i]));
        }
    }
}

// Acceptable types for InputNode
template class InputNode<int>;
template class InputNode<float>;
template class InputNode<double>;
template class InputNode<string>;

/**
 * @brief  Instantiates an InputFile
 * @note   Instantiates an InputFile.
 */
InputFile::InputFile()
{
}

/**
 * @brief  Copy data schema from an existing InputFile
 * @note   Copy the names and default values of all InputNodes in a
 * given InputFile, without copying the actual values
 * 
 * @param  schema: InputFile to use as a reference
 * @retval None
 */
void InputFile::copySchema(InputFile schema)
{
    // Copy string values
    string_values.clear();
    for (map<string, InputNode<string>>::iterator it = schema.string_values.begin();
         it != schema.string_values.end(); ++it)
    {
        InputNode<string> inode = it->second;
        inode.clear();
        string_values[it->first] = inode;
    }

    double_values.clear();
    for (map<string, InputNode<double>>::iterator it = schema.double_values.begin();
         it != schema.double_values.end(); ++it)
    {
        InputNode<double> inode = it->second;
        inode.clear();
        double_values[it->first] = inode;
    }

    int_values.clear();
    for (map<string, InputNode<int>>::iterator it = schema.int_values.begin();
         it != schema.int_values.end(); ++it)
    {
        InputNode<int> inode = it->second;
        inode.clear();
        int_values[it->first] = inode;
    }
}

void InputFile::parseFile(string path)
{
    string line;
    string key, value;
    vector<string> tokens;
    ifstream ifile(path);

    while (getline(ifile, line))
    {
        // Remove comments
        line = splitString(line, "#")[0];
        line = stripString(line);
        if (line.size() == 0)
            continue;

        tokens = splitString(line, ":");
        key = stripString(tokens[0]);
        value = stripString(tokens[1]);
        if (string_values.find(key) != string_values.end())
        {
            string_values[key].parseValue(value);
        }
        else if (double_values.find(key) != double_values.end())
        {
            double_values[key].parseValue(value);
        }
        else if (int_values.find(key) != int_values.end())
        {
            int_values[key].parseValue(value);
        }
    }
}

/**
 * @brief  Return a list of keys for string type values
 * @note   Return a list of keys for string type values.
 * 
 * @retval  Vector of string keys
 */
vector<string> InputFile::getStringKeys()
{
    vector<string> keys(0);

    for (map<string, InputNode<string>>::iterator it = string_values.begin();
         it != string_values.end(); ++it)
    {
        keys.push_back(it->first);
    }

    return keys;
}

/**
 * @brief  Return a list of keys for double type values
 * @note   Return a list of keys for double type values.
 * 
 * @retval  Vector of double keys
 */
vector<string> InputFile::getDoubleKeys()
{
    vector<string> keys(0);

    for (map<string, InputNode<double>>::iterator it = double_values.begin();
         it != double_values.end(); ++it)
    {
        keys.push_back(it->first);
    }

    return keys;
}

/**
 * @brief  Return a list of keys for integer type values
 * @note   Return a list of keys for integer type values.
 * 
 * @retval  Vector of integer keys
 */
vector<string> InputFile::getIntKeys()
{
    vector<string> keys(0);

    for (map<string, InputNode<int>>::iterator it = int_values.begin();
         it != int_values.end(); ++it)
    {
        keys.push_back(it->first);
    }

    return keys;
}

/**
 * @brief  Get a single string value
 * @note   Get a single string value from the given key
 * 
 * @param  k:   string value key
 * @param  i:   string value index (for vector values)
 * @retval      Value
 */
string InputFile::getStringValue(string k, int i)
{
    if (string_values.find(k) == string_values.end())
    {
        throw "Invalid key";
    }
    return string_values[k].getValue(i);
}

/**
 * @brief  Get a single double value
 * @note   Get a single double value from the given key
 * 
 * @param  k:   double value key
 * @param  i:   double value index (for vector values)
 * @retval      Value
 */
double InputFile::getDoubleValue(string k, int i)
{
    if (double_values.find(k) == double_values.end())
    {
        throw "Invalid key";
    }
    return double_values[k].getValue(i);
}

/**
 * @brief  Get a single integer value
 * @note   Get a single integer value from the given key
 * 
 * @param  k:   integer value key
 * @param  i:   integer value index (for vector values)
 * @retval      Value
 */
int InputFile::getIntValue(string k, int i)
{
    if (int_values.find(k) == int_values.end())
    {
        throw "Invalid key";
    }
    return int_values[k].getValue(i);
}

/**
 * @brief  Get a vector of string values
 * @note   Get a vector of string values for the given key.
 * 
 * @param  k:   string value key
 * @retval      Value
 */
vector<string> InputFile::getStringValues(string k)
{
    if (string_values.find(k) == string_values.end())
    {
        throw "Invalid key";
    }
    return string_values[k].getValues();
}

/**
 * @brief  Get a vector of double values
 * @note   Get a vector of double values for the given key.
 * 
 * @param  k:   double value key
 * @retval      Value
 */
vector<double> InputFile::getDoubleValues(string k)
{
    if (double_values.find(k) == double_values.end())
    {
        throw "Invalid key";
    }
    return double_values[k].getValues();
}

/**
 * @brief  Get a vector of integer values
 * @note   Get a vector of integer values for the given key.
 * 
 * @param  k:   integer value key
 * @retval      Value
 */
vector<int> InputFile::getIntValues(string k)
{
    if (int_values.find(k) == int_values.end())
    {
        throw "Invalid key";
    }
    return int_values[k].getValues();
}

/**
 * @brief  Define a string valued node
 * @note   Define an InputNode with a string value and the given name
 * 
 * @param  name: Name of the new InputNode
 * @param  node: Starting value of the InputNode
 * @retval None
 */
void InputFile::defineStringNode(string name, InputNode<string> node)
{
    string_values[name] = node;
}

/**
 * @brief  Define a double valued node
 * @note   Define an InputNode with a double value and the given name
 * 
 * @param  name: Name of the new InputNode
 * @param  node: Starting value of the InputNode
 * @retval None
 */
void InputFile::defineDoubleNode(string name, InputNode<double> node)
{
    double_values[name] = node;
}

/**
 * @brief  Define a integer valued node
 * @note   Define an InputNode with a integer value and the given name
 * 
 * @param  name: Name of the new InputNode
 * @param  node: Starting value of the InputNode
 * @retval None
 */
void InputFile::defineIntNode(string name, InputNode<int> node)
{
    int_values[name] = node;
}