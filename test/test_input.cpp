#include "../lib/input.hpp"
#include <stdio.h>
#include <iostream>
#include <fstream>
#include "datapath.h"

#include "../vendor/catch/catch.hpp"

using namespace std;

TEST_CASE("Single input nodes", "[InputNode]")
{
    // Parse properly all types
    {
        InputNode<int> inode = InputNode<int>(0);
        inode.parseValue("3");
        REQUIRE(inode.getValue() == 3);
    }
    {
        InputNode<float> inode = InputNode<float>(0);
        inode.parseValue("1.2");
        REQUIRE(inode.getValue() == 1.2f);
    }
    {
        InputNode<double> inode = InputNode<double>(0);
        inode.parseValue("1.2");
        REQUIRE(inode.getValue() == 1.2);
    }
    {
        InputNode<string> inode = InputNode<string>("");
        inode.parseValue("test");
        REQUIRE(inode.getValue() == "test");
    }
    // Try arrays
    {
        InputNode<int> inode = InputNode<int>(vector<int>{});
        inode.parseValue("3,4, 5, 2");
        CHECK(inode.getValues() == vector<int>{3, 4, 5, 2});
    }
    // Set values manually
    {
        InputNode<int> inode = InputNode<int>(0);
        inode.setValue(3);
        REQUIRE(inode.getValue() == 3);
        REQUIRE_THROWS(inode.setValues(vector<int>{3, 4}));
    }
    {
        InputNode<int> inode = InputNode<int>(vector<int>{0, 1});
        inode.setValues(vector<int>{3, 5});
        CHECK(inode.getValues() == vector<int>{3, 5});
        REQUIRE_THROWS(inode.setValue(2));
    }
    {
        InputNode<bool> inode = InputNode<bool>(false);
        inode.parseValue("TRUE");
        REQUIRE(inode.getValue());
        inode.parseValue("FALSE");
        REQUIRE(!inode.getValue());
    }
}

TEST_CASE("Input file", "[InputFile]")
{
    {
        InputFile ifile = InputFile();
        ifile.defineStringNode("street", InputNode<string>("Downing St."));
        ifile.defineIntNode("number", InputNode<int>(10));

        CHECK(ifile.getStringKeys() == vector<string>{"street"});
        CHECK(ifile.getIntKeys() == vector<string>{"number"});
        REQUIRE(ifile.getIntValue("number") == 10);
        REQUIRE_THROWS(ifile.getDoubleValue("height"));

        // Try schema
        InputFile ifile_deriv = InputFile();
        ifile_deriv.copySchema(ifile);
        REQUIRE(ifile_deriv.getIntValue("number") == 10);
    }

    {
        InputFile ifile = InputFile();
        ifile.defineIntNode("integer", InputNode<int>(0));
        ifile.defineIntNode("integers", InputNode<int>(vector<int>{0}));
        ifile.defineBoolNode("bool", InputNode<bool>(true));
        ifile.defineBoolNode("bools", InputNode<bool>(vector<bool>{false}, false));
        ifile.defineDoubleNode("double", InputNode<double>(0.0));
        ifile.defineDoubleNode("doubles", InputNode<double>(vector<double>{0.0}));
        ifile.defineStringNode("string", InputNode<string>(""));
        ifile.defineStringNode("strings", InputNode<string>(vector<string>{""}));
        ifile.parseFile(string(CURRENT_DATAPATH) + "/data/inputtest.in");

        REQUIRE(ifile.getStringValue("string") == "mostly harmless");
        REQUIRE(ifile.getBoolValue("bool") == true);
        REQUIRE(ifile.getIntValue("integer") == 42);
        REQUIRE(ifile.getDoubleValue("double") == 3.14159);

        CHECK(ifile.getStringValues("strings") == vector<string>{"life", "the", "universe", "and", "everything"});
        CHECK(ifile.getBoolValues("bools") == vector<bool>{true, true, false});
        CHECK(ifile.getIntValues("integers") == vector<int>{6, 7, 42});
        CHECK(ifile.getDoubleValues("doubles") == vector<double>{1.414, 2.718, 3.142});
    }
}