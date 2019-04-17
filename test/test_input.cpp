#include "../src/input.hpp"
#include <iostream>

#include "catch/catch.hpp"

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
    }
}