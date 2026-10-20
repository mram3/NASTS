#pragma once

#include <iostream>
#include <string>

class Tests
{
private:
    inline static int pass_ = 0;
    inline static int fail_ = 0;

public:

    //Expecting equality
    template <typename A, typename B>
    static void expect_equal
    (
        const A& actual,
        const B& real,
        const std::string& testName
    )
    {
        if(actual == real){
            std::cout << "\033[32m[ PASS ]\033[0m" << testName << "\n";
            pass_++;
        }

        else{
            std::cout << "\033[31m[ FAIL ]\033[0m" << testName 
                      << "(Expected : " << real << ", Got : " << actual << ")\n";
            fail_++;
        }
    }

    //Expecting true
    static void expect_true
    (
        bool condition,
        const std::string& testName
    )
    {
        if(condition){
            std::cout << "\033[32m[ PASS ]\033[0m" << testName << "\n";
            pass_++;
        }
        else{
            std::cout << "\033[31m[ FAIL ]\033[0m" << testName << "\n";
            fail_++;
        }
    }

    static int summary(){
        std::cout << "************************************************\n"
                  << "                 TEST SUMMARY                   \n";

        if(fail_ > 0){
            std::cout << "\033[31m[ FAILED ]\033[0m\n"
                      << "Tests Failed : " << fail_ << "\n"
                      << "Tests Passed : " << pass_ << "\n";
            return 1;
        }
        else{
            std::cout << "\033[32m[ PASSED ]\033[0m\n"
                      << "Tests Passed: " << pass_ << "\n";
            return 0;
        }
    }
};

#define EXPECT_EQ(actual, real, name) Tests::expect_equal(actual, real, name)
#define EXPECT_TRUE(condition, name) Tests::expect_true(condition, name)