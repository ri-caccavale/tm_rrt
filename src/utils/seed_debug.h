#include <iostream>

//text colors using ansi syntax... for linux terminal!
namespace ansi{
/*
         foreground background
black        30         40
red          31         41
green        32         42
yellow       33         43
blue         34         44
magenta      35         45
cyan         36         46
white        37         47
*/

//default style
const std::string end = "\033[0m";
//foreground colors
const std::string red     = "\033[31m";
const std::string green   = "\033[32m";
const std::string yellow  = "\033[33m";
const std::string blue    = "\033[34m";
const std::string magenta = "\033[35m";
const std::string cyan    = "\033[36m";
const std::string white   = "\033[37m";

//background colors
const std::string back_red       = "\033[41m";
const std::string back_green     = "\033[42m";
const std::string back_yellow    = "\033[43m";
const std::string back_blue      = "\033[44m";
const std::string back_magenta   = "\033[45m";
const std::string back_cyan      = "\033[46m";
const std::string back_white     = "\033[47m";


}

//convert string to number considering the std::locale::classic()
//  NOTE: the decimal are always separated by "." independently from the current system
//        this should be more reliable than atof and stod
inline double ston(std::string str){
    float number = 0.0f;
    std::istringstream istr(str);
    istr.imbue(std::locale::classic());
    istr >> number;
    return number;
}