/*
 * File:   seed_header.h
 * Author: hargalaten
 *
 * Created on 12 dicembre 2012, 15.45
 */

#ifndef LTM_SWIPL_H
#define	LTM_SWIPL_H

#include <iostream>
#include <vector>
#include <string>
#include <algorithm>    // std::replace
#include <sstream>      // std::stringstream
#include <unordered_map>

#include <SWI-cpp.h>
#include <SWI-Prolog.h>


//LTM class
class swipl_interface {
public:

    swipl_interface();

    swipl_interface(int &, char **);

    ~swipl_interface();

    static std::vector<std::string> instance2vector(std::string);

    static std::vector<std::string> list2vector(std::string);
    
    bool consult(std::string);
    
    std::string query(std::string);
    
private:
    
    PlTerm string2term(std::string str, std::unordered_map<std::string, PlTerm> &);
    
    std::string prolog_path;
    
    PlEngine *swi_engine;
    
};



#endif	/* LTM_SWIPL_H */

