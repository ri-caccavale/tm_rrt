/*
 * File:   seed_header.h
 * Author: hargalaten
 *
 * Created on 12 dicembre 2012, 15.45
 */

#ifndef ECLIPSECLP_INTERFACE_H
#define	ECLIPSECLP_INTERFACE_H


#include <iostream>
#include <string>
#include <sstream>
#include <cstring>
#include <map>
#include <math.h>
#include <fstream>
#include <vector>


#include "eclipseclass.h"



//Eclipseclp_interface class
class Eclipseclp_interface {
public:
    
    Eclipseclp_interface(std::string path_to_eclipseclp_directory = "/eclipseclp_interface");
    
    bool open(std::string);

    std::vector<std::string> query(std::string);
    
    bool close();
    
private:
    
    EC_word writeList(std::vector<std::string>);

    std::vector<std::string> readList(EC_word);

    void printList(EC_word);

    std::string functor2string(EC_word);
    
    std::string eclpse_path;
    
    std::string sourcefile_path;
    
};

std::vector<std::string> instance2vector(std::string schemaInstance);


#endif	/* ECLIPSECLP_INTERFACE_H */

