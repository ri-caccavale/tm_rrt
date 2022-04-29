/*
 * File:   main.cpp
 * Author: hargalaten
 *
 * Created on 4 dicembre 2012, 13.59
 */

#include "swipl_interface.h"

swipl_interface::swipl_interface(){
    
    swi_engine = new PlEngine("default");
    
}

swipl_interface::swipl_interface(int &argc, char **argv){
    
    swi_engine = new PlEngine(argv[0]);
    
}

swipl_interface::~swipl_interface(){
    
    delete swi_engine;
}

/**
 * Convert a prolog-like functor into a vector of arguments
 * 
 * @param functor: string representing a functor, e.g. "fun(x,y)"
 * @return vector of strings representing the arguments of the functor
 * 
 * This function is implemented as a pushdown automata and returns a
 * vector similar to the "univ" (=..) Prolog predicate:
 * EG.
 *      vec[0] -> "name of functor"
 *      vec[1] -> "first argument"
 *      vec[2] -> "second argument"
 *      etc.
 *
 * Nested functors are still returned as elements of the vector:
 * EG.
 *      functor2vector("fun1(fun2(x,y),z)")[1] == "fun2(x,y)"
 *
 * prolog-lists [] are also considered
 */
std::vector<std::string> swipl_interface::functor2vector(std::string functor){
    //assume by default this is an atom and not a string
    bool isAtom=true, isString=false;
    char c;
    std::string app;
    std::vector<std::string> result;
    std::stringstream ss(functor);
    int count=0;
    ss >> std::noskipws;
    //read the first character of the string
    ss>>c;
    //while inside the string
    while(!ss.eof())
    {
        // STRING CHECKING
        //if the character is a quotation mark
        if(c=='"' && !isString){
            //a string is detected
            isString=true;
            //store the mark to the argument
            app=app+c;
        }
        //if the character is a quotation mark but I'm into a string
        else if(c=='"' && isString){
            //the string is ended
            isString=false;
            //store the mark
            app=app+c;
        }
        //while reading a string
        else if(isString){
            //store the character without furhter checks
            app=app+c;
        }
        // ATOM CHECKING
        //if the character is an open bracket (the first one found)
        else if(c=='(' && isAtom){
            //this is not an atom anymore
            isAtom=false;
            //push this element as the name of the functor
            result.push_back(app);
            //clear the current string
            app="";
        }
        // NESTED FUNCTOR/LIST CHECKING
        //if the character is an open backet
        else if(c=='(' || c=='['){
            //this is either an inner functor or list so increase the pushdown stack
            count++;
            //store the character
            app=app+c;
        }
        //if the character is a closed bracket and the stack is not empty
        else if( ( c==')' || c==']' ) && count!=0){
            //decrease the stack
            count--;
            //store the character
            app=app+c;
        }
        // ARGUMENT CHECKING
        //if the character is a comma and the stack is not empty
        else if(c!=',' || count!=0)
            //store the character into the current argument
            app=app+c;
        //otherwise (the character is a comma with empty stack)
        else {
            //the argument is terminated, push it into the vector
            result.push_back(app);
            //clear the current string (in so skipping the comma)
            app="";
        }
        //read the next character
        ss>>c;
    }
    //if this functor is an atom (no arguments)
    if(isAtom) {
        //push the name into the vector
        result.push_back(app);
    }
    //otherwise, push the last argument and remove the final backet
    else{
        app.erase(app.size()-1);
        result.push_back(app);
    }
    //return the vector
    return result;
}

/**
 * Consult a Prolog file (database)
 * 
 * @param path: absolute path to the Prolog file to be consulted
 * @return true/false for success/fail
 */   
bool swipl_interface::consult(std::string path){
    
    prolog_path = path;

    //put the path as the first argument of the query
    PlTermv args(1);
    args[0] = prolog_path.c_str();
    //try to consult the database
    try{
        PlQuery pq("consult", args);

        if (pq.next_solution()) {
            std::cout << "SWIPL consult done"<<std::endl;
            return true;
        } else {
            std::cout << "SWIPL consult failed"<<std::endl;
            return false;
        }
    }
    catch(const PlException& e){
        std::cout << "SWIPL: cosult error!"<<std::endl;
        std::cout << (char*) e << std::endl;
        return false;
    }
}

/**
 * Convert a prolog-like list into a vector of elemets
 * 
 * @param list: string representing a list, e.g. "[x,y]"
 * @return vector of strings representing the elements of the list
 * 
 * The list-reading problem is reduced to a functor-reading problem
 * so functor2vector is used to preserve also possible inner functors
 * or lists, elements are returned as follows:
 * EG.
 *      vec[0] -> "first element"
 *      vec[2] -> "second element"
 *      etc.
 * The "name" of the list is not returned
 */
std::vector<std::string> swipl_interface::list2vector(std::string list){
    std::vector<std::string> lv;

    //rapidly check the syntax
    if(list.front() != '[' || list.back() != ']' )
        //return empty vector if the string is not a list
        return lv;

    //transform the list into a functor by replacing the brackets 
    std::replace( list.begin(), list.begin()+1, '[', '(' );
    std::replace( list.begin()+list.size()-1, list.end(), ']', ')' );
    
    //read the functor
    lv = functor2vector(list);
    //erase the first element (empty functor name)
    lv.erase(lv.begin());

    //for(auto i=0; i<lv.size(); i++)
    //    std::cout<<lv[i]<<std::endl;
    
    //return the vector
    return lv;
}

/**
 * Convert a prolog-like term string into a SWI term, non-ground terms
 * can be also converted
 * 
 * @param str: string representing the compound, e.g. "arg(P,f(a,b,c),a)"
 * @param swi_vars: map of variables already assocaited to open terms
 * @return SWI term
 */
PlTerm swipl_interface::string2term(std::string str, std::unordered_map<std::string, PlTerm> &swi_vars){
    std::vector<std::string> strv = functor2vector(str);
    
//    if(strv.size()<=0)
//        return PlCompound(); //not allowed!
    
    std::string fun_name = strv[0];
    
    PlTermv fun_args( strv.size()-1 );
    
    try {
        
        //if it is an atom
        if(strv.size()<=1){
            
            if(fun_name.at(0) != '_' && !isupper(fun_name.at(0)) ){
                //std::cout<<"\t\t term: "<<fun_name<<std::endl;
                return PlCompound(fun_name.c_str()); //PlTerm(fun_name.c_str());
            }
            else if(fun_name.at(0) != '_'){
                
                auto it = swi_vars.find(fun_name);
                
                if (it != swi_vars.end()){
                    //std::cout<<"\t\t old var: "<<fun_name<<std::endl;
                    return swi_vars[fun_name];
                }
                else {
                    //std::cout<<"\t\t new var: "<<fun_name<<std::endl;
                    swi_vars[fun_name] = PlTerm(); //PlCompound(fun_name.c_str());
                    return swi_vars[fun_name];
                }   
            }
        }
        
        for(auto i=1; i<strv.size(); i++){
            fun_args[i-1] = string2term(strv[i],swi_vars);
        }
        
    } catch (PlException& e) {
        std::cout << "SWIPL: Query error!"<<std::endl;
        std::cout << (char*) e << std::endl;
    }
    
    //std::cout<<"\t\t comp: "<<fun_name<<std::endl;
    return PlCompound(fun_name.c_str(),fun_args);
    
}


/**
 * Execute a query to the consulted database, it returns the first solution
 * as the instantiated query: having all variables unified to ground terms
 * 
 * @param str: string representing the query, e.g. "arg(P,f(a,b,c),a)"
 * @return string representing the instantiated query, e.g. "arg(1,f(a,b,c),a)"
 */
std::string swipl_interface::query(std::string request){
    //map associating variable names to open terms
    std::unordered_map<std::string, PlTerm> swi_vars;// = new std::unordered_map<std::string, term_t>();
    
    std::string from_prolog = "";
    
    std::vector<std::string> rv = functor2vector(request);
    
    std::stringstream ss;
    
    if( rv.size()<=1 )
        return "";
    
    //PL_thread_attach_engine(NULL);
    
    try {
        
        PlTermv inputs( rv.size()-1 );
        for(auto i=1; i<rv.size(); i++){
            //if(rv[i].at(0) != '_' && !isupper(rv[i].at(0)) )
                //inputs[i-1] = PlCompound(rv[i].c_str());
            //std::cout<<"\tstring2term: "<<rv[i]<<std::endl;
            inputs[i-1] = string2term(rv[i],swi_vars);
        }
        
        //send the request
        PlQuery pq(rv[0].c_str(), inputs);

        if( pq.next_solution() ) {
        
            ss<<rv[0]<<"(";
            
            bool not_first = false;
            
            for(auto i=0; i<rv.size()-1; i++){
                from_prolog = (char*) inputs[i];
                
                if(not_first)
                    ss<<",";
                else
                    not_first = true;
                
                ss<<from_prolog;
            }
            ss<<")";
            
        }
        
    } catch (PlException& e) {
        std::cout << "SWIPL: Query error!"<<std::endl;
        std::cout << (char*) e << std::endl;
    }
    
    //free SWI
    //PL_thread_destroy_engine();

    //return from_prolog;
    
    std::cout<<"RESULT: "<<ss.str()<<std::endl;

    return ss.str();
}

