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
 *
 * @param schemaInstance
 * @return vettore di stringhe rappresentanti i parametri dello schema
 * in versione prolog-like:
 * EG.
 *      vec[0]="nome schema"
 *      vec[1]="primo parametro"
 *      vec[2]="secondo parametro"
 *      etc.
 *
 * eventuali parametri che siano essi stessi funtori vengono
 * restituiti ugualmente come elemento del vettore
 * EG.
 *      vec[i]="fun1(fun2(x,y),z)"
 *
 * NOTE: This version also consider the [ ] as a list!
 */
std::vector<std::string> swipl_interface::instance2vector(std::string schemaInstance){
    bool isAtom=true, isString=false;
    char c;
    std::string app;
    std::vector<std::string> result;
    std::stringstream ss(schemaInstance);
    int count=0;
    ss >> std::noskipws;
    //leggi il primo carattere della stringa
    ss>>c;
    //mentre non sei a fine stringa
    while(!ss.eof())
    {
        //se il carattere è un doppio apice e non sono in una stringa
        if(c=='"' && !isString){
            //allora sono in una stringa
            isString=true;
            //aggiungo l'apice
            app=app+c;
        }
        //se il carattere è un doppio apice e sono in una stringa
        else if(c=='"' && isString){
            //la stringa è finita
            isString=false;
            //aggiungo l'apice
            app=app+c;
            //aggiungila come elemento del funtore
            //result.push_back(app);
        }
        //mentre sono in una stringa
        else if(isString){
            //aggiungi il carattere senza controllarlo
            app=app+c;
        }
        //se sono un atomo ed il carattere letto è una parentesi aperta
        else if(c=='(' && isAtom){
            //non sono più un atomo
            isAtom=false;
            //inserisco il nome come primo elemento del vettore
            result.push_back(app);
            //pulisco la stringa d'appoggio
            app="";
            //salto la parentesi
//            ss>>c;
        }
        else if(c=='(' || c=='['){
            count++;
            app=app+c;
        }
        else if( ( c==')' || c==']' ) && count!=0){
            count--;
            app=app+c;
        }
        //se il carattere letto non è una virgola
        else if(c!=',' || count!=0)
            //aggiungilo alla stringa d'appoggio
            app=app+c;
        //altrimenti (ie. il carattere è una virgola)
        else {
            //inserisci la stringa d'appoggio nel vettore risultato
            result.push_back(app);
            //pulisci la stringa d'appoggio
            app="";
            //ho saltato la virgola
        }
        //leggi il successivo carattere
        ss>>c;
    }
    //se lo schema non ha parametri aggiungi il solo nome (vec[0])
    if(isAtom) {
        //check the \ character and split by it (added 01/12/2020 in seed 4.0)
        if( app.find('\\') != std::string::npos ){
            std::stringstream ss2(app);
            std::string substr;
            //std::cout<<"INSTANCE TO VECTOR: "<<schemaInstance<<std::endl;
            while(std::getline(ss2, substr, '\\')){
                result.push_back(substr);
                //std::cout<<"split: "<<substr<<std::endl;
            }
        }
        else
            result.push_back(app);
    }
    //altrimenti aggiungi l'ultima stringa rimuovendo l'ultima parentesi
    else{
        app.erase(app.size()-1);
        result.push_back(app);
    }
    //ritorna il vettore calcolato
    return result;
}

    
bool swipl_interface::consult(std::string path){
    
    prolog_path = path;

    PlTermv args(1);
    args[0] = prolog_path.c_str();

    try{
        PlQuery pq("consult", args);

        if (pq.next_solution()) {
            std::cout << "SWIPL consult done"<<std::endl;
            return true;
        } else {
            std::cout << "SWIPL consult fail"<<std::endl;
            return false;
        }
    }
    catch(const PlException& e){
        return false;
    }
}

std::vector<std::string> swipl_interface::list2vector(std::string list){
    std::vector<std::string> lv;
    
    std::replace( list.begin(), list.begin()+1, '[', '(' );
    std::replace( list.begin()+list.size()-1, list.end(), ']', ')' );
    
    lv = instance2vector(list);
    lv.erase(lv.begin());

    for(auto i=0; i<lv.size(); i++)
        std::cout<<lv[i]<<std::endl;
    
    return lv;
}


PlTerm swipl_interface::string2term(std::string str, std::unordered_map<std::string, PlTerm> &swi_vars){
    std::vector<std::string> strv = instance2vector(str);
    
//    if(strv.size()<=0)
//        return PlCompound(); //not allowed!
    
    //std::cout<<"inside string2compound"<<std::endl;
    
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
 * post a generic goal to prolog and get the answer as a string 
 * containing a list of functors.
 * 
 * if no answer available, the returned string is empty
 * 
 */
std::string swipl_interface::query(std::string request){
    
    std::unordered_map<std::string, PlTerm> swi_vars;// = new std::unordered_map<std::string, term_t>();
    
    std::string from_prolog = "";
    
    std::vector<std::string> rv = instance2vector(request);
    
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

