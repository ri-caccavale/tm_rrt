/*
 * File:   main.cpp
 * Author: hargalaten
 *
 * Created on 4 dicembre 2012, 13.59
 */

#include "Eclipseclp_interface.h"

    
Eclipseclp_interface::Eclipseclp_interface(std::string path_to_eclipseclp_directory){
    
    eclpse_path = path_to_eclipseclp_directory;
    
    ec_set_option_ptr(EC_OPTION_ECLIPSEDIR, (void *) eclpse_path.c_str());
    
    ec_init();
}

bool Eclipseclp_interface::close(){
    ec_cleanup();
    
    return true;
}
    
bool Eclipseclp_interface::open(std::string path_to_file){
    
    sourcefile_path = path_to_file;

    //read long time memory
    std::string sourcefile_path_goal = "[\'" + sourcefile_path + "\']";

    ec_post_string(sourcefile_path_goal.c_str());
    EC_resume();
    
    return true; //cant fail?
}

//LTM
EC_word Eclipseclp_interface::writeList(std::vector<std::string> vec){
    EC_word l=nil();

    for(int i=0;i<vec.size();i++){
        char *cstr = new char[vec[i].length() + 1];
        strcpy(cstr, vec[i].c_str());
        l=list(EC_atom(cstr),l);
        delete [] cstr;
    }
    return l;
}

//LTM
std::vector<std::string> Eclipseclp_interface::readList(EC_word list){
    char *buf;
    std::vector<std::string> app;
    EC_word head,tail;
    if(list.is_nil()){
        list.is_list(head,tail);
        head.is_string(&buf);
        //std::cout<<"read: "<<buf<<"\n";
        app=readList(tail);
        app.push_back(std::string(buf));
    }
    return app;
}

//LTM
void Eclipseclp_interface::printList(EC_word list)
{
    char *buf;
    EC_word head,tail;

    list.is_list(head,tail);
    head.is_string(&buf);
    if(tail.is_nil())
    {
        printf(" %s", buf);
        printList(tail);
    }
    else printf(" %s\n", buf);
}

//LTM
//trasforma in stringa un generico funtore restituito da ECLIPSE
std::string Eclipseclp_interface::functor2string(EC_word f){
    char* buf;
    char* subbuf;
    char* st;
    long stl;
    std::string result,subresult;
    std::stringstream ss;
    EC_word arg, parHead,parTail;
    EC_functor fun,subfun;
    EC_atom atm;
    double d;
    int i;

    //se il funtore ha arietà maggiore di zero
    if(f.arity()!=0)
    {
        //estrai in nome del funtore
        f.functor(&fun);
        buf=fun.name();
        result.append(buf);

        //se è un funtore punto (ie. del tipo schema.attributo)
        if(result=="." && f.arity()==2)
        {
            //svuota la stringa
            result="";
            //prendi lo schema
            f.arg(1,arg);
            result.append(functor2string(arg));
            //inseriscilo nella stringa seguito dal punto
            result=result+".";
            //prendi il metodo
            f.arg(2,arg);
            result.append(functor2string(arg));
        }

        else if(result=="-" && f.arity()==1){
            //prendi il parametro ed aggiungilo alla stringa
            f.arg(1,arg);
            //arg.is_string(&buf);
            result.append(functor2string(arg));
        }
        //altrimenti è uno schema avente parametri (eg. goto(X))
        else
        {
            //std::cout<<"funtore: "<<result<<"\n";
            //aggiungi una parentesi aperta
            result.append("(");
            //per ogni argomento
            for(i=1;i<=f.arity();i++)
            {
                //inserisci l'argomento nella stringa
                f.arg(i,arg);
                //std::cout<<"isFun: "<<arg.functor(&subfun) << " isList:"<<!arg.is_list(parHead, parTail)<<"\n";
                //se l'i-esimo argomento è una lista ...e non è un funtore!!
                //  (il punto è sia funtore che lista -.-" quindi per scartarlo va controllato)
                if (arg.functor(&subfun) && !arg.is_list(parHead, parTail)) {
                    //mantieni la sintasi di lista
                    result.append("[");
                    result.append(functor2string(parHead));
                    //std::cout<<"head: "<<functor2string(parHead)<<"\n";
                    while (parTail.is_nil()) {
                        result.append(",");
                        parTail.is_list(parHead, parTail);
                        result.append(functor2string(parHead));
                        //std::cout<<"newhead: "<<functor2string(parHead)<<"\n";
                    }
                    result.append("],");
                    //std::cout << "list: " << result << "\n";
                }
                //altrimenti è un funtore
                else
                    result=result+functor2string(arg)+",";

                //sleep(0.1);
                //std::cout<<result<<"\n";
            }
            result.replace(result.size()-1,1,")");

            //std::cout<<"end functor: "<<result<<"\n";
        }
    }
    //altrimenti non ha argomenti
    else
    {
        //se è un atomo
        if(!f.is_atom(&atm)){
//            std::cout<<"è atomo: "<<atm.name()<<"\n";
            result.append(atm.name());
        }
        //altrimenti, se è una stringa
        else if(!f.is_string(&buf)){
//            std::cout<<"è stringa: "<<buf<<"\n";
            std::string element(buf);

            //se la stringa è un TRUE o FALSE aggiungilo
            if(element=="TRUE" || element=="FALSE")
                result.append(element);
            //altrimenti aggiungi gli apici
            else
                result.append("\""+element+"\"");

        }
//        else if(f.is_string(&buf)){
//            std::cout<<"dentro\n";
//            result.append(buf);
//            result="\""+result+"\"";
//        }
        else if(!f.is_double(&d)){
            ss<<d;
            result.append(ss.str());
        }
    }
    //restituisci al stringa
    return result;
}

std::vector<std::string> Eclipseclp_interface::query(std::string request){
    EC_ref elemFromEclipse;
    EC_word elem,elemList;

    std::vector<std::string> response;

    EC_resume();
    //send the request
    post_goal(request.c_str());
    EC_resume();

    std::cout<<"\t Query: "<<request<<std::endl;

    //get the response (list of EC_word)
    if (EC_resume(EC_atom("ok"), elemFromEclipse) == EC_yield) {
        elemList = EC_word(elemFromEclipse);
        //get the elements of the response
        while (elemList.is_nil()) {
            elemList.is_list(elem, elemList);
            std::cout<<"\t\t "<<functor2string(elem)<<std::endl;
            response.push_back(functor2string(elem));
        }
    }
    //dismiss ECLIPSE
    EC_resume();

    return response;
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
 */
std::vector<std::string> instance2vector(std::string schemaInstance){
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
        else if(c=='('){
            count++;
            app=app+c;
        }
        else if(c==')'&& count!=0){
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
    if(isAtom) result.push_back(app);
    //altrimenti aggiungi l'ultima stringa rimuovendo l'ultima parentesi
    else{
        app.erase(app.size()-1);
        result.push_back(app);
    }
    //ritorna il vettore calcolato
    return result;
}
