#include "../include/inicpp.h"

using namespace std;


vector<double> read_ini(string path){
    static int ini_section;
    vector<double> vars(18);
    vector<string> names(18);
       
    names= {"dm21", "dm31", "theta12", "theta13", "theta23",
            "param", "xi1", "xi2", "xi3", "xit1", "xit2", "xit3", "dn1", "dn2", "dn3", 
            "zon", "vs", "Tend"};

    ini::IniFile inif;
    inif.load(path);

    ini_section= 0;
    for(const auto &sectionPair : inif){
        const string &sectionName       = sectionPair.first;
        if (sectionName == "Ini"){
            ini_section++;
            for(const auto &fieldPair : sectionPair.second){
                const string &fieldName     = fieldPair.first;
                const ini::IniField &field  = fieldPair.second;
                for (int i=0; i<names.size(); i++){
                    if (fieldName == names[i]){
                        vars[i] = field.as<double>(); 
                    }
                }
            }
        }
    }


    if (ini_section != 0){
        cout << "succesfully read doubles of .ini file" << endl;
        return vars;
    }
    else{
        cout << "The .ini file needs to contain the section [Ini]" << endl;
        cout << "It is the only section which will be read" << endl;
        exit(0);
    }
    
}




vector<string> read_ini_filename(string path){
    static int ini_section;
    vector<string> vars(2);
    vector<string> names(2);
       
    names= {"filename", "progressbar"};

    ini::IniFile inif;
    inif.load(path);

    ini_section= 0;
    for(const auto &sectionPair : inif){
        const string &sectionName       = sectionPair.first;
        if (sectionName == "Ini"){
            ini_section++;
            for(const auto &fieldPair : sectionPair.second){
                const string &fieldName     = fieldPair.first;
                const ini::IniField &field  = fieldPair.second;
                for (int i=0; i<names.size(); i++){
                    if (fieldName == names[i]){
                        vars[i] = field.as<string>(); 
                    }
                }
            }
        }
    }


    if (ini_section != 0){
        cout << "succesfully read strings of .ini file" << endl;
        return vars;
    }
    else{
        cout << "The .ini file needs to contain the section [Ini]" << endl;
        cout << "It is the only section which will be read" << endl;
        exit(0);
    }
    
}

