#include <iostream>
#include <stdexcept>

#include "cSaturnModel.h"
#include "Config.h"
#include "tinyxml2.h"

void cSaturnModel::LoadConfig(const char *filename){
    XMLDocument doc;
    XMLError err = doc.LoadFile(filename);
    cout << "reading config-xml-file in cSaturnModel::LoadConfig ======= error message = err " << err << "\n\n";
    if(err){
        doc.PrintError();
        throw std::invalid_argument("couldn't load config file");
    }
    XMLElement* atsat = doc.FirstChildElement("atsat");
    if(!atsat){
        return;
    }
    XMLElement* elem_common = doc.FirstChildElement("atsat")->FirstChildElement("elem_common");
    if(!elem_common){
        return;
    }
    XMLElement* elem_atmosphere = doc.FirstChildElement("elem_common")->FirstChildElement("elem_saturn");
    if(!elem_saturn){
        return;
    }
    #include "SaturnLoadConfig.cpp.inc"
}
