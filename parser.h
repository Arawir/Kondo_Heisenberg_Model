#ifndef PARSER
#define PARSER

#include <iostream>
#include <iomanip>
#include <ctime>
#include <cstdlib>
#include <fstream>
#include <vector>
#include <string>
#include <cassert>

#include "operators.h"

std::vector<std::string> parseInitState(std::string str)
{
    int L = Args::global().getInt("L");

    std::vector<std::string> state{};
    int i=0;
    state.push_back("");
    while(str[i]!='\0'){
        if(str[i]=='-'){
            state.push_back("");
        } else {
            state.back().append(std::string{str[i]});
        }
        i++;
    }

    if(state.size()==0){
        std::cerr << "  Init state is not delcared!" << std::endl;
    }else if(state.size()<(uint)L){
        uint templateLength = state.size()-1;
        std::cout << "  Init state deduced" << std::endl;
        for(int i=state.size()-1; i<L; i++){
            state.push_back( state.at(i-templateLength) );
        }
    }
    return state;
}

MPS prepareInitState(BasicSiteSet<KondoHeisenbergSite> &sites)
{
    std::vector<std::string> initState = parseInitState(Args::global().getString("state"));
    int L = Args::global().getInt("L");

    auto state = InitState(sites);
    for(int i=1; i<=L; i++){
       state.set(i,initState[i-1]);
    }

    return MPS(state);
}

std::string paramName(std::string argv)
{
    size_t pos = argv.find("=");
    if(pos==std::string::npos){
        std::cerr << "WARNING: cannot find \"=\" in parameter (" << argv << ")" << std::endl;
    }
    return argv.substr(0,pos);
}

void setParam(std::string type, std::string argv)
{
    size_t pos = argv.find("=");
    if(pos!=std::string::npos){
        if(type=="bool"){
            Args::global().add(paramName(argv), (bool)atoi( argv.substr(pos+1).c_str() ));
        } else if(type=="string"){
            Args::global().add(paramName(argv), argv.substr(pos+1) );
        } else if(type=="int"){
            Args::global().add(paramName(argv), atoi( argv.substr(pos+1).c_str() ));
        } else if(type=="double"||type=="real"){
            Args::global().add(paramName(argv), atof( argv.substr(pos+1).c_str() ));
        } else {
            std::cerr << "WARNING: unknown type (" << argv << ")" << std::endl;
        }
    } else {
        std::cerr << "WARNING: cannot find \"=\" in parameter (" << argv << ")" << std::endl;
    }
}

void readExperimentParameters(int argc, char *argv[])
{
    for(int i=1; i<argc; i++){
        if( paramName(argv[i]) == "L" ){ setParam("int", argv[i]); }
        else if( paramName(argv[i]) == "PBC" ){ setParam("bool", argv[i]); }
        else if( paramName(argv[i]) == "thop" ){ setParam("double", argv[i]); }
        else if( paramName(argv[i]) == "U" ){ setParam("double", argv[i]); }
        else if( paramName(argv[i]) == "K" ){ setParam("double", argv[i]); }
        else if( paramName(argv[i]) == "Jh" ){ setParam("double", argv[i]); }
        else if( paramName(argv[i]) == "Mu" ){ setParam("double", argv[i]); }
        else if( paramName(argv[i]) == "time" ){ setParam("double", argv[i]); }
        else if( paramName(argv[i]) == "exp" ){ setParam("int", argv[i]); }
        else if( paramName(argv[i]) == "conSz" ){ setParam("bool", argv[i]); }
        else if( paramName(argv[i]) == "conN" ){ setParam("bool", argv[i]); }
        else if( paramName(argv[i]) == "state" ){ setParam("string", argv[i]); }
        else if( paramName(argv[i]) == "maxDim" ){ setParam("int", argv[i]); }
        else if( paramName(argv[i]) == "minDim" ){ setParam("int", argv[i]); }
        else if( paramName(argv[i]) == "sweeps" ){ setParam("int", argv[i]); }
        else if( paramName(argv[i]) == "cutoff" ){ setParam("double", argv[i]); }
        else if( paramName(argv[i]) == "Silent" ){ setParam("bool", argv[i]); }
        else {
            std::cerr << "WARNING: unknown argument (" << argv[i] << ")" << std::endl;
        }
    }
}

void writeExperimentParameters()
{
    std::cout << "---------------------------- System parameters ---------------------------" << std::endl;

    std::cout << "  L = " << Args::global().getInt("L",4) << std::endl;
    std::cout << "  thop = " << Args::global().getReal("thop",0) << std::endl;
    std::cout << "  U = " << Args::global().getReal("U",0) << std::endl;
    std::cout << "  K = " << Args::global().getReal("K",0) << std::endl;
    std::cout << "  Jh = " << Args::global().getReal("Jh",0) << std::endl;
    std::cout << "  Mu = " << Args::global().getReal("Mu",0) << std::endl;
    std::cout << "  PBC = " << Args::global().getBool("PBC",0) << std::endl;
    std::cout << "  Conserve Sz = " << Args::global().getBool("conserveSz",0) << std::endl;
    std::cout << "  Conserve N = " << Args::global().getBool("conserveN",0) << std::endl;
    std::cout << "  sweeps = " << Args::global().getInt("sweeps",4) << std::endl;
    std::cout << "  minDim = " << Args::global().getInt("minDim",1) << std::endl;
    std::cout << "  maxDim = " << Args::global().getInt("maxDim",100) << std::endl;
    std::cout << "  cutoff = " << Args::global().getReal("cutoff",1E-6) << std::endl;
    std::cout << "  Init state = " << Args::global().getString("state","") << std::endl;
}

#endif // PARSER

