#ifndef INTERFACE
#define INTERFACE

#include "itensor/all.h"
#include "kondo_heisenberg.h"
#include <ctime>
#include <functional>
#include <cstdio>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
#include <array>

typedef std::complex<double> cpx;
#define im std::complex<double>{0.0,1.0}


inline void massert(std::string info)
{
    std::cerr << info << std::endl;
    assert(false);
}

using namespace itensor;

double getD(std::string data)
{
    return Args::global().getReal(data);
}

int getI(std::string data)
{
    return Args::global().getInt(data);
}

bool getB(std::string data)
{
    return Args::global().getBool(data);
}

std::string getS(std::string data)
{
    return Args::global().getString(data);
}

struct Param{
    std::string name;
    std::string type;
    std::string value;
    bool wasSet=false;

    int getInt(){
        if(type!="int"){
            std::cerr << "WARNING: wrong type conversion (" + name + ")!" << std::endl;
        }
        return atoi(value.c_str());
    }

    double getDouble(){
        if(type!="double"){
            std::cerr << "WARNING: wrong type conversion (" + name + ")!" << std::endl;
        }
        return atof(value.c_str());
    }

    bool getBool(){
        if(type!="bool"){
            std::cerr << "WARNING: wrong type conversion (" + name + ")!" << std::endl;
        }
        return (bool)atoi(value.c_str());
    }
};

class Parameters{
private:
    std::vector<Param> params;

public:
    Parameters(){ }
    ~Parameters() = default;

    void add(std::string name, std::string type, std::string value){
        params.push_back( Param{name, type, value} );
    }

    void set(int argc, char *argv[]){
        for(int i=1; i<argc; i++){
            set(argv[i]);
        }
        addArgs();
        write();
    }

    bool getBool(std::string name){
        return param(name)->getBool();
    }
    int getInt(std::string name){
        return param(name)->getInt();
    }
    double getDouble(std::string name){
        return param(name)->getDouble();
    }

    void write(){
        std::cout << "---------------------------- System parameters ---------------------------" << std::endl;
        for(auto& param : params){
            std::cout << "  " << param.name << " -> ";
            std::cout << "  " << param.value << " ";
            if(!param.wasSet){ std::cout << "(default)"; }
            std::cout << std::endl;
        }
    }

private:
    void set(std::string data){
        size_t pos = data.find("=");
        param( data.substr(0,pos) )->value = data.substr(pos+1);
        param( data.substr(0,pos) )->wasSet = true;
    }

    void addArgs(){
        for(auto &param : params){
            if(param.type=="bool"){
                Args::global().add(param.name, param.getBool());
            } else if(param.type=="string"){
                Args::global().add(param.name, param.value);
            } else if(param.type=="int"){
                Args::global().add(param.name, param.getInt());
            } else if(param.type=="double"){
                Args::global().add(param.name, param.getDouble());
            }
        }
    }


    Param* param(std::string name){
        for(auto &param : params){
            if(param.name == name) return &param;
        }
        massert("ERROR: Wrong parameter name (" +name + ")!");
        return nullptr;
    }
} Params;


std::vector<std::string> separateSubBlocs(std::string str)
{
    std::vector<std::string> out{};
    int i=0;
    out.push_back("");
    while(str[i]!='\0'){
        if(str[i]=='-'){
            out.push_back("");
        } else {
            out.back().append(std::string{str[i]});
        }
        i++;
    }
    return out;
}

int identifyMultiplayer(std::string data)
{
    if(data.find("L/")!=std::string::npos){
        return getI("L")/atoi( data.substr(data.find("L/")+2).c_str() );
    }
    return atoi(data.c_str());
}

std::tuple<int,std::string> separateMultiplayer(std::string data)
{
    int multiplayer=1;
    std::string state=data;
    if(data.find("*")!=std::string::npos){
        multiplayer = identifyMultiplayer( data.substr(0,data.find("*")) );
        state = data.substr(data.find("*")+1);
    }

    return std::make_tuple( multiplayer,state );
}

std::vector<std::string> expandSubBlocs(std::vector<std::string> subBlocs)
{
    std::vector<std::string> out{};

    for(auto subBlock : subBlocs){
        auto [n,state] = separateMultiplayer(subBlock);
        for(int i=0; i<n; i++){
            out.push_back(state);
        }
    }
    return out;
}

std::vector<std::string> parseInitState(std::string str)
{
    int L = Params.getInt("L");

    std::vector<std::string> subBlocs = separateSubBlocs(str);
    std::vector<std::string> state = expandSubBlocs(subBlocs);


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



itensor::MPS prepareInitState(auto &sites)
{
    if( (Args::global().getString("state")=="random")||(Args::global().getString("state")=="r")){
        auto state = itensor::InitState(sites);
        auto s2 = itensor::randomMPS(state);
        s2.orthogonalize();
        s2.normalize();
        return s2;
    }

    std::vector<std::string> initState = parseInitState(Args::global().getString("state"));
    int L = Params.getInt("L");



    auto state = itensor::InitState(sites);
    for(int i=1;i<=L;i++){
        state.set(i,initState[i-1]);
        std::cout << initState[i-1] << "-";
    }
    std::cout << std::endl;

    return itensor::MPS(state);
}

struct Experiment
{
    std::string name;
    std::function<void(void)> experiment;

    Experiment(std::string nName) :
        name{nName}
    {

    }

    void operator = (std::function<void(void)> exp)
    {
        experiment = exp;
    }
};


std::string exec(const char* cmd) {
    std::array<char, 128> buffer;
    std::string result;
    std::unique_ptr<FILE, decltype(&pclose)> pipe(popen(cmd, "r"), pclose);
    if (!pipe) {
        throw std::runtime_error("popen() failed!");
    }
    while (fgets(buffer.data(), buffer.size(), pipe.get()) != nullptr) {
        result += buffer.data();
    }
    return result;
}

///////////////////////////////////////////

struct Observable
{
    std::string name;
    std::vector<MPO> matrix;
    bool wasGenerated=false;
    bool useMultMpos = false;
    std::function<MPO(const KH &sites)> generateMPO;
    std::function<std::vector<MPO>(const KH &sites)> generateMPOs;
    KH *sites = nullptr;

    void operator = (std::function<MPO(const KH &sites)> generate)
    {
        generateMPO = generate;
    }

    void operator = (std::function<std::vector<MPO>(const KH &sites)> generate)
    {
        useMultMpos=true;
        generateMPOs = generate;
    }

    void operator =  (MPO &nMatrix)
    {
        matrix.push_back( nMatrix );
        wasGenerated=true;
    }

    void generateIfNeeded()
    {
        if(wasGenerated==false){
            wasGenerated=true;
            if(!useMultMpos){
                matrix.push_back( generateMPO(*sites) );
            } else {
                matrix = generateMPOs(*sites);
            }
        }
    }

    std::vector<cpx> expectedValue(MPS psi)
    {
        std::vector<cpx> out;
        for(auto &m : matrix){
            out.push_back( innerC(psi,m,psi) );
        }
        return out;
    }
};

enum class oMode{
    a, b, c
};

class ObservableContainer
{
public:
    std::vector<Observable> observables;
    oMode mode = oMode::a;
    clock_t  time0, timeLast;

public:
    ObservableContainer() :
        time0{ clock() }
        , timeLast{ clock() }
    { }
    ~ObservableContainer() = default;

    Observable& operator () (std::string name)
    {
        for(auto& obs : observables){
            if(obs.name == name){ return obs; }
        }
        observables.push_back( Observable{name} );
        return observables.back();
    }

    template<typename T,typename... Args>
    void calc(MPS &psi, T val, Args... args)
    {
        std::cout << " ";
        writeObservableValue(val,psi);
        if(mode==oMode::c){ std::cout << std::endl; }
        calc(psi, args...) ;
    }

    template<typename... Args>
    void calc(MPS &psi, oMode val, Args... args)
    {
        mode = val;
        calc(psi, args...) ;
    }

    void calc(MPS &psi)
    {
        std::cout << std::endl;
    }

    void setSites(KH &sites)
    {
        for(auto& obs : observables){
            obs.sites = &sites;
        }
    }


private:
    Observable* observable(std::string name)
    {
        for(auto& obs : observables){
            if(obs.name == name){ return &obs; }
        }
        massert("ERROR: unknown observable name!");
        return nullptr;
    }
    bool exists(std::string name)
    {
        for(auto& obs : observables){
            if(obs.name == name){ return true; }
        }
        return false;
    }

    void initObservablesIfNeeded()
    {
        for(auto& obs : observables){
            obs.generateIfNeeded();
        }
    }

    void writeObservableValue(std::string name, MPS &psi)
    {
        if( exists(name) == true){
            observable(name)->generateIfNeeded();
            if(mode==oMode::b || mode==oMode::c){ std::cout << name << "= "; }
            for(auto val : observable(name)->expectedValue(psi)){
                std::cout << val.real() << " ";
            }
        } else if(name=="dim"){
            std::cout << "dim= " << maxLinkDim(psi);
        } else if(name=="mem"){
            if(getB("PBSenable")==true){
                std::string command = "qstat -fx " + std::to_string(getI("PBSjobid")) + " | grep  'vmem' | awk '{print $3}'";
                std::cout << "mem= " << exec(command.c_str());
            } else {
                std::cout << "mem= " << "unknown";
            }
        }  else if(name=="rtime"){
            std::cout << "rtime= " << (clock()-time0)/(double)CLOCKS_PER_SEC;
        } else {
            std::cout << name;
        }
    }
    void writeObservableValue(double val, MPS &psi){ std::cout << val; }
    void writeObservableValue(cpx val, MPS &psi){ std::cout << val; }
    void writeObservableValue(int val, MPS &psi){ std::cout << val; }
};

/////////////////////////////////////////////////
class Controller : public ObservableContainer
{

public:

    Controller() :
        ObservableContainer{ }
    {

    }
    ~Controller() = default;

    void addPoint(std::string text)
    {
        writePointInfo(text);
    }

private:
    void writePointInfo(std::string text)
    {
        std::cout << "--------------------";
        std::cout << text;
        for(int i=0; i<60-(int)text.size(); i++){
           std::cout << "-";
        }

        std::cout << "(" << (clock()-time0)/(double)CLOCKS_PER_SEC << " [s],"
                  << (clock()-timeLast)/(double)CLOCKS_PER_SEC << " [s])";
        std::cout << std::endl;
        timeLast = clock();
    }
} ExpCon;

////////////////////////////////////////////////
class ExperimentsClass
{
private:
    std::vector<Experiment> experiments;
public:
    ExperimentsClass() :
        experiments {  }
    {

    }
    ~ExperimentsClass() = default;

    Experiment& operator () (std::string name)
    {
        for(auto& exp : experiments){
            if(exp.name == name){ return exp; }
        }
        experiments.push_back( Experiment{name} );
        return experiments.back();
    }

    void run()
    {
        ExpCon.addPoint("Start");
        experiment(getS("exp"))->experiment();
        ExpCon.addPoint("Finish");
    }

    Experiment* experiment(std::string name)
    {
        for(auto& exp : experiments){
            if(exp.name == name){ return &exp; }
        }
        massert("ERROR: unknown experiment name!");
        return nullptr;
    }

} Experiments;

#endif // INTERFACE



