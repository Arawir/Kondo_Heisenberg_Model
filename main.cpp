#include "interface.h"
#include "model.h"
#include "tdvp.h"

int main(int argc, char *argv[])
{
    Experiments("Dmrg") = [](){
        ExpCon.addPoint("Initialization");

        seedRNG(1);
        auto sites = KH( getI("L") );
        auto psi = prepareInitState(sites);
        auto H = KHHamiltonian(sites,getI("L"),getD("thop"),getD("K"),getD("Jh"),getD("Mu"),getD("U"));
        auto sweeps = prepareSweepClass();

        std::cout << "  Energy: " << real(innerC(psi,H,psi)) << std::endl;
        std::cout << "  N: " << calculateN(sites, psi) << std::endl;
        std::cout << "  N dublon: " << calculateNd(sites, psi) << std::endl;
        std::cout << "  Sz_0: " << calculateSz0(sites, psi) << std::endl;
        std::cout << "  Sz_1: " << calculateSz1(sites, psi) << std::endl;
        std::cout << "  Sz_t: " << calculateSzt(sites, psi) << std::endl;

        ExpCon.addPoint("Starting DMRG");
        dmrg(psi,H,sweeps);

        ExpCon.addPoint("Output data");
        std::cout << "  Energy: " << real(innerC(psi,H,psi)) << std::endl;
        std::cout << "  N: " << calculateN(sites, psi) << std::endl;
        std::cout << "  N dublon: " << calculateNd(sites, psi) << std::endl;
        std::cout << "  Sz_0: " << calculateSz0(sites, psi) << std::endl;
        std::cout << "  Sz_1: " << calculateSz1(sites, psi) << std::endl;
        std::cout << "  Sz_t: " << calculateSzt(sites, psi) << std::endl;
    };

    Experiments("DmrgWithCorrelations") = [](){
        ExpCon.addPoint("Initialization");

        seedRNG(1);
        auto sites = KH( getI("L") );
        auto psi = prepareInitState(sites);
        auto H = KHHamiltonian(sites,getI("L"),getD("thop"),getD("K"),getD("Jh"),getD("Mu"),getD("U"));
        auto sweeps = prepareSweepClass();

        std::cout << "  Energy: " << std::real(innerC(psi,H,psi)) << std::endl;
        std::cout << "  N: " << calculateN(sites, psi) << std::endl;
        std::cout << "  N dublon: " << calculateNd(sites, psi) << std::endl;
        std::cout << "  Sz_0: " << calculateSz0(sites, psi) << std::endl;
        std::cout << "  Sz_1: " << calculateSz1(sites, psi) << std::endl;
        std::cout << "  Sz_t: " << calculateSzt(sites, psi) << std::endl;


        ExpCon.addPoint("Starting DMRG");
        dmrg(psi,H,sweeps,{"Silent",true});


        ExpCon.addPoint("Output data");
        std::cout << "  Energy: " << real(innerC(psi,H,psi)) << std::endl;
        std::cout << "  N: " << calculateN(sites, psi) << std::endl;
        std::cout << "  N dublon: " << calculateNd(sites, psi) << std::endl;
        std::cout << "  Sz_0: " << calculateSz0(sites, psi) << std::endl;
        std::cout << "  Sz_1: " << calculateSz1(sites, psi) << std::endl;
        std::cout << "  Sz_t: " << calculateSzt(sites, psi) << std::endl;
        calculateCorrelationMatrixSz(sites,psi, "Sz0");
        calculateCorrelationMatrixSz(sites,psi, "Sz1");
        calculateCorrelationMatrixSz(sites,psi, "SzSz");
        calculateCorrelationMatrixSz(sites,psi, "SmSp");
        calculateCorrelationMatrixSz(sites,psi, "SpSm");
    };


    Params.add("thop","double","0.0");
    Params.add("U","double","0.0");
    Params.add("K","double","0.0");
    Params.add("Jh","double","0.0");
    Params.add("Mu","double","0.0");

    Params.add("L","int","4");
    Params.add("PBC","bool","0");

    Params.add("Silent","bool","1");
    Params.add("cutoff","double","1E-6");
    Params.add("sweeps","int","4");
    Params.add("minDim","int","1");
    Params.add("maxDim","int","100");
    Params.add("niter","int","10");
    Params.add("state","string","ud-du");

    Params.add("ConserveN","bool","0");
    Params.add("ConserveSz","bool","0");

    Params.add("exp","string","1");

    Params.set(argc,argv);
    Experiments.run();

    return 0;
}
