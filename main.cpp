#include "parser.h"


void exp1()
{
    printfln("-------------------------------------EXP1---------------------------------");
    std::cout << "  Only DMRG" << std::endl;

    printfln("-------------------------------------PREAPRE OBJECTS----------------------");
    clock_t time0 = clock();

    seedRNG(1);
    auto sites = KondoHeisenberg( Args::global().getInt("L") );
    auto psi = prepareInitState(sites);
    auto H = prepareHamiltonian(sites);
    auto sweeps = prepareSweepClass();

    std::cout << "  Energy of initial state: " << real(innerC(psi,H,psi)) << std::endl;
    std::cout << "  N: " << calculateN(sites, psi) << std::endl;
    std::cout << "  N dublon: " << calculateNd(sites, psi) << std::endl;
    std::cout << "  Sz_0: " << calculateSz0(sites, psi) << std::endl;
    std::cout << "  Sz_1: " << calculateSz1(sites, psi) << std::endl;
    std::cout << "  Sz_01: " << calculateSz_01(sites, psi) << std::endl;
    std::cout << "  Time: " << (clock()-time0)/(double)CLOCKS_PER_SEC << " [s]" << std::endl;

    printfln("-------------------------------------DMRG---------------------------------");
    time0 = clock();

    dmrg(psi,H,sweeps);
    printfln("-------------------------------------DATA---------------------------------");
    std::cout << "  Energy after DMRG: " << real(innerC(psi,H,psi)) << std::endl;
    std::cout << "  N: " << calculateN(sites, psi) << std::endl;
    std::cout << "  N dublon: " << calculateNd(sites, psi) << std::endl;
    std::cout << "  Sz_0: " << calculateSz0(sites, psi) << std::endl;
    std::cout << "  Sz_1: " << calculateSz1(sites, psi) << std::endl;
    std::cout << "  Sz_01: " << calculateSz_01(sites, psi) << std::endl;
    std::cout << "  Time: " << (clock()-time0)/(double)CLOCKS_PER_SEC << " [s]" << std::endl;
}

void exp2()
{
    printfln("-------------------------------------EXP2---------------------------------");
    std::cout << "  Time evolution of init state" << std::endl;
    double t = Args::global().getReal("t");

    printfln("-------------------------------------PREAPRE OBJECTS----------------------");
    clock_t  time0 = clock();

    seedRNG(1);
    auto sites = KondoHeisenberg( Args::global().getInt("L") );
    auto psi = prepareInitState(sites);
    auto H = prepareHamiltonian(sites);
    auto sweeps = prepareSweepClass();

    std::cout << "Energy of initial state: " << std::real(innerC(psi,H,psi)) << std::endl;
    std::cout << "Time: " << (clock()-time0)/(double)CLOCKS_PER_SEC << " [s]" << std::endl;

    printfln("-------------------------------------TVDP in real time---------------------");
    time0 = clock();
    auto energy = tdvp(psi,H,-t*10.0,sweeps,{"DoNormalize",true,"Quiet",true,"NumCenter",2});
    std::cout << "Energy after TVDP: " << energy << std::endl;
    std::cout << "Time: " << (clock()-time0)/(double)CLOCKS_PER_SEC << " [s]" << std::endl;
}

void exp3()
{
    printfln("-------------------------------------EXP1---------------------------------");
    std::cout << "  DMRG with correlation" << std::endl;

    printfln("-------------------------------------PREAPRE OBJECTS----------------------");
    clock_t  time0 = clock();

    seedRNG(1);
    auto sites = KondoHeisenberg( Args::global().getInt("L") );
    auto psi = prepareInitState(sites);
    auto H = prepareHamiltonian(sites);
    auto sweeps = prepareSweepClass();

    std::cout << "  Energy of initial state: " << std::real(innerC(psi,H,psi)) << std::endl;

    std::cout << "  Time: " << (clock()-time0)/(double)CLOCKS_PER_SEC << " [s]" << std::endl;

    printfln("-------------------------------------DMRG---------------------------------");
    time0 = clock();

    dmrg(psi,H,sweeps,{"Silent",true});
    std::cout << "  Energy after DMRG: " << real(innerC(psi,H,psi)) << std::endl;
    std::cout << "  N: " << calculateN(sites, psi) << std::endl;
    std::cout << "  N dublon: " << calculateNd(sites, psi) << std::endl;
    std::cout << "  Sz_0: " << calculateSz0(sites, psi) << std::endl;
    std::cout << "  Sz_1: " << calculateSz1(sites, psi) << std::endl;
    std::cout << "  Sz_01: " << calculateSz_01(sites, psi) << std::endl;
    calculateCorrelationMatrixSz(sites,psi,"0");
    calculateCorrelationMatrixSz(sites,psi,"1");
    calculateCorrelationMatrixSz(sites,psi,"01");
    std::cout << "  Time: " << (clock()-time0)/(double)CLOCKS_PER_SEC << " [s]" << std::endl;
}

void run()
{
    switch (Args::global().getInt("exp",-1)){
    case 1 :
        exp1();
        break;
    case 2 :
        exp2();
        break;
    case 3 :
        exp3();
        break;
    default :
        std::cerr << "ERROR: epxeriment was not selected!" << std::endl;
        break;
    }
}

int main(int argc, char *argv[])
{
    readExperimentParameters(argc,argv);
    writeExperimentParameters();
    run();

    return 0;
}
