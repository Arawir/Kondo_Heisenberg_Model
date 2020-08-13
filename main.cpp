#include "parser.h"


void exp1()
{
    printfln("------------------------------ Initial State -----------------------------");
    clock_t timecpu = clock();

    seedRNG(1);
    double thop = Args::global().getReal("thop",0);
    double K = Args::global().getReal("K",0);
    double U = Args::global().getReal("U",0);
    double Jh = Args::global().getReal("Jh",0);
    double Mu = Args::global().getReal("Mu",0);
    int L = Args::global().getInt("L",4);

    auto sites = KondoHeisenberg( Args::global().getInt("L") );
    auto psi = prepareInitState(sites);
    auto H = prepareHamiltonianKH(sites,L,thop,K,U,Jh,Mu);
    auto sweeps = prepareSweepClass();

    std::cout << "  Energy: " << real(innerC(psi,H,psi)) << std::endl;
    std::cout << "  N: " << calculateN(sites, psi) << std::endl;
    std::cout << "  N dublon: " << calculateNd(sites, psi) << std::endl;
    std::cout << "  Sz_0: " << calculateSz0(sites, psi) << std::endl;
    std::cout << "  Sz_1: " << calculateSz1(sites, psi) << std::endl;
    std::cout << "  Sz_t: " << calculateSzt(sites, psi) << std::endl;
    std::cout << "  Time: " << (clock()-timecpu)/(double)CLOCKS_PER_SEC << " [s]" << std::endl;

    printfln("------------------------- Ground State from DMRG -------------------------");
    timecpu = clock();

    dmrg(psi,H,sweeps);
    std::cout << "  Energy: " << real(innerC(psi,H,psi)) << std::endl;
    std::cout << "  N: " << calculateN(sites, psi) << std::endl;
    std::cout << "  N dublon: " << calculateNd(sites, psi) << std::endl;
    std::cout << "  Sz_0: " << calculateSz0(sites, psi) << std::endl;
    std::cout << "  Sz_1: " << calculateSz1(sites, psi) << std::endl;
    std::cout << "  Sz_t: " << calculateSzt(sites, psi) << std::endl;
    std::cout << "  Time: " << (clock()-timecpu)/(double)CLOCKS_PER_SEC << " [s]" << std::endl;
}

void exp2()
{
    double time = Args::global().getReal("time");

    printfln("------------------------------ Initial State -----------------------------");
    clock_t  timecpu = clock();

    seedRNG(1);
    double thop = Args::global().getReal("thop",0);
    double K = Args::global().getReal("K",0);
    double U = Args::global().getReal("U",0);
    double Jh = Args::global().getReal("Jh",0);
    double Mu = Args::global().getReal("Mu",0);
    int L = Args::global().getInt("L",4);

    auto sites = KondoHeisenberg( Args::global().getInt("L") );
    auto psi = prepareInitState(sites);
    auto H = prepareHamiltonianKH(sites,L,thop,K,U,Jh,Mu);
    auto sweeps = prepareSweepClass();

    std::cout << "Energy: " << std::real(innerC(psi,H,psi)) << std::endl;
    std::cout << "Time: " << (clock()-timecpu)/(double)CLOCKS_PER_SEC << " [s]" << std::endl;

    printfln("------------------ Time evolution of initial state ----------------");
    timecpu = clock();
    auto energy = tdvp(psi,H,-time*10.0,sweeps,{"DoNormalize",true,"Quiet",true,"NumCenter",2});
    std::cout << "Energy after TVDP: " << energy << std::endl;
    std::cout << "Time: " << (clock()-timecpu)/(double)CLOCKS_PER_SEC << " [s]" << std::endl;
}

void exp3()
{
    double time = Args::global().getReal("time");

    printfln("------------------------------ Initial State -----------------------------");

    seedRNG(1);
    double thop = Args::global().getReal("thop",0);
    double K = Args::global().getReal("K",0);
    double U = Args::global().getReal("U",0);
    double Jh = Args::global().getReal("Jh",0);
    double Mu = Args::global().getReal("Mu",0);

    double thop2 = Args::global().getReal("thop2",thop);
    double K2 = Args::global().getReal("K2",K);
    double U2 = Args::global().getReal("U2",U);
    double Jh2 = Args::global().getReal("Jh2",Jh);
    double Mu2 = Args::global().getReal("Mu2",Mu);

    int L = Args::global().getInt("L",4);

    auto sites = KondoHeisenberg( Args::global().getInt("L") );
    auto psi = prepareInitState(sites);
    auto H1 = prepareHamiltonianKH(sites,L,thop,K,U,Jh,Mu);
    auto H2 = prepareHamiltonianKH(sites,L,thop2,K2,U2,Jh2,Mu2);
    auto sweeps = prepareSweepClass();

    printfln("------------------ Calculate H0 ground state ----------------");
    dmrg(psi,H1,sweeps);

    printfln("------------------ Time evolution of H0 ground state with H1 hamiltonian ----------------");
    tdvp(psi,H2,-time*10.0,sweeps,{"DoNormalize",true,"Quiet",true,"NumCenter",2});


}

void exp4()
{
    printfln("------------------------------ Initial State -----------------------------");
    clock_t  timecpu = clock();

    seedRNG(1);
    double thop = Args::global().getReal("thop",0);
    double K = Args::global().getReal("K",0);
    double U = Args::global().getReal("U",0);
    double Jh = Args::global().getReal("Jh",0);
    double Mu = Args::global().getReal("Mu",0);
    int L = Args::global().getInt("L",4);

    auto sites = KondoHeisenberg( Args::global().getInt("L") );
    auto psi = prepareInitState(sites);
    auto H = prepareHamiltonianKH(sites,L,thop,K,U,Jh,Mu);
    auto sweeps = prepareSweepClass();

    std::cout << "  Energy: " << std::real(innerC(psi,H,psi)) << std::endl;

    std::cout << "  Time: " << (clock()-timecpu)/(double)CLOCKS_PER_SEC << " [s]" << std::endl;

    printfln("------------------- Ground State correlations from DMRG ------------------");
    timecpu = clock();

    dmrg(psi,H,sweeps,{"Silent",true});
    std::cout << "  Energy: " << real(innerC(psi,H,psi)) << std::endl;
    std::cout << "  N: " << calculateN(sites, psi) << std::endl;
    std::cout << "  N dublon: " << calculateNd(sites, psi) << std::endl;
    std::cout << "  Sz_0: " << calculateSz0(sites, psi) << std::endl;
    std::cout << "  Sz_1: " << calculateSz1(sites, psi) << std::endl;
    std::cout << "  Sz_t: " << calculateSzt(sites, psi) << std::endl;
    calculateCorrelationMatrixSz(sites,psi,"0");
    calculateCorrelationMatrixSz(sites,psi,"1");
    calculateCorrelationMatrixSz(sites,psi,"t");
    std::cout << "  Time: " << (clock()-timecpu)/(double)CLOCKS_PER_SEC << " [s]" << std::endl;
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
