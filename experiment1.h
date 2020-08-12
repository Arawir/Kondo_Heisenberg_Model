#ifndef EXPERIMENT1
#define EXPERIMENT1

#include "parser.h"

namespace Experiment1
{
    MPO prepareHamiltonian(BasicSiteSet<KondoHeisenbergSite> &sites)
    {
        double t00 = Args::global().getReal("t00",0);
        double K = Args::global().getReal("K",0);
        double U = Args::global().getReal("U",0);
        double Jh = Args::global().getReal("Jh",0);
        double Mu = Args::global().getReal("Mu",0);
        int L = Args::global().getInt("L",4);

        auto ampo = AutoMPO(sites);
        for(int j=1; j<L; j++){
            ampo += t00,"cT_0,u",j,"c_0,u",j+1;
            ampo += t00,"c_0,u",j,"cT_0,u",j+1;
            ampo += t00,"cT_0,d",j,"c_0,d",j+1;
            ampo += t00,"c_0,d",j,"cT_0,d",j+1;

            ampo += K/2,"s+_1",j,"s-_1",j+1;
            ampo += K/2,"s-_1",j,"s+_1",j+1;
            ampo += K,"sz_1",j,"sz_1",j+1;
        }

        if(Args::global().getBool("PBC",0)){
            ampo += t00,"cT_0,u",L,"c_0,u",1;
            ampo += t00,"c_0,u",L,"cT_0,u",1;
            ampo += t00,"cT_0,d",L,"c_0,d",1;
            ampo += t00,"c_0,d",L,"cT_0,d",1;

            ampo += K/2,"s+_1",L,"s-_1",1;
            ampo += K/2,"s-_1",L,"s+_1",1;
            ampo += K,"sz_1",L,"sz_1",1;
        }

        for(int j=1; j<=L; j++){
            ampo += U,"n_0,ud",j;
            ampo += -2.0*Jh,"s_01",j;
            ampo += Mu,"n_01",j;
        }

        return toMPO(ampo);
    }

    void calculateCorrelationMatrixSz(const BasicSiteSet<KondoHeisenbergSite> &sites, const MPS &psi, std::string type)
    {
        int L = Args::global().getInt("L");

        if(type=="0"){ std::cout << "Correlation matrix <Sz0_i Sz0_j>" << std::endl; }
        if(type=="1"){ std::cout << "Correlation matrix <Sz1_i Sz1_j>" << std::endl; }
        if(type=="01"){ std::cout << "Correlation matrix <Sz01_i Sz01_j>" << std::endl; }
        for(int i=1; i<=L; i++){
            std::cout << "  ";
            std::cout.width(3);
            std::cout << i << " | ";
            for(int j=1; j<=L; j++){
                std::cout.width(8);
                std::cout.precision(4);
                if(type=="0"){ std::cout << std::fixed << calculateCorrelationSz0(sites,psi,i,j) << " | "; }
                if(type=="1"){ std::cout << std::fixed << calculateCorrelationSz1(sites,psi,i,j) << " | "; }
                if(type=="01"){ std::cout << std::fixed << calculateCorrelationSz01(sites,psi,i,j) << " | "; }
            }
            std::cout << std::endl;
        }
    }

    Sweeps prepareSweepClass()
    {
        auto sweeps = Sweeps(Args::global().getInt("sweeps",4));
        sweeps.maxdim() = Args::global().getInt("maxDim",100);
        sweeps.mindim() = Args::global().getInt("maxDim",1);
        sweeps.cutoff() = Args::global().getReal("cutoff",1E-6);
        sweeps.niter() = 10;
        return sweeps;
    }

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
        std::cout << "  Sz_0: " << calculateSz_0(sites, psi) << std::endl;
        std::cout << "  Sz_1: " << calculateSz_1(sites, psi) << std::endl;
        std::cout << "  Sz_01: " << calculateSz_01(sites, psi) << std::endl;
        std::cout << "  Time: " << (clock()-time0)/(double)CLOCKS_PER_SEC << " [s]" << std::endl;

        printfln("-------------------------------------DMRG---------------------------------");
        time0 = clock();

        dmrg(psi,H,sweeps);
        printfln("-------------------------------------DATA---------------------------------");
        std::cout << "  Energy after DMRG: " << real(innerC(psi,H,psi)) << std::endl;
        std::cout << "  N: " << calculateN(sites, psi) << std::endl;
        std::cout << "  N dublon: " << calculateNd(sites, psi) << std::endl;
        std::cout << "  Sz_0: " << calculateSz_0(sites, psi) << std::endl;
        std::cout << "  Sz_1: " << calculateSz_1(sites, psi) << std::endl;
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
        std::cout << "  Sz_0: " << calculateSz_0(sites, psi) << std::endl;
        std::cout << "  Sz_1: " << calculateSz_1(sites, psi) << std::endl;
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

}



#endif // EXPERIMENT1

