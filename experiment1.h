#ifndef EXPERIMENT1
#define EXPERIMENT1

#include <iostream>
#include <iomanip>
#include <ctime>
#include <cstdlib>
#include <fstream>
#include <vector>
#include "operators.h"


typedef std::complex<double> cpx;
#define im std::complex<double>{0.0,1.0}


namespace Experiment1
{
    double t00 = 0.0;
    double U   = 0.0;
    double K   = 0.0;
    double Jh  = 0.0;
    double Mu  = 0.0;

    bool PBC = false;
    int L = 4;
    auto t = 0.1;

    bool conserveSz = false;
    bool conserveN = false;

    std::vector<std::string> initState{};

    void parseInitState(const char* str)
    {
        int i=0;
        initState.push_back("");
        while(str[i]!='\0'){
            if(str[i]=='-'){
                initState.push_back("");
            } else {
                initState.back().append(std::string{str[i]});
            }
            i++;
        }
    }


    void readExperimentParameters(int argc, char *argv[])
    {
        for(uint i=1; i<argc; i++){
            if((argv[i][0]=='L')&&(argv[i][1]=='=')){
                L = atoi(argv[i]+2);
            } else if((argv[i][0]=='P')&&(argv[i][1]=='B')&&(argv[i][2]=='C')&&(argv[i][3]=='=')){
                PBC = (bool)atoi(argv[i]+4);
            } else if((argv[i][0]=='t')&&(argv[i][1]=='0')&&(argv[i][2]=='0')&&(argv[i][3]=='=')){
                t00 = atof(argv[i]+4);
                Args::global().add("t00",t00);
            } else if((argv[i][0]=='U')&&(argv[i][1]=='=')){
                U = atof(argv[i]+2);
                Args::global().add("U",U);
            } else if((argv[i][0]=='K')&&(argv[i][1]=='=')){
                K = atof(argv[i]+2);
                Args::global().add("K",K);
            } else if((argv[i][0]=='J')&&(argv[i][1]=='h')&&(argv[i][2]=='=')){
                Jh = atof(argv[i]+3);
                Args::global().add("Jh",Jh);
            } else if((argv[i][0]=='M')&&(argv[i][1]=='u')&&(argv[i][2]=='=')){
                Mu = atof(argv[i]+3);
            } else if((argv[i][0]=='t')&&(argv[i][1]=='=')){
                t = atof(argv[i]+2);
            } else if((argv[i][0]=='e')&&(argv[i][1]=='x')&&(argv[i][2]=='p')){
                Args::global().add(argv[i], true);
            } else if((argv[i][0]=='c')&&(argv[i][1]=='o')&&(argv[i][2]=='n')&&(argv[i][3]=='S')&&(argv[i][4]=='z')&&(argv[i][5]=='=')){
                conserveSz = (bool)atoi(argv[i]+6);
                Args::global().add("ConserveSz",conserveSz);
            } else if((argv[i][0]=='c')&&(argv[i][1]=='o')&&(argv[i][2]=='n')&&(argv[i][3]=='N')&&(argv[i][4]=='=')){
                conserveN = (bool)atoi(argv[i]+5);
                Args::global().add("ConserveN",conserveN);
            } else if((argv[i][0]=='s')&&(argv[i][1]=='t')&&(argv[i][2]=='a')&&(argv[i][3]=='t')&&(argv[i][4]=='e')&&(argv[i][5]=='=')){
                parseInitState(argv[i]+6);
            } else {
                std::cerr << "WARNING: unknown argument (" << argv[i] << ")" << std::endl;
            }
        }

        if(initState.size()==0){
            for(int i=0; i<L; i+=2){
               initState.push_back("U_D");
               initState.push_back("D_U");
            }
        } else if(initState.size()<L){
            for(int i=initState.size()-1; i<L; i++){
                std::cerr << "WARNING: wrong init state" << std::endl;
                initState.push_back("U_D");
            }
        }
    }

    std::string convertInitStateString(std::string initState)
    {
        if(initState.back()=='\0'){ initState.pop_back(); }
        if(initState=="0") return "0_U"; //1
        if(initState=="1") return "0_D"; //2
        if(initState=="2") return "U_U"; //3
        if(initState=="3") return "U_D"; //4
        if(initState=="4") return "D_U"; //5
        if(initState=="5") return "D_D"; //6
        if(initState=="6") return "UD_U";//7
        if(initState=="7") return "UD_D";//8
        return initState;
    }

    void writeExperimentParameters()
    {
        std::cout << "Experiment parameters: " << std::endl;
        std::cout << "  L = " << L << std::endl;
        std::cout << "  t00 = " << t00 << std::endl;
        std::cout << "  U = " << U << std::endl;
        std::cout << "  K = " << K << std::endl;
        std::cout << "  Jh = " << Jh << std::endl;
        std::cout << "  Mu = " << Mu << std::endl;
        std::cout << "  PBC = " << PBC << std::endl;
        std::cout << "  Conserve Sz = " << conserveSz << std::endl;
        std::cout << "  Conserve N = " << conserveN << std::endl;

        std::cout << "  Init state = ";
        for(int i=0; i<L-1; i++){
           std::cout << convertInitStateString(initState[i]) << "--";
        }
        std::cout << convertInitStateString(initState[L-1]) << std::endl;

    }

    MPS prepareInitState(BasicSiteSet<KondoHeisenbergSite> &sites)
    {
        auto state = InitState(sites);
        for(int i=1; i<=L; i++){
           state.set(i,initState[i-1]);
        }

        return MPS(state);
    }

    MPO prepareHamiltonian(BasicSiteSet<KondoHeisenbergSite> &sites)
    {
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

        if(PBC){
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
        auto sweeps = Sweeps(4);
        sweeps.maxdim() = 100;
        sweeps.cutoff() = 1E-6;
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
        auto sites = KondoHeisenberg(L);
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

        dmrg(psi,H,sweeps,{"Silent",false});
        //dmrg(psi,H,sweeps);
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
        printfln("-------------------------------------PREAPRE OBJECTS----------------------");
        clock_t  time0 = clock();

        seedRNG(1);
        auto sites = KondoHeisenberg(L);
        auto psi = prepareInitState(sites);
        auto H = prepareHamiltonian(sites);
        auto sweeps = prepareSweepClass();

        std::cout << "Energy of initial state: " << real(innerC(psi,H,psi)) << std::endl;
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
        auto sites = KondoHeisenberg(L);
        auto psi = prepareInitState(sites);
        auto H = prepareHamiltonian(sites);
        auto sweeps = prepareSweepClass();

        std::cout << "  Energy of initial state: " << real(innerC(psi,H,psi)) << std::endl;

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

//    void exp4()
//    {
//        printfln("-------------------------------------EXP4---------------------------------");
//        std::cout << "  iDMRG" << std::endl;

//        printfln("-------------------------------------PREAPRE OBJECTS----------------------");
//        clock_t  time0 = clock();

//        int Nuc = 2;
//        int N = 2*Nuc;

//        auto sites = KondoHeisenberg(L);
//        auto H = Kondo_Heisenberg(sites,{"Infinite=",true});
//        auto sweeps = prepareSweepClass();

//        auto state = InitState(sites);
//        for(int i = 1; i <= N; ++i)
//            {
//            if(i%2 == 1)
//                state.set(i,"U_D");
//            else
//                state.set(i,"D_U");
//            }
//        auto psi = MPS(state);

//        std::cout << "  Energy of initial state: " << real(innerC(psi,H,psi)) << std::endl;

//        auto res = idmrg(psi,H,sweeps,{"OutputLevel",1});
//        printfln("\nGround state energy after DMRG / site = %.20f",res.energy/N);
//    }



    void run()
    {
        if( Args::global().getBool("exp1",false)){
            exp1();
        } else if( Args::global().getBool("exp2",false)){
            exp2();
        } else if( Args::global().getBool("exp3",false)){
            exp3();
        } else if( Args::global().getBool("exp4",false)){
            //exp4();
        } else {
            std::cerr << "ERROR: epxeriment was not selected!" << std::endl;
        }
    }

}



#endif // EXPERIMENT1

