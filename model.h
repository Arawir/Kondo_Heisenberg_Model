#ifndef MODEL
#define MODEL

#include "itensor/all.h"
#include "kondo_heisenberg.h"
#include "interface.h"

using namespace itensor;



Sweeps prepareSweepClass()
{
    auto sweeps = Sweeps(Args::global().getInt("sweeps"));
    sweeps.maxdim() = Args::global().getInt("maxDim");
    sweeps.mindim() = Args::global().getInt("minDim");
    sweeps.cutoff() = Args::global().getReal("cutoff");
    sweeps.niter() = Args::global().getReal("niter");
    return sweeps;
}

MPO KHHamiltonian(KH &sites,
                       int L, double thop, double K, double Jh, double Mu, double U)
{
    auto ampo = AutoMPO(sites);
    for(int j=1; j<L; j++){
        ampo += +thop,"Cdagup",j,"Cup",j+1;
       // ampo += +thop,"Cdagup",j+1,"Cup",j;
        ampo += -thop,"Cup",j,"Cdagup",j+1;
        ampo += +thop,"Cdagdn",j,"Cdn",j+1;
       // ampo += +thop,"Cdagdn",j+1,"Cdn",j;
        ampo += -thop,"Cdn",j,"Cdagdn",j+1;

        ampo += K/2,"S+1",j,"S-1",j+1;
        ampo += K/2,"S-1",j,"S+1",j+1;
        ampo += K,"Sz1",j,"Sz1",j+1;
    }

    if(Args::global().getBool("PBC")){
        ampo += +thop,"Cdagup",L,"Cup",1;
       // ampo += +thop,"Cdagup",1,"Cup",L;
        ampo += -thop,"Cup",L,"Cdagup",1;
        ampo += +thop,"Cdagdn",L,"Cdn",1;
       // ampo += +thop,"Cdagdn",1,"Cdn",L;
        ampo += -thop,"Cdn",L,"Cdagdn",1;

        ampo += K/2,"S+1",L,"S-1",1;
        ampo += K/2,"S-1",L,"S+1",1;
        ampo += K,"Sz1",L,"Sz1",1;
    }

    for(int j=1; j<=L; j++){
        ampo += +U,"Nupdn",j;
        ampo += +Jh,"S01",j;
    }

    return toMPO(ampo);
}

void prepareObservables()
{
    Obs("N") = [](const BasicSiteSet<KHSite> &sites){
        auto N = AutoMPO(sites);

        for(int i=1; i<=getI("L"); i++){
            N += 1,"Nupdn",i;
        }

        return toMPO(N);
    };
}

std::tuple<KH,MPS,MPO,Sweeps> prepareExpBasic()
{
    seedRNG(1);
    auto sites = KH( getI("L") );
    auto psi = prepareInitState(sites);
    auto H = KHHamiltonian(sites,getI("L"),getD("thop"),getD("K"),getD("Jh"),getD("Mu"),getD("U"));
    auto sweeps = prepareSweepClass();

    return std::make_tuple( sites,psi,H,sweeps );
}

double calculateN(const BasicSiteSet<KHSite> &sites, const MPS &psi)
{
    auto N = AutoMPO(sites);

    for(int i=1; i<=psi.length(); i++){
        N += 1,"Ntot",i;
    }

    return innerC(psi,toMPO(N),psi).real();
}



double calculateNd(const BasicSiteSet<KHSite> &sites, const MPS &psi)
{
    auto N = AutoMPO(sites);

    for(int i=1; i<=psi.length(); i++){
        N += 1,"Nupdn",i;
    }

    return innerC(psi,toMPO(N),psi).real();
}

double calculateSz0(const BasicSiteSet<KHSite> &sites, const MPS &psi)
{
    auto N = AutoMPO(sites);

    for(int i=1; i<=psi.length(); i++){
        N += 1,"Sz0",i;
    }

    return innerC(psi,toMPO(N),psi).real();
}

double calculateSz1(const BasicSiteSet<KHSite> &sites, const MPS &psi)
{
    auto N = AutoMPO(sites);

    for(int i=1; i<=psi.length(); i++){
        N += 1,"Sz1",i;
    }

    return innerC(psi,toMPO(N),psi).real();
}

double calculateSzt(const BasicSiteSet<KHSite> &sites, const MPS &psi)
{
    auto N = AutoMPO(sites);

    for(int i=1; i<=psi.length(); i++){
        N += 1,"Sz0",i;
        N += 1,"Sz1",i;
    }

    return innerC(psi,toMPO(N),psi).real();
}


double calculateCorrelation(const BasicSiteSet<KHSite> &sites, const MPS &psi, int i, int j, std::string type)
{
    auto Si = AutoMPO(sites);
    auto Sj = AutoMPO(sites);

    if(type=="Sz0"){
        Si += 1,"Sz0",i;
        Sj += 1,"Sz0",j;
    } else if(type=="Sz1"){
        Si += 1,"Sz1",i;
        Sj += 1,"Sz1",j;
    } else if(type=="SzSz"){
        Si += 1,"Sz0",i;
        Si += 1,"Sz1",i;
        Sj += 1,"Sz0",j;
        Sj += 1,"Sz1",j;
    } else if(type=="SpSm"){
        Si += 1,"S+0",i;
        Si += 1,"S+1",i;
        Sj += 1,"S-0",j;
        Sj += 1,"S-1",j;
    } else if(type=="SmSp"){
        Si += 1,"S-0",i;
        Si += 1,"S-1",i;
        Sj += 1,"S+0",j;
        Sj += 1,"S+1",j;
    }

    MPO SiSj = nmultMPO(toMPO(Si),prime(toMPO(Sj)));

    return inner(psi,SiSj,psi);
}

void calculateCorrelationMatrixSz(const BasicSiteSet<KHSite> &sites, const MPS &psi, std::string type)
{
    int L = Args::global().getInt("L");

    if(type=="Sz0") std::cout << "Correlation matrix <Sz0_i Sz0_j>" << std::endl;
    if(type=="Sz1") std::cout << "Correlation matrix <Sz1_i Sz1_j>" << std::endl;
    if(type=="SzSz") std::cout << "Correlation matrix <Sz_i Sz_j>" << std::endl;
    if(type=="SpSm") std::cout << "Correlation matrix <Sp_i Sm_j>" << std::endl;
    if(type=="SmSp") std::cout << "Correlation matrix <Sm_i Sp_j>" << std::endl;
    for(int i=1; i<=L; i++){
        for(int j=1; j<=L; j++){
            std::cout.width(8);
            std::cout.precision(4);
            std::cout << std::fixed << calculateCorrelation(sites,psi,i,j,type) << " ";
        }
        std::cout << std::endl;
    }
}

#endif // MODEL

