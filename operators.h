#ifndef OPERATORS
#define OPERATORS

#include "itensor/all.h"
#include "tdvp.h"
#include "kondo_heisenberg.h"

using namespace itensor;

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
        ampo += U,"nD",j;
        ampo += -2.0*Jh,"s_01",j;
        ampo += Mu,"n",j;
    }

    return toMPO(ampo);
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


double calculateN(const BasicSiteSet<KondoHeisenbergSite> &sites, const MPS &psi)
{
    auto N = AutoMPO(sites);

    for(int i=1; i<=psi.length(); i++){
        N += 1,"n",i;
    }

    return inner(psi,toMPO(N),psi);
}

double calculateNd(const BasicSiteSet<KondoHeisenbergSite> &sites, const MPS &psi)
{
    auto N = AutoMPO(sites);

    for(int i=1; i<=psi.length(); i++){
        N += 1,"nd_01",i;
    }

    return inner(psi,toMPO(N),psi);
}

double calculateSz_0(const BasicSiteSet<KondoHeisenbergSite> &sites, const MPS &psi)
{
    auto N = AutoMPO(sites);

    for(int i=1; i<=psi.length(); i++){
        N += 1,"sz_0",i;
    }

    return inner(psi,toMPO(N),psi);
}

double calculateSz_1(const BasicSiteSet<KondoHeisenbergSite> &sites, const MPS &psi)
{
    auto N = AutoMPO(sites);

    for(int i=1; i<=psi.length(); i++){
        N += 1,"sz_1",i;
    }

    return inner(psi,toMPO(N),psi);
}

double calculateSz_01(const BasicSiteSet<KondoHeisenbergSite> &sites, const MPS &psi)
{
    auto N = AutoMPO(sites);

    for(int i=1; i<=psi.length(); i++){
        N += 1,"sz_01",i;
    }

    return inner(psi,toMPO(N),psi);
}

double calculateCorrelationSz01(const BasicSiteSet<KondoHeisenbergSite> &sites, const MPS &psi, int i, int j)
{
    auto Sz_i = AutoMPO(sites);
    auto Sz_j = AutoMPO(sites);

    Sz_i += 1,"sz_01",i;
    Sz_j += 1,"sz_01",j;

    MPO SziSzj = nmultMPO(toMPO(Sz_i),prime(toMPO(Sz_j)));

    return inner(psi,SziSzj,psi);
}

double calculateCorrelationSz0(const BasicSiteSet<KondoHeisenbergSite> &sites, const MPS &psi, int i, int j)
{
    auto Sz_i = AutoMPO(sites);
    auto Sz_j = AutoMPO(sites);

    Sz_i += 1,"sz_0",i;
    Sz_j += 1,"sz_0",j;

    MPO SziSzj = nmultMPO(toMPO(Sz_i),prime(toMPO(Sz_j)));

    return inner(psi,SziSzj,psi);
}

double calculateCorrelationSz1(const BasicSiteSet<KondoHeisenbergSite> &sites, const MPS &psi, int i, int j)
{
    auto Sz_i = AutoMPO(sites);
    auto Sz_j = AutoMPO(sites);

    Sz_i += 1,"sz_1",i;
    Sz_j += 1,"sz_1",j;

    MPO SziSzj = nmultMPO(toMPO(Sz_i),prime(toMPO(Sz_j)));

    return inner(psi,SziSzj,psi);
}

double calculateMagnetizationSz(const BasicSiteSet<KondoHeisenbergSite> &sites, const MPS &psi)
{
    auto N = AutoMPO(sites);

    for(int i=1; i<=psi.length(); i++){
        N += 1,"sz_01",i;
    }

    return inner(psi,toMPO(N),psi);
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

#endif // OPERATORS

