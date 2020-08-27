#ifndef MODEL
#define MODEL

#include "itensor/all.h"

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

MPO KHHamiltonian(KondoHeisenberg &sites,
                       int L, double thop, double K, double Jh, double Mu, double U)
{
    auto ampo = AutoMPO(sites);
    for(int j=1; j<L; j++){
        ampo += +thop,"Cdagup",j,"Cup",j+1;
        ampo += +thop,"Cup",j,"Cdagup",j+1;
        ampo += +thop,"Cdagdn",j,"Cdn",j+1;
        ampo += +thop,"Cdn",j,"Cdagdn",j+1;

        ampo += +K/2,"splus1",j,"sminus1",j+1;
        ampo += +K/2,"sminus1",j,"splus1",j+1;
        ampo += +K,"sz1",j,"sz1",j+1;
    }

    if(Args::global().getBool("PBC")){
        ampo += +thop,"Cdagup",L,"Cup",1;
        ampo += +thop,"Cup",L,"Cdagup",1;
        ampo += +thop,"Cdagdn",L,"Cdn",1;
        ampo += +thop,"Cdn",L,"Cdagdn",1;

        ampo += +K/2,"splus1",L,"sminus1",1;
        ampo += +K/2,"sminus1",L,"splus1",1;
        ampo += +K,"sz1",L,"sz1",1;
    }

    for(int j=1; j<=L; j++){
        ampo += +U,"nD",j;
        ampo += -2*Jh,"shund",j;
        ampo += +Mu,"n",j;
    }

    return toMPO(ampo);
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
        N += 1,"nD",i;
    }

    return inner(psi,toMPO(N),psi);
}

double calculateSz0(const BasicSiteSet<KondoHeisenbergSite> &sites, const MPS &psi)
{
    auto N = AutoMPO(sites);

    for(int i=1; i<=psi.length(); i++){
        N += 1,"sz0",i;
    }

    return inner(psi,toMPO(N),psi);
}

double calculateSz1(const BasicSiteSet<KondoHeisenbergSite> &sites, const MPS &psi)
{
    auto N = AutoMPO(sites);

    for(int i=1; i<=psi.length(); i++){
        N += 1,"sz1",i;
    }

    return inner(psi,toMPO(N),psi);
}

double calculateSzt(const BasicSiteSet<KondoHeisenbergSite> &sites, const MPS &psi)
{
    auto N = AutoMPO(sites);

    for(int i=1; i<=psi.length(); i++){
        N += 1,"sztot",i;
    }

    return inner(psi,toMPO(N),psi);
}

double calculateCorrelationSzt(const BasicSiteSet<KondoHeisenbergSite> &sites, const MPS &psi, int i, int j)
{
    auto Sz_i = AutoMPO(sites);
    auto Sz_j = AutoMPO(sites);

    Sz_i += 1,"sztot",i;
    Sz_j += 1,"sztot",j;

    MPO SziSzj = nmultMPO(toMPO(Sz_i),prime(toMPO(Sz_j)));

    return inner(psi,SziSzj,psi);
}

double calculateCorrelationSz0(const BasicSiteSet<KondoHeisenbergSite> &sites, const MPS &psi, int i, int j)
{
    auto Sz_i = AutoMPO(sites);
    auto Sz_j = AutoMPO(sites);

    Sz_i += 1,"sz0",i;
    Sz_j += 1,"sz0",j;

    MPO SziSzj = nmultMPO(toMPO(Sz_i),prime(toMPO(Sz_j)));

    return inner(psi,SziSzj,psi);
}

double calculateCorrelationSz1(const BasicSiteSet<KondoHeisenbergSite> &sites, const MPS &psi, int i, int j)
{
    auto Sz_i = AutoMPO(sites);
    auto Sz_j = AutoMPO(sites);

    Sz_i += 1,"sz1",i;
    Sz_j += 1,"sz1",j;

    MPO SziSzj = nmultMPO(toMPO(Sz_i),prime(toMPO(Sz_j)));

    return inner(psi,SziSzj,psi);
}

double calculateMagnetizationSz(const BasicSiteSet<KondoHeisenbergSite> &sites, const MPS &psi)
{
    auto N = AutoMPO(sites);

    for(int i=1; i<=psi.length(); i++){
        N += 1,"sztot",i;
    }

    return inner(psi,toMPO(N),psi);
}

void calculateCorrelationMatrixSz(const BasicSiteSet<KondoHeisenbergSite> &sites, const MPS &psi, std::string type)
{
    int L = Args::global().getInt("L");

    if(type=="0"){ std::cout << "Correlation matrix <Sz0_i Sz0_j>" << std::endl; }
    if(type=="1"){ std::cout << "Correlation matrix <Sz1_i Sz1_j>" << std::endl; }
    if(type=="t"){ std::cout << "Correlation matrix <Szt_i Szt_j>" << std::endl; }
    for(int i=1; i<=L; i++){
        std::cout << "  ";
        std::cout.width(3);
        std::cout << i << " | ";
        for(int j=1; j<=L; j++){
            std::cout.width(8);
            std::cout.precision(4);
            if(type=="0"){ std::cout << std::fixed << calculateCorrelationSz0(sites,psi,i,j) << " | "; }
            if(type=="1"){ std::cout << std::fixed << calculateCorrelationSz1(sites,psi,i,j) << " | "; }
            if(type=="t"){ std::cout << std::fixed << calculateCorrelationSzt(sites,psi,i,j) << " | "; }
        }
        std::cout << std::endl;
    }
}

#endif // MODEL

