#ifndef MODEL
#define MODEL

#include "itensor/all.h"
#include "kondo_heisenberg.h"

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
        ampo += +thop,"Cdagup",j+1,"Cup",j;
        ampo += +thop,"Cdagdn",j,"Cdn",j+1;
        ampo += +thop,"Cdagdn",j+1,"Cdn",j;

        ampo += K/2,"S+1",j,"S-1",j+1;
        ampo += K/2,"S-1",j,"S+1",j+1;
        ampo += K,"Sz1",j,"Sz1",j+1;
    }

    if(Args::global().getBool("PBC")){
        ampo += +thop,"Cdagup",L,"Cup",1;
        ampo += +thop,"Cdagup",1,"Cup",L;
        ampo += +thop,"Cdagdn",L,"Cdn",1;
        ampo += +thop,"Cdagdn",1,"Cdn",L;

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

double calculateN(const BasicSiteSet<KHSite> &sites, const MPS &psi)
{
    auto N = AutoMPO(sites);

    for(int i=1; i<=psi.length(); i++){
        N += 1,"Ntot",i;
    }

    return inner(psi,toMPO(N),psi);
}

double calculateNd(const BasicSiteSet<KHSite> &sites, const MPS &psi)
{
    auto N = AutoMPO(sites);

    for(int i=1; i<=psi.length(); i++){
        N += 1,"Nupdn",i;
    }

    return inner(psi,toMPO(N),psi);
}

double calculateSz0(const BasicSiteSet<KHSite> &sites, const MPS &psi)
{
    auto N = AutoMPO(sites);

    for(int i=1; i<=psi.length(); i++){
        N += 1,"Sz0",i;
    }

    return inner(psi,toMPO(N),psi);
}

double calculateSz1(const BasicSiteSet<KHSite> &sites, const MPS &psi)
{
    auto N = AutoMPO(sites);

    for(int i=1; i<=psi.length(); i++){
        N += 1,"Sz1",i;
    }

    return inner(psi,toMPO(N),psi);
}

double calculateSzt(const BasicSiteSet<KHSite> &sites, const MPS &psi)
{
    auto N = AutoMPO(sites);

    for(int i=1; i<=psi.length(); i++){
        N += 1,"Sz0",i;
        N += 1,"Sz1",i;
    }

    return inner(psi,toMPO(N),psi);
}

double calculateCorrelationSzt(const BasicSiteSet<KHSite> &sites, const MPS &psi, int i, int j)
{
    auto Sz_i = AutoMPO(sites);
    auto Sz_j = AutoMPO(sites);

    Sz_i += 1,"Sz0",i;
    Sz_i += 1,"Sz1",i;
    Sz_j += 1,"Sz0",j;
    Sz_j += 1,"Sz1",j;

    MPO SziSzj = nmultMPO(toMPO(Sz_i),prime(toMPO(Sz_j)));

    return inner(psi,SziSzj,psi);
}

//double calculateCorrelationSz0(const BasicSiteSet<KHSite> &sites, const MPS &psi, int i, int j)
//{
//    auto Sz_i = AutoMPO(sites);
//    auto Sz_j = AutoMPO(sites);

//    Sz_i += 1,"sz0",i;
//    Sz_j += 1,"sz0",j;

//    MPO SziSzj = nmultMPO(toMPO(Sz_i),prime(toMPO(Sz_j)));

//    return inner(psi,SziSzj,psi);
//}

//double calculateCorrelationSz1(const BasicSiteSet<KHSite> &sites, const MPS &psi, int i, int j)
//{
//    auto Sz_i = AutoMPO(sites);
//    auto Sz_j = AutoMPO(sites);

//    Sz_i += 1,"sz1",i;
//    Sz_j += 1,"sz1",j;

//    MPO SziSzj = nmultMPO(toMPO(Sz_i),prime(toMPO(Sz_j)));

//    return inner(psi,SziSzj,psi);
//}

//double calculateMagnetizationSz(const BasicSiteSet<KHSite> &sites, const MPS &psi)
//{
//    auto N = AutoMPO(sites);

//    for(int i=1; i<=psi.length(); i++){
//        N += 1,"sztot",i;
//    }

//    return inner(psi,toMPO(N),psi);
//}

void calculateCorrelationMatrixSz(const BasicSiteSet<KHSite> &sites, const MPS &psi)
{
    int L = Args::global().getInt("L");

    std::cout << "Correlation matrix <Szt_i Szt_j>" << std::endl;
    for(int i=1; i<=L; i++){
        std::cout << "  ";
        std::cout.width(3);
        std::cout << i << " | ";
        for(int j=1; j<=L; j++){
            std::cout.width(8);
            std::cout.precision(4);
            std::cout << std::fixed << calculateCorrelationSzt(sites,psi,i,j) << " | ";
        }
        std::cout << std::endl;
    }
}

#endif // MODEL

