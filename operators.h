#ifndef OPERATORS
#define OPERATORS

#include "itensor/all.h"
#include "tdvp.h"
#include "kondo_heisenberg.h"

using namespace itensor;

double calculateN(const BasicSiteSet<KondoHeisenbergSite> &sites, const MPS &psi)
{
    auto N = AutoMPO(sites);

    for(int i=1; i<=psi.length(); i++){
        N += 1,"n_0",i;
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

#endif // OPERATORS

