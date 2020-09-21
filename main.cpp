#include "interface.h"
#include "model.h"
#include "tdvp.h"
#include "itensor/util/print_macro.h"
#include "basisextension.h"

void tdvpStepWithBasisExtensionIfNeeded(MPS &psi, MPO &H, double dTime, Sweeps &sweeps)
{
    if(maxLinkDim(psi)<10){
       std::vector<Real> epsilonK = {getD("cutoff"),getD("cutoff"),getD("cutoff")};
       addBasis(psi,H,epsilonK,{"Cutoff",getD("cutoff"),"Method","DensityMatrix","KrylovOrd",4,"DoNormalize",true,"Quiet",true});
    }
    tdvp(psi,H,im*dTime,sweeps,{"DoNormalize",true,"Quiet",true,"NumCenter",2});
}


int main(int argc, char *argv[])
{
    Experiments("Dmrg") = [](){
        auto [sites,psi,H,sweeps] = prepareExpBasic();
        ExpCon.setSites(sites);ExpCon("E") = H;
        ExpCon.calc(psi,oMode::b,"E","N","Nd","Sz0","Sz1","Szt");

        ExpCon.addPoint("Starting DMRG");
        dmrg(psi,H,sweeps);

        ExpCon.addPoint("Output data");
        ExpCon.calc(psi,oMode::b,"E","N","Nd","Sz0","Sz1","Szt","dim");
    };

    Experiments("DmrgWithCorrelations") = [](){
        auto [sites,psi,H,sweeps] = prepareExpBasic();
        ExpCon.setSites(sites);ExpCon("E") = H;
        ExpCon.calc(psi, oMode::a, "E","N","Nd","Sz0","Sz1","Szt");

        ExpCon.addPoint("Starting DMRG");
        dmrg(psi,H,sweeps);

        ExpCon.addPoint("Output data");
        ExpCon.calc(psi,"E","N","Nd","Sz0","Sz1","Szt","dim");
        calculateCorrelationMatrixSz(sites,psi, "Sz0");
        calculateCorrelationMatrixSz(sites,psi, "Sz1");
        calculateCorrelationMatrixSz(sites,psi, "SzSz");
        calculateCorrelationMatrixSz(sites,psi, "SmSp");
        calculateCorrelationMatrixSz(sites,psi, "SpSm");
    };


    Experiments("timeEv") = [](){
        auto [sites,psi,H,sweeps] = prepareExpBasic();
        ExpCon.setSites(sites);ExpCon("E") = H;

        ExpCon.addPoint("Starting TDVP");

        for(double time=0; time<=getD("maxtime")+getD("dtime")+0.001; time+=getD("dtime")){
            ExpCon.calc(psi,oMode::b,"t:",time,"rtime","mem","E","N","Nd","Sz0","Sz1","Szt","dim","Sz1_1:L","Sz0_1:L","N1:L");
            tdvpStepWithBasisExtensionIfNeeded(psi,H,getD("dtime"),sweeps);
        }
    };


    Params.add("thop","double","0.5");
    Params.add("U","double","2.1");
    Params.add("Jh","double","-1.05");
    Params.add("K","double","0.042857143");

    Params.add("Mu","double","0.0");

    Params.add("L","int","24");
    Params.add("PBC","bool","0");

    Params.add("Silent","bool","1");
    Params.add("cutoff","double","1E-8");
    Params.add("sweeps","int","4");
    Params.add("minDim","int","1");
    Params.add("maxDim","int","100");
    Params.add("niter","int","30");
    Params.add("state","string","L/2*uu-L/2*dd");

    Params.add("ConserveN","bool","1");
    Params.add("ConserveSz","bool","0");
    Params.add("maxtime","double","48");
    Params.add("dtime","double","0");

    Params.add("exp","string","timeEv");
    Params.add("basisExtSteps","int","2");

    Params.add("PBSenable","bool","0");
    Params.add("PBSjobid","int","0");

    Params.set(argc,argv);
    prepareObservables();
    Experiments.run();

    return 0;
}
