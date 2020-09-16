#include "interface.h"
#include "model.h"
#include "tdvp.h"
#include "itensor/util/print_macro.h"
#include "basisextension.h"



int main(int argc, char *argv[])
{
    Experiments("Dmrg") = [](){

        ExpCon.addPoint("Initialization");
        auto [sites,psi,H,sweeps] = prepareExpBasic();
        ExpCon.setSites(sites);ExpCon("E") = H;
        ExpCon.calc(psi,oMode::b,"E","N","Nd","Sz0","Sz1","Szt");

        ExpCon.addPoint("Starting DMRG");
        dmrg(psi,H,sweeps);

        ExpCon.addPoint("Output data");
        ExpCon.calc(psi,oMode::b,"E","N","Nd","Sz0","Sz1","Szt",maxLinkDim(psi));
    };

    Experiments("DmrgWithCorrelations") = [](){

        ExpCon.addPoint("Initialization");
        auto [sites,psi,H,sweeps] = prepareExpBasic();
        ExpCon.setSites(sites);ExpCon("E") = H;
        ExpCon.calc(psi, oMode::a, "E","N","Nd","Sz0","Sz1","Szt");

        ExpCon.addPoint("Starting DMRG");
        dmrg(psi,H,sweeps);

        ExpCon.addPoint("Output data");
        ExpCon.calc(psi,"E","N","Nd","Sz0","Sz1","Szt",maxLinkDim(psi));
        calculateCorrelationMatrixSz(sites,psi, "Sz0");
        calculateCorrelationMatrixSz(sites,psi, "Sz1");
        calculateCorrelationMatrixSz(sites,psi, "SzSz");
        calculateCorrelationMatrixSz(sites,psi, "SmSp");
        calculateCorrelationMatrixSz(sites,psi, "SpSm");
    };

    Experiments("timeEv") = [](){

        ExpCon.addPoint("Initialization");
        auto [sites,psi,H,sweeps] = prepareExpBasic();
        ExpCon.setSites(sites);ExpCon("E") = H;
        ExpCon.calc(psi,oMode::b,"E","N","Nd","Sz0","Sz1","Szt");

        ExpCon.addPoint("Starting TDVP");

        for(double time=0; time<=getD("maxtime")+getD("dtime")+0.001; time+=getD("dtime")){
            ExpCon.calc(psi,oMode::b,"t:",time,"E","N","Nd","Sz0","Sz1","Szt",maxLinkDim(psi),"Sz1_1:L","Sz0_1:L","N1:L");

            if(time<getI("basisExtSteps")*getD("dtime")){
               std::vector<Real> epsilonK = {getD("cutoff"),getD("cutoff"),getD("cutoff")};
               addBasis(psi,H,epsilonK,{"Cutoff",getD("cutoff"),"Method","DensityMatrix","KrylovOrd",4,"DoNormalize",true,"Quiet",true});
            }
            tdvp(psi,H,im*getD("dtime"),sweeps,{"DoNormalize",true,"Quiet",true,"NumCenter",2});
        }
        ExpCon.addPoint("Finish");
    };


    Params.add("thop","double","0.0");
    Params.add("U","double","0.0");
    Params.add("K","double","0.0");
    Params.add("Jh","double","0.0");
    Params.add("Mu","double","0.0");

    Params.add("L","int","4");
    Params.add("PBC","bool","0");

    Params.add("Silent","bool","1");
    Params.add("cutoff","double","1E-6");
    Params.add("sweeps","int","4");
    Params.add("minDim","int","1");
    Params.add("maxDim","int","100");
    Params.add("niter","int","10");
    Params.add("state","string","ud-du");

    Params.add("ConserveN","bool","0");
    Params.add("ConserveSz","bool","0");
    Params.add("maxtime","double","0");
    Params.add("dtime","double","0");

    Params.add("exp","string","1");
    Params.add("basisExtSteps","int","2");

    Params.set(argc,argv);
    prepareObservables();
    Experiments.run();

    return 0;
}
