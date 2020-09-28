#ifndef CUSTOMELECTRON
#define CUSTOMELECTRON

#pragma once

#include "itensor/mps/siteset.h"

namespace itensor {

class KHSite;
using KH = BasicSiteSet<KHSite>;


class KHSite
    {
    Index s;
    public:

    KHSite(Index I) : s(I) { }

    KHSite(Args const& args = Args::global())
        {
        auto ts = TagSet("Site,Elec");

        if(args.defined("SiteNumber")){
            ts.addTags("n="+str(args.getInt("SiteNumber")));
        }

        if(args.getBool("ConserveN") && args.getBool("ConserveSz")) {
            s = Index(QN({"Sz", 1},{"N",0,-1}),1,
                      QN({"Sz", 2},{"N",1,-1}),1,
                      QN({"Sz", 0},{"N",1,-1}),1,
                      QN({"Sz", 1},{"N",2,-1}),1,
                      QN({"Sz",-1},{"N",0,-1}),1,
                      QN({"Sz", 0},{"N",1,-1}),1,
                      QN({"Sz",-2},{"N",1,-1}),1,
                      QN({"Sz",-1},{"N",2,-1}),1,Out,ts);
            } else if(args.getBool("ConserveN")) {
                s = Index(QN({"N",0,-1}),1, //0u
                          QN({"N",1,-1}),1, //uu
                          QN({"N",1,-1}),1, //du
                          QN({"N",2,-1}),1, //2u
                          QN({"N",0,-1}),1, //0d
                          QN({"N",1,-1}),1, //ud
                          QN({"N",1,-1}),1, //dd
                          QN({"N",2,-1}),1,Out,ts); //2d
            } else if(args.getBool("ConserveSz")) {
                s = Index(QN({"Sz", 1}),1,
                          QN({"Sz", 2}),1,
                          QN({"Sz", 0}),1,
                          QN({"Sz", 1}),1,
                          QN({"Sz",-1}),1,
                          QN({"Sz", 0}),1,
                          QN({"Sz",-2}),1,
                          QN({"Sz",-1}),1,Out,ts);
            } else {
                s = Index(8,ts);
            }
        }

    Index
    index() const { return s; }

    IndexVal
    state(std::string const& state)
    {
        if(     state == "1" || state == "0u"){ return s(1); }
        else if(state == "2" || state == "uu"){ return s(2);}
        else if(state == "3" || state == "du"){ return s(3);}
        else if(state == "4" || state == "2u"){ return s(4); }
        else if(state == "5" || state == "0d"){ return s(5);}
        else if(state == "6" || state == "ud"){ return s(6);}
        else if(state == "7" || state == "dd"){ return s(7);}
        else if(state == "8" || state == "2d"){ return s(8);}
        else{ Error("State " + state + " not recognized"); }

        return IndexVal{};
    }

        ITensor
        op(std::string const& opname,
           Args const& args) const
        {
        auto sP = prime(s);

        IndexVal eu(s(1)),
                 euP(sP(1)),
                 UpU(s(2)),
                 UpPU(sP(2)),
                 DnU(s(3)),
                 DnPU(sP(3)),
                 UDU(s(4)),
                 UDPU(sP(4)),
                EmD(s(5)),
                EmPD(sP(5)),
                UpD(s(6)),
                UpPD(sP(6)),
                DnD(s(7)),
                DnPD(sP(7)),
                UDD(s(8)),
                UDPD(sP(8));

        ITensor Op(dag(s),sP);

        if(opname == "Nup")
            {
            Op.set(UpU,UpPU,1);
            Op.set(UDU,UDPU,1);
            Op.set(UpD,UpPD,1);
            Op.set(UDD,UDPD,1);
            }
        else
        if(opname == "Ndn")
            {
            Op.set(DnU,DnPU,1);
            Op.set(UDU,UDPU,1);
            Op.set(DnD,DnPD,1);
            Op.set(UDD,UDPD,1);
            }
        else
        if(opname == "Nupdn")
            {
            Op.set(UDU,UDPU,1);
            Op.set(UDD,UDPD,1);
            }
        else
        if(opname == "Ntot")
            {
            Op.set(UpU,UpPU,1);
            Op.set(DnU,DnPU,1);
            Op.set(UDU,UDPU,2);

            Op.set(UpD,UpPD,1);
            Op.set(DnD,DnPD,1);
            Op.set(UDD,UDPD,2);
            }
        else
        if(opname == "Cup")
            {
            Op.set(UpU,euP,1);
            Op.set(UDU,DnPU,1);

            Op.set(UpD,EmPD,1);
            Op.set(UDD,DnPD,1);
            }
        else
        if(opname == "Cdagup")
            {
            Op.set(eu,UpPU,1);
            Op.set(DnU,UDPU,1);

            Op.set(EmD,UpPD,1);
            Op.set(DnD,UDPD,1);
            }
        else
        if(opname == "Cdn")
            {
            Op.set(DnU,euP,1);
            Op.set(UDU,UpPU,-1);

            Op.set(DnD,EmPD,1);
            Op.set(UDD,UpPD,-1);
            }
        else
        if(opname == "Cdagdn")
            {
            Op.set(eu,DnPU,1);
            Op.set(UpU,UDPU,-1);

            Op.set(EmD,DnPD,1);
            Op.set(UpD,UDPD,-1);
            }
        else
        if(opname == "Aup")
            {
            Op.set(UpU,euP,1);
            Op.set(UDU,DnPU,1);

            Op.set(UpD,EmPD,1);
            Op.set(UDD,DnPD,1);
            }
        else
        if(opname == "Adagup")
            {
            Op.set(eu,UpPU,1);
            Op.set(DnU,UDPU,1);

            Op.set(EmD,UpPD,1);
            Op.set(DnD,UDPD,1);
            }
        else
        if(opname == "Adn")
            {
            Op.set(DnU,euP,1);
            Op.set(UDU,UpPU,1);

            Op.set(DnD,EmPD,1);
            Op.set(UDD,UpPD,1);
            }
        else
        if(opname == "Adagdn")
            {
            Op.set(eu,DnPU,1);
            Op.set(UpU,UDPU,1);

            Op.set(EmD,DnPD,1);
            Op.set(UpD,UDPD,1);
            }
        else
        if(opname == "FermiPhase" || opname == "F")
            {
            Op.set(eu,euP,+1);
            Op.set(UpU,UpPU,-1);
            Op.set(DnU,DnPU,-1);
            Op.set(UDU,UDPU,+1);

            Op.set(EmD,EmPD,+1);
            Op.set(UpD,UpPD,-1);
            Op.set(DnD,DnPD,-1);
            Op.set(UDD,UDPD,+1);
            }
        else
        if(opname == "Fup")
            {
            Op.set(eu,euP,+1);
            Op.set(UpU,UpPU,-1);
            Op.set(DnU,DnPU,+1);
            Op.set(UDU,UDPU,-1);

            Op.set(EmD,EmPD,+1);
            Op.set(UpD,UpPD,-1);
            Op.set(DnD,DnPD,+1);
            Op.set(UDD,UDPD,-1);
            }
        else
        if(opname == "Fdn")
            {
            Op.set(eu,euP,+1);
            Op.set(UpU,UpPU,+1);
            Op.set(DnU,DnPU,-1);
            Op.set(UDU,UDPU,-1);

            Op.set(EmD,EmPD,+1);
            Op.set(UpD,UpPD,+1);
            Op.set(DnD,DnPD,-1);
            Op.set(UDD,UDPD,-1);
            }
        else
        if(opname == "Sz0")
            {
            Op.set(UpU,UpPU,+0.5);
            Op.set(DnU,DnPU,-0.5);

            Op.set(UpD,UpPD,+0.5);
            Op.set(DnD,DnPD,-0.5);
            }
        else
        if(opname == "Sz1")
        {
            Op.set(eu,euP,+0.5);
            Op.set(UpU,UpPU,+0.5);
            Op.set(DnU,DnPU,+0.5);
            Op.set(UDU,UDPU,+0.5);

            Op.set(EmD,EmPD,-0.5);
            Op.set(UpD,UpPD,-0.5);
            Op.set(DnD,DnPD,-0.5);
            Op.set(UDD,UDPD,-0.5);
        }
        else
        if(opname == "S01")
        {
            Op.set(UpU,UpPU,+0.25);
            Op.set(DnU,DnPU,-0.25);
            Op.set(UpD,UpPD,-0.25);
            Op.set(DnD,DnPD,+0.25);

            Op.set(UpD,DnPU,0.5);
            Op.set(DnU,UpPD,0.5);
        }
        else
        if(opname == "S+1")
        {
            Op.set(EmD,euP,1);
            Op.set(UpD,UpPU,1);
            Op.set(DnD,DnPU,1);
            Op.set(UDD,UDPU,1);
        }
        else
        if(opname == "S-1")
        {
            Op.set(eu,EmPD,1);
            Op.set(UpU,UpPD,1);
            Op.set(DnU,DnPD,1);
            Op.set(UDU,UDPD,1);
        }
        else
        if(opname == "S+0")
            {
            Op.set(DnU,UpPU,1);
            Op.set(DnD,UpPD,1);
            }
        else
        if(opname == "S-0")
            {
            Op.set(UpU,DnPU,1);

            Op.set(UpD,DnPD,1);
            }
        else
        if(opname == "S20")
            {
            //S dot S on-site
            Op.set(UpU,UpPU,0.75);
            Op.set(DnU,DnPU,0.75);

            Op.set(UpD,UpPD,0.75);
            Op.set(DnD,DnPD,0.75);
            }
        else
            {
            Error("Operator \"" + opname + "\" name not recognized");
            }

        return Op;
        }

    KHSite(int n, Args const& args = Args::global())
        {
        *this = KHSite({args,"SiteNumber=",n});
        }

    };

} //namespace itensor

#endif // CUSTOMELECTRON

