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
        op(std::string const& opname, Args const& args) const
        {
            auto sP = prime(s);

            IndexVal eu(s(1)),
                     euP(sP(1)),
                     uu(s(2)),
                     uuP(sP(2)),
                     du(s(3)),
                     duP(sP(3)),
                     Du(s(4)),
                     DuP(sP(4)),
                    ed(s(5)),
                    edP(sP(5)),
                    ud(s(6)),
                    udP(sP(6)),
                    dd(s(7)),
                    ddP(sP(7)),
                    Dd(s(8)),
                    DdP(sP(8));

            ITensor Op(dag(s),sP);

            if(opname == "Nup"){
                Op.set(uu,uuP,1);
                Op.set(Du,DuP,1);
                Op.set(ud,udP,1);
                Op.set(Dd,DdP,1);
            } else
            if(opname == "Ndn"){
                Op.set(du,duP,1);
                Op.set(Du,DuP,1);
                Op.set(dd,ddP,1);
                Op.set(Dd,DdP,1);
            } else
            if(opname == "Nupdn")
                {
                Op.set(Du,DuP,1);
                Op.set(Dd,DdP,1);
            } else
            if(opname == "Ntot")
                {
                Op.set(uu,uuP,1);
                Op.set(du,duP,1);
                Op.set(Du,DuP,2);

                Op.set(ud,udP,1);
                Op.set(dd,ddP,1);
                Op.set(Dd,DdP,2);
            } else
            if(opname == "Cup")
                {
                Op.set(uu,euP,1);
                Op.set(Du,duP,1);

                Op.set(ud,edP,1);
                Op.set(Dd,ddP,1);
            } else
            if(opname == "Cdagup")
                {
                Op.set(eu,uuP,1);
                Op.set(du,DuP,1);

                Op.set(ed,udP,1);
                Op.set(dd,DdP,1);
            } else
            if(opname == "Cdn")
                {
                Op.set(du,euP,1);
                Op.set(Du,uuP,-1);

                Op.set(dd,edP,1);
                Op.set(Dd,udP,-1);
            } else
            if(opname == "Cdagdn")
                {
                Op.set(eu,duP,1);
                Op.set(uu,DuP,-1);

                Op.set(ed,ddP,1);
                Op.set(ud,DdP,-1);
            } else
            if(opname == "Aup")
                {
                Op.set(uu,euP,1);
                Op.set(Du,duP,1);

                Op.set(ud,edP,1);
                Op.set(Dd,ddP,1);
            } else
            if(opname == "Adagup")
                {
                Op.set(eu,uuP,1);
                Op.set(du,DuP,1);

                Op.set(ed,udP,1);
                Op.set(dd,DdP,1);
            } else
            if(opname == "Adn")
                {
                Op.set(du,euP,1);
                Op.set(Du,uuP,1);

                Op.set(dd,edP,1);
                Op.set(Dd,udP,1);
            } else
            if(opname == "Adagdn")
                {
                Op.set(eu,duP,1);
                Op.set(uu,DuP,1);

                Op.set(ed,ddP,1);
                Op.set(ud,DdP,1);
            } else
            if(opname == "FermiPhase" || opname == "F")
                {
                Op.set(eu,euP,+1);
                Op.set(uu,uuP,-1);
                Op.set(du,duP,-1);
                Op.set(Du,DuP,+1);

                Op.set(ed,edP,+1);
                Op.set(ud,udP,-1);
                Op.set(dd,ddP,-1);
                Op.set(Dd,DdP,+1);
            } else
            if(opname == "Fup")
                {
                Op.set(eu,euP,+1);
                Op.set(uu,uuP,-1);
                Op.set(du,duP,+1);
                Op.set(Du,DuP,-1);

                Op.set(ed,edP,+1);
                Op.set(ud,udP,-1);
                Op.set(dd,ddP,+1);
                Op.set(Dd,DdP,-1);
            } else
            if(opname == "Fdn")
                {
                Op.set(eu,euP,+1);
                Op.set(uu,uuP,+1);
                Op.set(du,duP,-1);
                Op.set(Du,DuP,-1);

                Op.set(ed,edP,+1);
                Op.set(ud,udP,+1);
                Op.set(dd,ddP,-1);
                Op.set(Dd,DdP,-1);
            } else
            if(opname == "Sz0")
                {
                Op.set(uu,uuP,+0.5);
                Op.set(du,duP,-0.5);

                Op.set(ud,udP,+0.5);
                Op.set(dd,ddP,-0.5);
            } else
            if(opname == "Sz1")
            {
                Op.set(eu,euP,+0.5);
                Op.set(uu,uuP,+0.5);
                Op.set(du,duP,+0.5);
                Op.set(Du,DuP,+0.5);

                Op.set(ed,edP,-0.5);
                Op.set(ud,udP,-0.5);
                Op.set(dd,ddP,-0.5);
                Op.set(Dd,DdP,-0.5);
            } else
            if(opname == "S01")
            {
                Op.set(uu,uuP,+0.25);
                Op.set(du,duP,-0.25);
                Op.set(ud,udP,-0.25);
                Op.set(dd,ddP,+0.25);

                Op.set(ud,duP,0.5);
                Op.set(du,udP,0.5);
            } else
            if(opname == "S+1")
            {
                Op.set(ed,euP,1);
                Op.set(ud,uuP,1);
                Op.set(dd,duP,1);
                Op.set(Dd,DuP,1);
            } else
            if(opname == "S-1"){
                Op.set(eu,edP,1);
                Op.set(uu,udP,1);
                Op.set(du,ddP,1);
                Op.set(Du,DdP,1);
            } else
            if(opname == "S+0"){
                Op.set(du,uuP,1);
                Op.set(dd,udP,1);
            } else
            if(opname == "S-0"){
                Op.set(uu,duP,1);
                Op.set(ud,ddP,1);
            } else
            if(opname == "S20"){
                Op.set(uu,uuP,0.75);
                Op.set(du,duP,0.75);

                Op.set(ud,udP,0.75);
                Op.set(dd,ddP,0.75);
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

