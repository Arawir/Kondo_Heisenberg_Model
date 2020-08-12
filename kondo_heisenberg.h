#ifndef __ITENSOR_KONDO_HEISENBERG_H
#define __ITENSOR_KONDO_HEISENBERG_H
#include "itensor/mps/siteset.h"

namespace itensor {

    class KondoHeisenbergSite;
    using KondoHeisenberg = BasicSiteSet<KondoHeisenbergSite>;

    class KondoHeisenbergSite
    {
    private:
        Index s;

    public:
        KondoHeisenbergSite() { }
        KondoHeisenbergSite(Index I) : s(I) { }

        KondoHeisenbergSite(Args const& args = Args::global())
        {
            auto ts = TagSet("Site,KH");
            if( args.defined("SiteNumber") )
              ts.addTags("n="+str(args.getInt("SiteNumber")));

             auto conserveSz = args.getBool("ConserveSz",false);
             auto conserveN = args.getBool("ConserveN",false);

             if(conserveSz && conserveN){
                 s = Index{ QN( {"Sz=", +1 }, {"N=",0} ),1,  // 0U
                            QN( {"Sz=", -1 }, {"N=",0} ),1,  // 0D
                            QN( {"Sz=", +2 }, {"N=",1} ),1,  // uU
                            QN( {"Sz=", +0 }, {"N=",1} ),2,  // uD dU
                            QN( {"Sz=", -2 }, {"N=",1} ),1,  // dD
                            QN( {"Sz=", +1 }, {"N=",2} ),1,  // 2U
                            QN( {"Sz=", -1 }, {"N=",2} ),1,  // 2D
                            Out, ts};
             } else if(conserveSz && (!conserveN)){
                 s = Index{ QN( {"Sz=", +1 } ),1,  // 0U
                            QN( {"Sz=", -1 } ),1,  // 0D
                            QN( {"Sz=", +2 } ),1,  // uU
                            QN( {"Sz=", +0 } ),2,  // uD dU
                            QN( {"Sz=", -2 } ),1,  // dD
                            QN( {"Sz=", +1 } ),1,  // 2U
                            QN( {"Sz=", -1 } ),1,  // 2D
                            Out, ts};
             } else if((!conserveSz) && conserveN){
                 s = Index{ QN( {"N=",0} ),2,  // 0U 0D
                            QN( {"N=",1} ),4,  // uU uD dU dD
                            QN( {"N=",2} ),2,  // 2U 2D
                            Out, ts};
             } else {
                s = Index{ 8,ts};
             }
        }

        Index index() const
        {
            return s;
        }

        IndexVal state(std::string const& state)
        {
            if(state == "0" || state == "0U"){ return s(1); }
            else if(state == "1" || state == "0D"){ return s(2); }
            else if(state == "2" || state == "uU"){ return s(3); }
            else if(state == "3" || state == "uD"){ return s(4); }
            else if(state == "4" || state == "dU"){ return s(5); }
            else if(state == "5" || state == "dD"){ return s(6); }
            else if(state == "6" || state == "2U"){ return s(7); }
            else if(state == "7" || state == "2D"){ return s(8); }
            else {
                Error("State " + state + " not recognized");
            }

            return IndexVal{};
        }


        ITensor op(std::string const& opname, Args const& args) const
        {
            auto sP = prime(s);

            auto Op = ITensor{ dag(s),sP };

            if(opname == "Cup"){
                Op.set(s(3),sP(1),1);
                Op.set(s(4),sP(2),1);
                Op.set(s(7),sP(5),1);
                Op.set(s(8),sP(6),1);
            } else if(opname == "Cdagup"){
                Op.set(s(1),sP(3),1);
                Op.set(s(2),sP(4),1);
                Op.set(s(5),sP(7),1);
                Op.set(s(6),sP(8),1);
            } else if(opname == "Cdn"){
                Op.set(s(5),sP(1),1);
                Op.set(s(6),sP(2),1);
                Op.set(s(7),sP(3),-1);
                Op.set(s(8),sP(4),-1);
            } else if(opname == "Cdagdn"){
                Op.set(s(1),sP(5),1);
                Op.set(s(2),sP(6),1);
                Op.set(s(3),sP(7),-1);
                Op.set(s(4),sP(8),-1);

            } else if(opname == "Aup"){
                Op.set(s(3),sP(1),1);
                Op.set(s(4),sP(2),1);
                Op.set(s(7),sP(5),1);
                Op.set(s(8),sP(6),1);
            } else if(opname == "Adagup"){
                Op.set(s(1),sP(3),1);
                Op.set(s(2),sP(4),1);
                Op.set(s(5),sP(7),1);
                Op.set(s(6),sP(8),1);
            } else if(opname == "Adn"){
                Op.set(s(5),sP(1),1);
                Op.set(s(6),sP(2),1);
                Op.set(s(7),sP(3),1);
                Op.set(s(8),sP(4),1);
            } else if(opname == "Adagdn"){
                Op.set(s(1),sP(5),1);
                Op.set(s(2),sP(6),1);
                Op.set(s(3),sP(7),1);
                Op.set(s(4),sP(8),1);

            } else if(opname == "Fup"){
                Op.set(s(1),sP(1),1);
                Op.set(s(2),sP(2),1);
                Op.set(s(3),sP(3),-1);
                Op.set(s(4),sP(4),-1);
                Op.set(s(5),sP(5),1);
                Op.set(s(6),sP(6),1);
                Op.set(s(7),sP(7),-1);
                Op.set(s(8),sP(8),-1);
            } else if(opname == "Fdn"){
                Op.set(s(1),sP(1),1);
                Op.set(s(2),sP(2),1);
                Op.set(s(3),sP(3),1);
                Op.set(s(4),sP(4),1);
                Op.set(s(5),sP(5),-1);
                Op.set(s(6),sP(6),-1);
                Op.set(s(7),sP(7),-1);
                Op.set(s(8),sP(8),-1);
            } else if(opname == "FermiPhase" || opname == "F" ){
                Op.set(s(1),sP(1),1);
                Op.set(s(2),sP(2),1);
                Op.set(s(3),sP(3),-1);
                Op.set(s(4),sP(4),-1);
                Op.set(s(5),sP(5),-1);
                Op.set(s(6),sP(6),-1);
                Op.set(s(7),sP(7),1);
                Op.set(s(8),sP(8),1);

            } else if(opname == "nD"){
                Op.set(s(7),sP(7),1);
                Op.set(s(8),sP(8),1);

            } else if(opname == "n"){
                Op.set(s(3),sP(3),1);
                Op.set(s(4),sP(4),1);
                Op.set(s(5),sP(5),1);
                Op.set(s(6),sP(6),1);
                Op.set(s(7),sP(7),2);
                Op.set(s(8),sP(8),2);

            } else if(opname == "splus1"){
                Op.set(s(2),sP(1),1);
                Op.set(s(4),sP(3),1);
                Op.set(s(6),sP(5),1);
                Op.set(s(8),sP(7),1);
            } else if(opname == "sminus1"){
                Op.set(s(1),sP(2),1);
                Op.set(s(3),sP(4),1);
                Op.set(s(5),sP(6),1);
                Op.set(s(7),sP(8),1);
            } else if(opname == "sz1"){
                Op.set(s(1),sP(1),0.5);
                Op.set(s(3),sP(3),0.5);
                Op.set(s(5),sP(5),0.5);
                Op.set(s(7),sP(7),0.5);
                Op.set(s(2),sP(2),-0.5);
                Op.set(s(4),sP(4),-0.5);
                Op.set(s(6),sP(6),-0.5);
                Op.set(s(8),sP(8),-0.5);

            } else if(opname == "sz0"){
                Op.set(s(3),sP(3),0.5);
                Op.set(s(4),sP(4),0.5);
                Op.set(s(5),sP(5),-0.5);
                Op.set(s(6),sP(6),-0.5);

            } else if(opname == "sztot"){
                Op.set(s(1),sP(1),0.5);
                Op.set(s(2),sP(2),-0.5);
                Op.set(s(3),sP(3),1.0);
                Op.set(s(6),sP(6),-1.0);
                Op.set(s(7),sP(7),0.5);
                Op.set(s(8),sP(8),-0.5);

            } else if(opname == "shund"){
                Op.set(s(3),sP(3),0.25);
                Op.set(s(6),sP(6),0.25);
                Op.set(s(4),sP(5),0.5);
                Op.set(s(5),sP(4),0.5);
                Op.set(s(4),sP(4),-0.25);
                Op.set(s(5),sP(5),-0.25);

            } else {
                Error("Operator \"" + opname + "\" name not recognized");
            }

            return Op;
        }
    };


} //namespace itensor

#endif
