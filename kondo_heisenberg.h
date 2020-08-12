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
                 s = Index{ QN( {"Sz=", +1 }, {"N=",0} ),1,  //0u
                            QN( {"Sz=", -1 }, {"N=",0} ),1,  //0d
                            QN( {"Sz=", +2 }, {"N=",1} ),1,  //uu
                            QN( {"Sz=", +0 }, {"N=",1} ),2,  //ud du
                            QN( {"Sz=", -2 }, {"N=",1} ),1,  //dd
                            QN( {"Sz=", +1 }, {"N=",2} ),1,  //Du
                            QN( {"Sz=", -1 }, {"N=",2} ),1,  //Dd
                            Out, ts};
             } else if(conserveSz && (!conserveN)){
                 s = Index{ QN( {"Sz=", +1 } ),1,  //0u
                            QN( {"Sz=", -1 } ),1,  //0d
                            QN( {"Sz=", +2 } ),1,  //uu
                            QN( {"Sz=", +0 } ),2,  //ud du
                            QN( {"Sz=", -2 } ),1,  //dd
                            QN( {"Sz=", +1 } ),1,  //Du
                            QN( {"Sz=", -1 } ),1,  //Dd
                            Out, ts};
             } else if((!conserveSz) && conserveN){
                 s = Index{ QN( {"N=",0} ),2,  //0u 0d
                            QN( {"N=",1} ),4,  //uu ud du dd
                            QN( {"N=",2} ),2,  //Du Dd
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
            if(state == "0" || state == "0u"){ return s(1); }
            else if(state == "1" || state == "0d"){ return s(2); }
            else if(state == "2" || state == "uu"){ return s(3); }
            else if(state == "3" || state == "ud"){ return s(4); }
            else if(state == "4" || state == "du"){ return s(5); }
            else if(state == "5" || state == "dd"){ return s(6); }
            else if(state == "6" || state == "Du"){ return s(7); }
            else if(state == "7" || state == "Dd"){ return s(8); }
            else {
                Error("State " + state + " not recognized");
            }

            return IndexVal{};
        }


        ITensor op(std::string const& opname, Args const& args) const
        {
            auto sP = prime(s);

            auto Op = ITensor{ dag(s),sP };

            if(opname == "cup"){
                Op.set(s(1),sP(3),1);
                Op.set(s(2),sP(4),1);
                Op.set(s(5),sP(7),-1);
                Op.set(s(6),sP(8),-1);
            } else if(opname == "cdagup"){
                Op.set(s(3),sP(1),1);
                Op.set(s(4),sP(2),1);
                Op.set(s(7),sP(5),-1);
                Op.set(s(8),sP(6),-1);
            } else if(opname == "cdn"){
                Op.set(s(1),sP(5),1);
                Op.set(s(2),sP(6),1);
                Op.set(s(3),sP(7),-1);
                Op.set(s(4),sP(8),-1);
            } else if(opname == "cdagdn"){
                Op.set(s(5),sP(1),1);
                Op.set(s(6),sP(2),1);
                Op.set(s(7),sP(3),-1);
                Op.set(s(8),sP(4),-1);

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
                Op.set(s(1),sP(2),1);
                Op.set(s(3),sP(4),1);
                Op.set(s(5),sP(6),1);
                Op.set(s(7),sP(8),1);
            } else if(opname == "sminus1"){
                Op.set(s(2),sP(1),1);
                Op.set(s(4),sP(3),1);
                Op.set(s(6),sP(5),1);
                Op.set(s(8),sP(7),1);
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

            } else if(opname == "szhund"){
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
