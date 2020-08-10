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
                 s = Index{ QN( {"Sz=", +1 }, {"Nf=",0} ),1,  //0_u
                            QN( {"Sz=", -1 }, {"Nf=",0} ),1,  //0_d
                            QN( {"Sz=", +2 }, {"Nf=",1} ),1,  //u_u
                            QN( {"Sz=", +0 }, {"Nf=",1} ),2,  //u_d d_u
                            QN( {"Sz=", -2 }, {"Nf=",1} ),1,  //d_d
                            QN( {"Sz=", +1 }, {"Nf=",2} ),1,  //ud_u
                            QN( {"Sz=", -1 }, {"Nf=",2} ),1,  //ud_d
                            Out, ts};
             } else if(conserveSz && (~conserveN)){
                 s = Index{ QN( {"Sz=", +1 } ),1,  //0_u
                            QN( {"Sz=", -1 } ),1,  //0_d
                            QN( {"Sz=", +2 } ),1,  //u_u
                            QN( {"Sz=", +0 } ),2,  //u_d d_u
                            QN( {"Sz=", -2 } ),1,  //d_d
                            QN( {"Sz=", +1 } ),1,  //ud_u
                            QN( {"Sz=", -1 } ),1,  //ud_d
                            Out, ts};
             } else if((~conserveSz) && conserveN){
                 s = Index{ QN( {"Nf=",0} ),2,  //0_u 0_d
                            QN( {"Nf=",1} ),4,  //u_u u_d d_u d_d
                            QN( {"Nf=",2} ),2,  //ud_u ud_d
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
            if(state == "0" || state == "0_U"){ return s(1); }
            else if(state == "1" || state == "0_D"){ return s(2); }
            else if(state == "2" || state == "U_U"){ return s(3); }
            else if(state == "3" || state == "U_D"){ return s(4); }
            else if(state == "4" || state == "D_U"){ return s(5); }
            else if(state == "5" || state == "D_D"){ return s(6); }
            else if(state == "6" || state == "UD_U"){ return s(7); }
            else if(state == "7" || state == "UD_D"){ return s(8); }
            else {
                Error("State " + state + " not recognized");
            }

            return IndexVal{};
        }


        ITensor op(std::string const& opname, Args const& args) const
        {
            auto sP = prime(s);

            auto Op = ITensor{ dag(s),sP };

            if(opname == "c_0,u"){
                Op.set(s(1),sP(3),1);
                Op.set(s(2),sP(4),1);
                Op.set(s(5),sP(7),-1);
                Op.set(s(6),sP(8),-1);
            } else if(opname == "cT_0,u"){
                Op.set(s(3),sP(1),1);
                Op.set(s(4),sP(2),1);
                Op.set(s(7),sP(5),-1);
                Op.set(s(8),sP(6),-1);
            } else if(opname == "c_0,d"){
                Op.set(s(1),sP(5),1);
                Op.set(s(2),sP(6),1);
                Op.set(s(3),sP(7),-1);
                Op.set(s(4),sP(8),-1);
            } else if(opname == "cT_0,d"){
                Op.set(s(5),sP(1),1);
                Op.set(s(6),sP(2),1);
                Op.set(s(7),sP(3),-1);
                Op.set(s(8),sP(4),-1);

            } else if(opname == "n_0,ud"){
                Op.set(s(7),sP(7),1);
                Op.set(s(8),sP(8),1);

            } else if(opname == "n_0"){
                Op.set(s(3),sP(3),1);
                Op.set(s(4),sP(4),1);
                Op.set(s(5),sP(5),1);
                Op.set(s(6),sP(6),1);
                Op.set(s(7),sP(7),2);
                Op.set(s(8),sP(8),2);

            } else if(opname == "s+_1"){
                Op.set(s(1),sP(2),1);
                Op.set(s(3),sP(4),1);
                Op.set(s(5),sP(6),1);
                Op.set(s(7),sP(8),1);
            } else if(opname == "s-_1"){
                Op.set(s(2),sP(1),1);
                Op.set(s(4),sP(3),1);
                Op.set(s(6),sP(5),1);
                Op.set(s(8),sP(7),1);
            } else if(opname == "sz_1"){
                Op.set(s(1),sP(1),0.5);
                Op.set(s(3),sP(3),0.5);
                Op.set(s(5),sP(5),0.5);
                Op.set(s(7),sP(7),0.5);
                Op.set(s(2),sP(2),-0.5);
                Op.set(s(4),sP(4),-0.5);
                Op.set(s(6),sP(6),-0.5);
                Op.set(s(8),sP(8),-0.5);

            } else if(opname == "sz_0"){
                Op.set(s(3),sP(3),0.5);
                Op.set(s(4),sP(4),0.5);
                Op.set(s(5),sP(5),-0.5);
                Op.set(s(6),sP(6),-0.5);

            } else if(opname == "sz_01"){
                Op.set(s(1),sP(1),0.5);
                Op.set(s(2),sP(2),-0.5);
                Op.set(s(3),sP(3),1.0);
                Op.set(s(6),sP(6),-1.0);
                Op.set(s(7),sP(7),0.5);
                Op.set(s(8),sP(8),-0.5);

            } else if(opname == "s_01"){
                Op.set(s(3),sP(3),0.25);
                Op.set(s(6),sP(6),0.25);
                Op.set(s(4),sP(5),0.5);
                Op.set(s(5),sP(4),0.5);
                Op.set(s(4),sP(4),-0.25);  //-0.25
                Op.set(s(5),sP(5),-0.25);  //-0.25

            } else if(opname == "n_01"){
                Op.set(s(1),sP(1), 1);
                Op.set(s(2),sP(2), 1);
                Op.set(s(3),sP(3), 2);
                Op.set(s(4),sP(4), 2);
                Op.set(s(5),sP(5), 2);
                Op.set(s(6),sP(6), 2);
                Op.set(s(7),sP(7), 3);
                Op.set(s(8),sP(8), 3);

            } else if(opname == "nd_01"){
                Op.set(s(7),sP(7), 1);
                Op.set(s(8),sP(8), 1);

            } else {
                Error("Operator \"" + opname + "\" name not recognized");
            }

            return Op;
        }
    };


} //namespace itensor

#endif
