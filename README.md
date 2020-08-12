# Kondo_Heisenberg_Model

## Kompilacja

TODO 

## Uruchomienie

Aby wybrać współczynniki opisujące eksperyment należy uruchomić program z odpowiednimi flagami, na przykład:

$ ./Kondo_Heisenberg U=1E-9 K=1.4 Jh=0.002 t00=1 L=4 state=0-7-0-3-6-6 Mu=2

##Dostępne flagi (nazwy w cudzysłowach):
        if( paramName(argv[i]) == "L" ){ setParam("int", argv[i]); }
        else if( paramName(argv[i]) == "PBC" ){ setParam("bool", argv[i]); }
        else if( paramName(argv[i]) == "t00" ){ setParam("double", argv[i]); }
        else if( paramName(argv[i]) == "U" ){ setParam("double", argv[i]); }
        else if( paramName(argv[i]) == "K" ){ setParam("double", argv[i]); }
        else if( paramName(argv[i]) == "Jh" ){ setParam("double", argv[i]); }
        else if( paramName(argv[i]) == "Mu" ){ setParam("double", argv[i]); }
        else if( paramName(argv[i]) == "t" ){ setParam("double", argv[i]); }
        else if( paramName(argv[i]) == "exp" ){ setParam("int", argv[i]); }
        else if( paramName(argv[i]) == "conSz" ){ setParam("bool", argv[i]); }
        else if( paramName(argv[i]) == "conN" ){ setParam("bool", argv[i]); }
        else if( paramName(argv[i]) == "state" ){ setParam("string", argv[i]); }
        else if( paramName(argv[i]) == "maxDim" ){ setParam("int", argv[i]); }
        else if( paramName(argv[i]) == "minDim" ){ setParam("int", argv[i]); }
        else if( paramName(argv[i]) == "sweeps" ){ setParam("int", argv[i]); }
        else if( paramName(argv[i]) == "cutoff" ){ setParam("double", argv[i]); }
        else if( paramName(argv[i]) == "Silent" ){ setParam("bool", argv[i]); }

##Lista potencjalnie dostępnych stanów węzła. Można wstawić je do argumentu state=S-S-S-S-S...-S w miejsce "S". Można wstawiać zarówno numery (0-7) oraz postać opisową (0u, ..., Dd).
0 = 0u
1 = 0d
2 = uu
3 = ud
4 = du
5 = dd
6 = Du
7 = Dd



