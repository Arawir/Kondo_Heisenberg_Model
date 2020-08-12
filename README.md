# Kondo_Heisenberg_Model

## Kompilacja

TODO 

## Uruchomienie

Aby wybrać współczynniki opisujące eksperyment należy uruchomić program z odpowiednimi flagami, na przykład:

$ ./Kondo_Heisenberg U=1E-9 K=1.4 Jh=0.002 t00=1 L=4 state=0-7-0-3-6-6 Mu=2

##Dostępne flagi (nazwy w cudzysłowach):
"L"
"PBC"
"t00"
"U"
"K"
"Jh"
"Mu"
"t"
"exp"
"conSz"
"conN"
"state"
"maxDim"
"minDim"
"sweeps"
"cutoff" 
"Silent" 

##Lista potencjalnie dostępnych stanów węzła. Można wstawić je do argumentu state=S-S-S-S-S...-S w miejsce "S". Można wstawiać zarówno numery (0-7) oraz postać opisową (0u, ..., Dd).
0 = 0u
1 = 0d
2 = uu
3 = ud
4 = du
5 = dd
6 = Du
7 = Dd



