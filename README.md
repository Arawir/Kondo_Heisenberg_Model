# Kondo_Heisenberg_Model

## Kompilacja

TODO 

## Uruchomienie

Aby wybrać współczynniki opisujące eksperyment należy uruchomić program z odpowiednimi flagami, na przykład:

$ ./Kondo_Heisenberg U=1E-9 K=1.4 Jh=0.002 t00=1 L=4 state=0-7-0-3-6-6 Mu=2

##Dostępne flagi:
1) L=x - ustawia ilość węzłów w modelu
2) PBC=x - można włączyć lub wyłączyć PBC
3) t00=x , U=x , K=x , Jh=x , Mu=x - ustawia parametry modelu (Mu to potencjał chemiczny)
4) t=x - czas urojonej ewolucji czasowej po zakonczeniu DMRG (na razie zostawiam)
5) conSz=x , conN=x - wymuszamy zachowanie spinu lub liczby cząstek (wartości takie same jak dla stanu początkowego)
6) state=S-S-S-S - ustawiamy stan początkowy; jako S możemy podać numer stanu (0-7) lub jego opis (0_U, 0_D, ... , UD_D). Kolejne węzły oddzielone są myślnikiem bez spacji

##Lista potencjalnie dostępnych stanów węzła. Można wstawić je do argumentu state=S-S-S-S-S...-S w miejsce "S". Można wstawiać zarówno numery (0-7) oraz postać opisową (0_U, ..., UD_D).
0 = 0_U
1 = 0_D
2 = U_U
3 = U_D
4 = D_U
5 = D_D
6 = UD_U
7 = UD_D



