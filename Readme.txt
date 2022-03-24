Kompilace pomoci cmake (v prikazovem radku):

1) Vytvorit adresar, kde bude probihat kompilace, napr.
   jsme v adresari FluidFlow, zadame

   mkdir -p Build/Release

  Vytvori se nam adresar Build a v nem podadresar Release.

2) Prejdeme do adresare Release

   cd Build/Release

3) Zadame prikaz

   cmake -D CMAKE_BUILD_TYPE=Release ../../src/

   Tim ziskame Makefile

4) Zadame prikaz

   make (pripadne make -j #proc, kde #proc je pocet procesoru pouzitych pro
   kompilaci)

   Tim ziskame spustitelny soubor fluidFlow
   

Poznamka: pro kompilaci v debug modu zadame Debug misto Release v prikazu cmake
