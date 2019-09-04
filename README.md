# Nieprzecałkowane rozkłady partonów dla procesów w LHC

Oprogramowanie powstało na potrzeby pracy magisterskiej z zakresu fizyki czastek elementarnych.

## Na start

Sklonuj repozytorium i ściągnij na swój komputer. Jeśli chcesz, aby projekt od razu pokazał swoje możliwości, dobrze posiadać konto w infrastrukturze PL-GRID, co oszczędzi czas na instalację LHA PDF oraz GSL. Plik makefile do kompilacji na klastrze PROMETHEUS w repozytorium zapusany jest jako makefile_na_klastrze.make

### Potrzebne do uruchomienia

Biblioteka LHA PDF w wersji conjamniej 6.2.1: https://lhapdf.hepforge.org/

siatki obsługiwane przez LHA: https://lhapdf.hepforge.org/pdfsets.html

Biblioteka GNU GSL w wersji conajmniej 2.4: https://www.gnu.org/software/gsl/ 

### Uruchomienie

Jeśli masz już działający makefile (detale zależą od konkretnego systemu, lokalizacji plików itp).  

Wpisz w konsolę

$ ./main.out CT10nlo

a utworzy się plik CT10nlo.dat zawierający 

x, Q2, F2=F2_q + FT + FL, F2_q, FL, FT 

w pliku main zawarte są również inne przykłady. Zamiast siatki CT10nlo można użyć dowolnie innej, dostępnej lokalnie.

Skrypt srun.sh z repozytorium uruchamia jedno zadanie na klastrze. Aby uruchomić powyższe, należy wpisać:
$ srun.sh CT14nlo 0 0 0

pamiętając o załadowaniu wcześniej wszystkich potrzebnych modułów.

## Autor 

Dominik Kasperski

Projekt dostępny na zasadzie dostępu do prac magisterskich AGH.

## Podziękowania

Dziękuję promotorowi prof. Krzystofowi Kutakowi, oraz jego współpracownikom, prof. Maciejowi Skrzypkowi, dr. Marcinowi Buremu oraz prof. Sebastianowi Sapecie, reprezentujących IFJ PAN za pomoc i dyskusję. Praca zrealizowana dzięki grantowi PL-GRID

