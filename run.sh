#!/bin/bash -l
## Nazwa zlecenia
##SBATCH -J testjob
## Liczba alokowanych węzłów
#SBATCH -N 1
## Liczba zadań per węzeł (domyślnie jest to liczba alokowanych rdzeni na węźle)
#SBATCH --ntasks-per-node=1
## Ilość pamięci przypadającej na jeden rdzeń obliczeniowy (domyślnie 4GB na rdzeń)
#SBATCH --mem-per-cpu=2GB
## Maksymalny czas trwania zlecenia (format HH:MM:SS)
#SBATCH --time=8:00:00 
## Nazwa grantu do rozliczenia zużycia zasobów
##SBATCH -A <grant_id>
## Specyfikacja partycji
#SBATCH -p plgrid
## Plik ze standardowym wyjściem błędów
#SBATCH --error="error.err"

module load plgrid/tools/gcc/4.9.2
module load plgrid/tools/lhapdf/6.2.1
module load plgrid/tools/gcc/4.9.2
module load plgrid/libs/gsl/2.4

srun -N 1 --ntasks-per-node=1 -p plgrid -t 22:00 ./main.out $1 $2 $3 $4

##dddd

