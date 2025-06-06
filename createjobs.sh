#!/bin/bash


filename=2dIsing_L5040
jobname=L5040
#touch $filename$1.slurm

cat <<EOF >> jobs/$filename-$(($1*$2)).slurm
#!/bin/bash
#SBATCH --job-name=logs/$jobname 			# Nombre del trabajo
#SBATCH --output=logs/$jobname.log   		# Archivo de registro de salida
#SBATCH --error=logs/$jobname.err    		# Archivo de registro de errores
#SBATCH --partition=QuantPhysMC   	# Nombre de la partición o cola de trabajos
#SBATCH --nodes=$1     				# Número de nodos a utilizar (puedes cambiarlo)
#SBATCH --ntasks-per-node=$2   		# Número de tareas por nodo (1 para ejecución serial)
#SBATCH --cpus-per-task=1 		# Número de CPUs por tarea (puedes cambiarlo)
#SBATCH --mem=4G      			# Memoria RAM necesaria (puedes cambiarlo)

module load lamod/coarrays/2.10 
cd ~/ising_parallel
# Comando para ejecutar tu programa
./run.sh $1 $2

EOF
