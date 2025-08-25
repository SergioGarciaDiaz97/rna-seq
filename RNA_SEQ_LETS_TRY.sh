#!/bin/bash
#SBATCH --job-name=RNAseq_full
#SBATCH --output=RNAseq_full_%j.out
#SBATCH --error=RNAseq_full_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=08:00:00
#SBATCH --mem=80G

# Exportar PATH para incluir donde esté singularity/apptainer
export PATH=$PATH:/usr/bin

# Carpeta temporal para contenedor
export APPTAINER_TMPDIR=/mnt/beegfs/home/segardia/apptainer_tmp/${SLURM_JOB_ID}
mkdir -p $APPTAINER_TMPDIR

# Define una carpeta personal para instalar paquetes de R
export R_LIBS_USER=/mnt/beegfs/home/segardia/R_libs_personal
mkdir -p $R_LIBS_USER

# Usar singularity como comando principal (puedes cambiarlo a apptainer si es necesario)
export APPTAINER_CMD=singularity

# Rutas a tus archivos principales
CONFIG_JSON=/mnt/beegfs/home/segardia/config.json
SCRIPT_PY=/mnt/beegfs/home/segardia/rnaseq_pipeline.py

echo "🚀 Iniciando pipeline RNA-seq..."
python "$SCRIPT_PY" -c "$CONFIG_JSON"
echo "✅ Pipeline terminado."