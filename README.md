Pipeline de Análisis de ARN-Seq Automatizado 🧬
Este repositorio contiene un pipeline completo y automatizado para realizar análisis de expresión diferencial de datos de ARN-Seq, desde los archivos FASTQ crudos hasta los resultados y visualizaciones finales.

El proyecto ha sido desarrollado por Sergio García Díaz con fines didácticos y para demostrar un flujo de trabajo bioinformático reproducible y escalable. La base de este proyecto se sustenta en los conocimientos y el conjunto de datos del curso "Case Studies in Functional Genomics" (PH525.6x) de HarvardX, utilizando el experimento público GSE52778.
--------------------------------------------------------------------
✨ Características Principales
-Automatización de Extremo a Extremo: Ejecuta el flujo de trabajo completo con un solo comando.
-Reproducibilidad Garantizada: Utiliza Singularity/Apptainer para encapsular el software en contenedores, asegurando que el análisis sea replicable.
-Configuración Flexible: Todos los parámetros se gestionan a través de un único archivo config.json.
-Diseño Modular: Orquestado con Python y análisis estadístico con R (DESeq2).
-Preparado para HPC: Incluye un script de envío para clústeres con Slurm.
-Agnóstico al Organismo: El pipeline es versátil y puede ser adaptado para analizar datos de cualquier organismo (ej. Homo sapiens, Arabidopsis thaliana), simplemente modificando el archivo de configuración.
--------------------------------------------------------------------
workflow Flujo de Trabajo
El pipeline ejecuta los siguientes pasos de manera secuencial:

1) Descarga de Datos: Obtiene los genomas de referencia y archivos FASTQ.

2) Control de Calidad: Analiza las lecturas con FastQC.

3) Preprocesamiento: Limpia y recorta adaptadores con Trimmomatic.

4) Alineamiento: Genera el índice del genoma y alinea las lecturas con STAR.

5) Cuantificación: Genera la matriz de conteos con featureCounts.

6) Análisis Diferencial: Ejecuta DESeq2 para obtener genes diferencialmente expresados y genera gráficos (PCA, Volcano, Heatmap).
--------------------------------------------------------------------

USO
1. Prerrequisitos
Entorno Linux (preferiblemente un clúster HPC).

Singularity o Apptainer instalado.

Python 3.

Una carpeta R_CODES con el script de análisis DESeq2_analysis.R.

2. Configuración
Modifica el archivo config.json para ajustar las rutas, muestras y parámetros de tu análisis.

3. Ejecución
Lanza el pipeline en un clúster con Slurm: sbatch run_pipeline.sh

--------------------------------------------------------------------
✨ Características Principales
Automatización de Extremo a Extremo: Ejecuta el flujo de trabajo completo con un solo comando.

Reproducibilidad Garantizada: Utiliza Singularity/Apptainer para encapsular el software en contenedores, asegurando que el análisis sea replicable.

Configuración Flexible: Todos los parámetros se gestionan a través de un único archivo config.json.

Diseño Modular: Orquestado con Python y análisis estadístico con R (DESeq2).

Preparado para HPC: Incluye un script de envío para clústeres con Slurm.

Agnóstico al Organismo: El pipeline es versátil y puede ser adaptado para analizar datos de cualquier organismo (ej. Homo sapiens, Arabidopsis thaliana), simplemente modificando el archivo de configuración.
--------------------------------------------------------------------
⚠️ Nota Técnica Importante: Parámetro sjdbOverhang
Para un alineamiento correcto con STAR, el valor de sjdbOverhang debe ser igual a la longitud de tus lecturas menos uno (longitud_lectura - 1).

Para verificar la longitud de tus lecturas (después del trimming), usa este comando:
zcat TU_MUESTRA_1.trimmed.fastq.gz | head -n 2 | tail -n 1 | wc -c ##Es un ejemplo
Si el resultado es 101, la longitud es 100pb, y el valor de sjdbOverhang debe ser 99.
--------------------------------------------------------------------
📈 Futuras Implementaciones
Este proyecto está en desarrollo activo. Las próximas mejoras planeadas incluyen:

Análisis de Enriquecimiento Funcional: Integración de análisis de Gene Ontology (GO) y rutas KEGG para la interpretación biológica de los resultados.
--------------------------------------------------------------------

Para dudas o sugerencias de mejora, no dudes en escribirme a: sergio120897@gmail.com. ¡Estaré encantado de ayudar!