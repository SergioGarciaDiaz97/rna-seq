Pipeline de Análisis de ARN-Seq Automatizado 🧬

Este repositorio contiene un pipeline completo y automatizado para realizar análisis de expresión diferencial de datos de ARN-Seq, desde los archivos FASTQ crudos hasta los resultados de enriquecimiento funcional y las visualizaciones finales.

El proyecto ha sido desarrollado por Sergio García Díaz con fines didácticos y para demostrar un flujo de trabajo bioinformático reproducible, escalable y agnóstico al organismo. La base de este proyecto se sustenta en los conocimientos y el conjunto de datos del curso "Case Studies in Functional Genomics" (PH525.6x) de HarvardX, utilizando varias muestras del experimento público GSE52778 (usar human.json para dicho experimento). 

Al ser un código que puede emplearse en otros organismos, se han utilizado también muestras del experimento CARA (Characterizing Arabidopsis Root Attractions); Organismo: Arabidopsis thaliana, ecotipo estándar Columbia-0 (Col-0); Tejido: Puntas de raíz (root tips); Condición: Crecimiento bajo luz ambiental; Variable Principal: Entorno de Vuelo (ISS) vs. Control en Tierra (Ground Control). Para correr este experimento usar arabidopsis.json y configurarlo en el sh para que sea utilizado.

IMPORTANTE: Lo primero, crear una carpeta R_CODES y dentro copiar los códigos en R DESeq2_analysis.R y PLOTING_ONTOGENESIS.R, como se verá siguiendo el flujo de trabajo, al final del todo, se genera un txt de los análisis de enriquecimiento, generándose una copia en R_CODES, allí, simplemente como paso final de todo, cargar PLOTING_ONTOGENESIS.R, que lee los txt generados y genera un pdf con del enrriquecimiento (este paso no está aún automatizado, pero funciona usando el código manualmente).

--------------------------------------------------------------------
✨ Características Principales
-Automatización de Extremo a Extremo: Ejecuta el flujo de trabajo completo con un solo comando, desde la descarga de datos hasta el reporte final.

-Reproducibilidad Garantizada: Utiliza Singularity/Apptainer para encapsular todo el software en contenedores, asegurando que el análisis sea replicable en cualquier sistema.

-Configuración Flexible: Todos los parámetros, rutas, URLs de datos y diseños experimentales se gestionan a través de un único archivo config.json.

-Agnóstico al Organismo: El pipeline es totalmente versátil y ha sido probado en distintos organismos (ej. Homo sapiens, Arabidopsis thaliana). Se puede adaptar a cualquier especie simplemente modificando el archivo de configuración.

-Análisis Funcional Integrado: Realiza automáticamente un análisis de enriquecimiento de rutas Gene Ontology (GO), KEGG y Reactome para los genes significativos usando gprofiler2.

-Reporte Agregado de Calidad: Genera un informe interactivo final con MultiQC que resume los resultados de calidad de todos los pasos y muestras.

-Soporte para Diseños Complejos: Permite especificar diseños experimentales complejos en DESeq2 (ej. ~ donante + condicion) para controlar variables de confusión.

-Preparado para HPC: Incluye un script de envío run_pipeline.sh listo para ser usado en clústeres con el gestor de trabajos Slurm.

--------------------------------------------------------------------
workflow Flujo de Trabajo
El pipeline ejecuta los siguientes pasos de manera secuencial:

1) Descarga de Datos: Obtiene los genomas de referencia (FASTA y GTF) y los archivos FASTQ crudos desde las URLs especificadas.

2) Control de Calidad Inicial: Analiza la calidad de las lecturas crudas con FastQC.

3) Preprocesamiento: Limpia las lecturas de adaptadores y baja calidad utilizando Trimmomatic.

4) Alineamiento: Construye el índice del genoma y alinea las lecturas procesadas contra la referencia usando STAR.

5) Cuantificación: Genera una matriz de conteo gen-muestra con featureCounts.

6) Análisis de Expresión Diferencial:

7) Utiliza DESeq2 en R para identificar genes diferencialmente expresados.

8) Genera automáticamente los contrastes a comparar basándose en el archivo de metadatos.

9) Crea visualizaciones clave: gráficos PCA, Volcano Plots y Heatmaps.

10) Análisis de Enriquecimiento Funcional: Para cada contraste con suficientes genes significativos, identifica rutas biológicas (GO, KEGG, Reactome) sobrerrepresentadas para facilitar la interpretación biológica.

11) Generación de Reporte Final: Ejecuta MultiQC para agregar los logs y reportes de FastQC, STAR y featureCounts en un único informe HTML interactivo.
--------------------------------------------------------------------

USO
1. Prerrequisitos
Entorno Linux (preferiblemente un clúster HPC).

Singularity o Apptainer instalado.

Python 3.

Una carpeta R_CODES con el script de análisis DESeq2_analysis.R y PLOTING_ONTOGENESIS.R (ojo, este código cargarlo aparte, sobre los txt que se almacenan en sus carpetas).

2. Configuración
Modifica el archivo config.json para ajustar las rutas, muestras y parámetros de tu análisis. Presta especial atención a:

base_dir: El directorio principal de tu proyecto.

fastq_urls y genome_urls: Las fuentes de tus datos.

metadata_path: La ruta a tu archivo de metadatos.

design_formula y control_group: Para definir el análisis estadístico.

container_images: Las rutas a las imágenes de Singularity en tu sistema.

La base de datos de anotación: qué organismos son los de la base de datos que se van a emplear. Así como en los análisis de enriquecimiento.

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
--------------------------------------------------------------------

Para dudas o sugerencias de mejora, no dudes en escribirme a: sergio120897@gmail.com. ¡Estaré encantado de ayudar!
