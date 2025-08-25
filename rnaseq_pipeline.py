import os
import subprocess
import argparse
import json
import shutil
import gzip
import logging
import pandas as pd

def setup_logging(log_file):
    """Configura el sistema de logging para guardar en un archivo y mostrar en consola."""
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] - %(message)s",
        handlers=[
            logging.FileHandler(log_file, mode='w'),
            logging.StreamHandler()
        ]
    )

def create_directory(path):
    """Crea un directorio si no existe."""
    if not os.path.exists(path):
        os.makedirs(path)
        logging.info(f"📁 Carpeta creada: {path}")
    else:
        logging.info(f"📁 Carpeta ya existía: {path}")

def validate_container(image_path, binary):
    """Verifica si un binario está disponible dentro de un contenedor."""
    container_cmd = os.environ.get("APPTAINER_CMD", "singularity")
    logging.info(f"Usando contenedor con: {container_cmd}")
    try:
        subprocess.run([container_cmd, "exec", image_path, "which", binary],
                       check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        logging.info(f"✅ {binary} disponible en {os.path.basename(image_path)}")
    except (subprocess.CalledProcessError, FileNotFoundError):
        logging.error(f"❌ {binary} NO disponible o contenedor no encontrado en: {image_path}")
        exit(1)

def download_files(urls, destination_dir):
    """Descarga una lista de archivos desde URLs a un directorio de forma inteligente."""
    logging.info("🚀 Descargando archivos...")
    for url in urls:
        gz_filename = os.path.basename(url)
        fastq_filename = gz_filename.replace(".gz", "")
        gz_filepath = os.path.join(destination_dir, gz_filename)
        fastq_filepath = os.path.join(destination_dir, fastq_filename)
        
        if os.path.exists(fastq_filepath):
            logging.info(f"✅ {fastq_filename} ya existe. Saltando descarga.")
            if os.path.exists(gz_filepath):
                os.remove(gz_filepath)
            continue
        
        if os.path.exists(gz_filepath):
            logging.info(f"📦 {gz_filename} ya existe. Saltando descarga.")
            continue
            
        logging.info(f"⬇️ Descargando {gz_filename}...")
        try:
            subprocess.run(["wget", "-P", destination_dir, url], check=True)
        except subprocess.CalledProcessError:
            logging.warning(f"⚠️ Error descargando {gz_filename}. Ignorando...")
            continue

def unzip_files(directory):
    """Descomprime todos los archivos .gz en un directorio."""
    logging.info(f"📦 Descomprimiendo archivos en {directory}...")
    for filename in os.listdir(directory):
        if filename.endswith(".gz"):
            gz_path = os.path.join(directory, filename)
            out_path = gz_path.replace(".gz", "")
            if os.path.exists(out_path):
                logging.info(f"✅ {os.path.basename(out_path)} ya existe descomprimido. Eliminando .gz...")
                os.remove(gz_path)
                continue
            logging.info(f"🧨 Descomprimiendo {filename}...")
            try:
                with gzip.open(gz_path, 'rb') as f_in:
                    with open(out_path, 'wb') as f_out: shutil.copyfileobj(f_in, f_out)
                os.remove(gz_path)
            except Exception as e:
                logging.error(f"❌ Error descomprimiendo {filename}: {e}")
    logging.info("✅ Descompresión terminada.")

def fastqc(directory_fastq, output_dir, container_img):
    """Ejecuta FastQC en archivos FASTQ."""
    create_directory(output_dir)
    fastq_files = [os.path.join(directory_fastq, f) for f in os.listdir(directory_fastq) if f.endswith((".fastq", ".fq"))]
    if not fastq_files: return

    report_files = [f for f in os.listdir(output_dir) if f.endswith("_fastqc.html")]
    if len(report_files) >= len(fastq_files):
        logging.info("⏩ FastQC ya ejecutado para todos los archivos. Se omite.")
        return
        
    container_cmd = os.environ.get("APPTAINER_CMD", "singularity")
    command = [container_cmd, "exec", container_img, "fastqc", "-o", output_dir] + fastq_files
    logging.info("🚀 Ejecutando FastQC...")
    try:
        subprocess.run(command, check=True, capture_output=True, text=True)
        logging.info("✅ FastQC completado.")
    except subprocess.CalledProcessError as e:
        logging.error(f"❌ Error durante la ejecución de FastQC: {e.stderr}")

def download_genome_files(reference_dir, genome_urls):
    """Descarga y descomprime los archivos del genoma de referencia."""
    logging.info("🚀 Descargando genoma de referencia...")
    download_files(genome_urls, reference_dir)
    unzip_files(reference_dir)

def prepare_adapters(adapters_file, adapter_url):
    """Descarga el archivo de adaptadores si no existe."""
    if not os.path.exists(adapters_file):
        create_directory(os.path.dirname(adapters_file))
        logging.info(f"⬇️ Descargando adaptadores desde {adapter_url}...")
        try:
            subprocess.run(["wget", "-O", adapters_file, adapter_url], check=True)
            logging.info(f"✅ Adaptadores descargados en {adapters_file}")
        except Exception as e:
            logging.error(f"❌ Error al descargar adaptadores: {e}")
            return None
    else:
        logging.info(f"✅ Adaptadores ya existen: {adapters_file}")
    return adapters_file

def run_trimmomatic(input_dir, output_dir, container_img, adapters_file, seq_type, threads):
    """Ejecuta Trimmomatic para limpiar las lecturas."""
    create_directory(output_dir)
    container_cmd = os.environ.get("APPTAINER_CMD", "singularity")
    if seq_type == "paired-end":
        logging.info("Trimmomatic se ejecutará en modo Paired-End (PE).")
        samples = {}
        for f in sorted(os.listdir(input_dir)):
            if f.endswith(("_1.fastq", "_R1.fastq")):
                sample = f.replace("_1.fastq", "").replace("_R1.fastq", "")
                r1 = os.path.join(input_dir, f)
                r2 = os.path.join(input_dir, f.replace("_1.fastq", "_2.fastq").replace("_R1.fastq", "_R2.fastq"))
                if os.path.exists(r2): samples[sample] = (r1, r2)
        
        for sample, (r1, r2) in samples.items():
            out1 = os.path.join(output_dir, f"{sample}_1.trimmed.fastq.gz")
            out2 = os.path.join(output_dir, f"{sample}_2.trimmed.fastq.gz")
            if os.path.exists(out1) and os.path.exists(out2):
                logging.info(f"⏩ Trimmomatic ya hecho para {sample}, saltando.")
                continue
            
            out1_unp = os.path.join(output_dir, f"{sample}_1.unpaired.fastq.gz")
            out2_unp = os.path.join(output_dir, f"{sample}_2.unpaired.fastq.gz")
            cmd = [container_cmd, "exec", container_img, "trimmomatic", "PE", "-threads", str(threads), r1, r2, out1, out1_unp, out2, out2_unp, f"ILLUMINACLIP:{adapters_file}:2:30:10", "LEADING:3", "TRAILING:3", "SLIDINGWINDOW:4:15", "MINLEN:36"]
            
            logging.info(f"🚀 Ejecutando Trimmomatic para {sample}...")
            try:
                subprocess.run(cmd, check=True, capture_output=True, text=True)
                logging.info(f"✅ Trimmomatic completado: {sample}")
                if os.path.exists(out1_unp): os.remove(out1_unp)
                if os.path.exists(out2_unp): os.remove(out2_unp)
            except subprocess.CalledProcessError as e:
                logging.error(f"❌ Error en Trimmomatic para {sample}: {e.stderr}")
                exit(1)
    else: # single-end
        pass

def get_reference_files(reference_dir):
    """Obtiene las rutas de los archivos FASTA y GTF de un directorio."""
    fasta_file, gtf_file = None, None
    for fname in os.listdir(reference_dir):
        if fname.endswith((".fa", ".fna", ".fasta")):
            fasta_file = os.path.join(reference_dir, fname)
        elif fname.endswith(".gtf"):
            gtf_file = os.path.join(reference_dir, fname)
    return fasta_file, gtf_file

def build_genome_index(fasta_file, gtf_file, index_dir, star_container, threads, sjdb_overhang):
    """Construye el índice del genoma para STAR."""
    container_cmd = os.environ.get("APPTAINER_CMD", "singularity")
    if os.path.exists(os.path.join(index_dir, "SA")):
        logging.info(f"⏩ Índice STAR ya existe en {index_dir}. Saltando.")
        return
    cmd = [container_cmd, "exec", star_container, "STAR", "--runThreadN", str(threads), "--runMode", "genomeGenerate", "--genomeDir", index_dir, "--genomeFastaFiles", fasta_file, "--sjdbGTFfile", gtf_file, "--sjdbOverhang", str(sjdb_overhang)]
    logging.info(f"Construyendo índice STAR en {index_dir} ...")
    try:
        subprocess.run(cmd, check=True, capture_output=True, text=True)
        logging.info(f"✅ Índice STAR generado en {index_dir}")
    except subprocess.CalledProcessError as e:
        logging.error(f"❌ Error al construir el índice STAR: {e.stderr}")
        exit(1)

def alignment(trimmed_dir, star_container, reference_dir, output_dir, seq_type, threads):
    """Realiza el alineamiento de las lecturas con STAR."""
    create_directory(output_dir)
    container_cmd = os.environ.get("APPTAINER_CMD", "singularity")
    if seq_type == "paired-end":
        samples = {}
        for f in sorted(os.listdir(trimmed_dir)):
            if f.endswith(("_1.trimmed.fastq.gz", "_R1.trimmed.fastq.gz")):
                sample = f.replace("_1.trimmed.fastq.gz", "").replace("_R1.trimmed.fastq.gz", "")
                r1 = os.path.join(trimmed_dir, f)
                r2 = os.path.join(trimmed_dir, f.replace("_1.trimmed.fastq.gz", "_2.trimmed.fastq.gz").replace("_R1.trimmed.fastq.gz", "_R2.trimmed.fastq.gz"))
                if os.path.exists(r2): samples[sample] = (r1, r2)
        
        for sample, (r1_path, r2_path) in samples.items():
            final_bam_file = os.path.join(output_dir, f"{sample}_Aligned.sortedByCoord.out.bam")
            if os.path.exists(final_bam_file):
                logging.info(f"⏩ El archivo BAM final para {sample} ya existe. Saltando alineamiento.")
                continue
            
            logging.info(f"Procesando muestra pareada: {sample}")
            out_prefix = os.path.join(output_dir, sample + "_")
            cmd = [container_cmd, "exec", star_container, "STAR", "--runThreadN", str(threads), "--genomeDir", reference_dir, "--outFileNamePrefix", out_prefix, "--outSAMtype", "BAM", "SortedByCoordinate", "--readFilesCommand", "zcat", "--readFilesIn", r1_path, r2_path]
            logging.info(f"Ejecutando STAR para {sample}...")
            try:
                subprocess.run(cmd, check=True, text=True, capture_output=True)
                logging.info(f"✅ STAR finalizó para {sample}")
            except subprocess.CalledProcessError as e:
                logging.error(f"❌ ERROR: STAR falló para {sample}\n{e.stderr}")
                exit(1)
    else: # single-end
        pass

def generate_count_matrix(bam_dir, gtf_file, output_file, container_img, seq_type, threads):
    """Genera la matriz de conteo con featureCounts."""
    if os.path.exists(output_file):
        logging.info(f"⏩ La matriz de conteo '{os.path.basename(output_file)}' ya existe. Saltando este paso.")
        return
        
    container_cmd = os.environ.get("APPTAINER_CMD", "singularity")
    bam_files = sorted([os.path.join(bam_dir, f) for f in os.listdir(bam_dir) if f.endswith(".bam")])
    if not bam_files:
        logging.warning("⚠️ No se encontraron archivos BAM.")
        return
        
    cmd_list = [container_cmd, "exec", container_img, "featureCounts", "-T", str(threads), "-a", gtf_file, "-o", output_file]
    if seq_type == "paired-end":
        logging.info("featureCounts se ejecutará en modo Paired-End.")
        cmd_list.extend(["-p"])
    else:
        logging.info("featureCounts se ejecutará en modo Single-End.")
    
    cmd_list.extend(bam_files)
    logging.info("Ejecutando featureCounts...")
    try:
        subprocess.run(cmd_list, check=True, text=True, capture_output=True)
        logging.info(f"✅ Matriz de conteo generada: {output_file}")
    except subprocess.CalledProcessError as e:
        logging.error(f"❌ Error en featureCounts:\n{e.stderr}")
        exit(1)

def generate_contrasts_from_metadata(metadata_path, control_group):
    """Genera automáticamente contrastes a partir del archivo de metadatos."""
    try:
        df = pd.read_csv(metadata_path)
        all_conditions = df['condition'].unique()
        contrasts = []
        for condition in all_conditions:
            if condition != control_group:
                contrasts.append([condition, control_group])
        logging.info(f"Contraste(s) generado(s) automáticamente: {contrasts}")
        return contrasts
    except Exception as e:
        logging.error(f"❌ No se pudo generar los contrastes desde {metadata_path}: {e}")
        return None

def run_deseq2(config, base_dir, counts_file_path, alignments_dir, gtf_file):
    """Ejecuta el script de R para el análisis con DESeq2."""
    logging.info("--- PASO 8: Análisis de Expresión Diferencial con R/DESeq2 ---")
    setup_config = config.get("project_setup", {})
    deseq2_config = config.get("deseq2_experiment", {})
    scripts_config = config.get("scripts", {})
    images_config = config.get("container_images", {})
    annotation_config = config.get("annotation", {})
    tool_params = config.get("tool_parameters", {})
    analysis_thresholds = tool_params.get("analysis_thresholds", {})
    
    metadata_file = deseq2_config.get("metadata_path")
    control_group = deseq2_config.get("control_group")
    
    contrasts = generate_contrasts_from_metadata(metadata_file, control_group)
    
    if not contrasts:
        logging.warning("⚠️ No se pudieron generar los contrastes. Saltando análisis DESeq2.")
        return

    counting_method = setup_config.get("counting_method", "featureCounts").lower()
    r_script = scripts_config.get("r_script_path")
    r_container = images_config.get("r_deseq2")
    output_dir = os.path.join(base_dir, setup_config.get("deseq2_results_dir", "DESEQ2_RESULTS"))
    
    for i, contrast in enumerate(contrasts):
        treatment_group = contrast[0]
        control_group_contrast = contrast[1]
        
        logging.info(f"🚀 Ejecutando DESeq2 (Comparación {i+1}/{len(contrasts)}): {treatment_group} vs {control_group_contrast}")

        cmd = [
            os.environ.get("APPTAINER_CMD", "singularity"), "exec", r_container, "Rscript", r_script,
            "--counting_method", counting_method,
            "--counts_file", counts_file_path if counts_file_path else "NULL",
            "--bam_dir", alignments_dir,
            "--gtf_file", gtf_file,
            "--metadata_file", metadata_file,
            "--design_formula", deseq2_config.get("design_formula"),
            "--control_group", control_group,
            "--organism_db", annotation_config.get("organism_db"),
            "--key_type", annotation_config.get("key_type"),
            "--output_dir", output_dir,
            "--padj_threshold", str(analysis_thresholds.get("padj", 0.05)),
            "--log2fc_threshold", str(analysis_thresholds.get("log2fc", 1.0)),
            "--treatment_group", treatment_group,
            "--control_group_contrast", control_group_contrast
        ]
        
        try:
            result = subprocess.run(cmd, check=True, capture_output=True, text=True)
            logging.info(f"✅ Análisis para {treatment_group} vs {control_group_contrast} completado.")
            logging.debug(f"--- Salida de R ---\n{result.stdout}")
            if result.stderr:
                logging.warning(f"--- Mensajes de R (stderr) ---\n{result.stderr}")
        except subprocess.CalledProcessError as e:
            logging.error(f"❌ Error en DESeq2 para el contraste {treatment_group} vs {control_group_contrast}:")
            logging.error(f"--- Salida de R (stdout) ---\n{e.stdout}")
            logging.error(f"--- Error de R (stderr) ---\n{e.stderr}")
            exit(1)

def main():
    """Función principal que orquesta todo el pipeline."""
    parser = argparse.ArgumentParser(description="Pipeline de RNA-seq de extremo a extremo")
    parser.add_argument("-c", "--config", required=True, help="Archivo de configuración JSON.")
    args = parser.parse_args()
    
    try:
        with open(args.config, 'r') as f: config = json.load(f)
    except FileNotFoundError:
        print(f"❌ No se encontró el archivo de configuración: {args.config}")
        return
    except json.JSONDecodeError as e:
        print(f"❌ Error leyendo JSON: {e}")
        return
    
    base_dir = config.get("project_setup", {}).get("base_dir", ".")
    log_file = os.path.join(base_dir, "pipeline.log")
    if not os.path.exists(base_dir):
        os.makedirs(base_dir)
    setup_logging(log_file)
    
    logging.info("======================================================")
    logging.info("🚀 INICIANDO NUEVA EJECUCIÓN DEL PIPELINE DE RNA-SEQ 🚀")
    logging.info("======================================================")
    
    setup_params = config.get("project_setup", {})
    source_params = config.get("source_data", {})
    images = config.get("container_images", {})
    tool_params = config.get("tool_parameters", {})
    
    seq_type = setup_params.get("sequencing_type", "paired-end").lower()
    threads = 8 
    counting_method = setup_params.get("counting_method", "featureCounts").lower()

    fastq_dir = os.path.join(base_dir, "FASTQ_FILES")
    fastqc_dir = os.path.join(base_dir, "FASTQC")
    reference_dir = os.path.join(base_dir, "REFERENCE_GENOMES_FILES")
    trimmed_dir = os.path.join(base_dir, "TRIMMED_READS")
    alignments_dir = os.path.join(base_dir, "ALIGMENTS")
    adapters_file = os.path.join(base_dir, "adapters/TruSeq3-PE.fa")
    
    for d in [fastq_dir, fastqc_dir, reference_dir, trimmed_dir, alignments_dir]:
        create_directory(d)
    
    logging.info(f"Iniciando pipeline en modo '{counting_method}'...")
    logging.info("--- PASO 1-6: Procesamiento Bioinformático ---")
    download_files(source_params.get("fastq_urls", []), fastq_dir)
    unzip_files(fastq_dir)
    if images.get("fastqc"):
        validate_container(images.get("fastqc"), "fastqc")
        fastqc(fastq_dir, fastqc_dir, images.get("fastqc"))
    download_genome_files(reference_dir, source_params.get("genome_urls", []))
    fasta_file, gtf_file = get_reference_files(reference_dir)
    if not fasta_file or not gtf_file:
        logging.error("❌ No se encontraron archivos de referencia FASTA o GTF. Saliendo.")
        return
        
    trimmomatic_img = images.get("trimmomatic")
    adapter_url = tool_params.get("trimmomatic", {}).get("adapter_fasta_url")
    if trimmomatic_img and adapter_url:
        adapters_path = prepare_adapters(adapters_file, adapter_url)
        if adapters_path:
            validate_container(trimmomatic_img, "trimmomatic")
            run_trimmomatic(fastq_dir, trimmed_dir, trimmomatic_img, adapters_path, seq_type, threads)
            
    star_img = images.get("star")
    sjdb_overhang = tool_params.get("star", {}).get("sjdbOverhang", 100)
    if star_img:
        validate_container(star_img, "STAR")
        build_genome_index(fasta_file, gtf_file, reference_dir, star_img, threads, sjdb_overhang)
        alignment(trimmed_dir, star_img, reference_dir, alignments_dir, seq_type, threads)
    
    counts_file_path = os.path.join(base_dir, "counts.txt")
    logging.info("--- PASO 7: Generación de Matriz de Conteo ---")
    if counting_method == "featurecounts":
        featurecounts_img = images.get("featurecounts")
        if featurecounts_img:
            validate_container(featurecounts_img, "featureCounts")
            generate_count_matrix(alignments_dir, gtf_file, counts_file_path, featurecounts_img, seq_type, threads)
    elif counting_method == "summarizeoverlaps":
        logging.info("⏩ Saltando conteo con featureCounts. R se encargará del conteo.")
        counts_file_path = None
    else:
        logging.error(f"❌ Método de conteo '{counting_method}' no reconocido.")
        return

    run_deseq2(config, base_dir, counts_file_path, alignments_dir, gtf_file)

    logging.info("✅ ENHORABUENA. Pipeline de RNA-seq completado exitosamente.")

if __name__ == "__main__":
    main()