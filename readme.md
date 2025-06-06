# RNAseq con genoma de referencia

Este pipeline realiza un análisis transcriptómico completo a partir de lecturas crudas en formato FASTQ. Está diseñado para ser reproducible, modular y ejecutarse tanto en local como en entornos HPC con Slurm.

## 1. Requisitos e instalación

Este proyecto usa [Snakemake](https://snakemake.readthedocs.io) y entornos `conda` para manejar dependencias. Si no tienes Snakemake instalado:

### Instalar Miniconda (si no lo tienes)
```bash
# En Linux
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
```bash
### Crear entorno con snakemake

```bash
conda create -n snakemake_env -c conda-forge -c bioconda snakemake
conda activate snakemake_env
```

### Descarga de proyecto

Primero, crea un directorio para tu proyecto (si aún no lo has hecho) y clona el repositorio desde GitHub dentro de él:

```bash
mkdir -p mi_proyecto
cd mi_proyecto
git clone https://github.com/trippv/2025_snakemake_tests.git
cd 2025_snakemake_tests
```


## 2. Proceso
El flujo de trabajo automatiza las siguientes etapas:

1. Control de calidad y recorte: fastp

2. Extracción de sitios de empalme: hisat2_extract_splice_sites.py

3. Indexado del genoma: hisat2-build

4. Alineamiento de lecturas: hisat2

5. Estadísticas de alineamiento: samtools

6. Cuantificación por ensamblaje: stringtie

7. Unión y comparación con anotación: stringtie --merge, gffcompare

8. Cuantificación final y generación de matrices: prepDE.py

9. Generación de metadata y análisis PCA: script en Rmarkdown

10. Todos los pasos se ejecutan dentro de entornos conda definidos en la carpeta envs/.




## Estructura de archivos
```bash
project/
├── config/
│   ├── config.yaml        # Rutas de genoma y archivo de muestras
│   └── samples.tsv        # Información de muestras
├── data/                  # Carpeta de trabajo
├── testdata/              # (Ejemplo) Datos de prueba
├── results/               # Salida del análisis
├── rules/                 # Reglas individuales de Snakemake
├── envs/                  # Entornos conda por herramienta
├── scripts/               # Scripts auxiliares (e.g., PCA)
├── profiles/slurm/        # Configuración de ejecución en SLURM
└── Snakefile              # Punto de entrada del workflow
```

# Definir objetos

## lecturas

Crea el archivo samples.tsv dentro de `config/`. Ejemplo de contenido 


### Ejemplo de `samples.tsv`

| sample_id | group | fastq1                               | fastq2                               | include |
|-----------|--------|--------------------------------------|--------------------------------------|---------|
| batch1    | batch | testdata/raw/batch1_chrI_1.fastq     | testdata/raw/batch1_chrI_2.fastq     | 1       |
| chem1     | chem  | testdata/raw/chem1_chrI_1.fastq      | testdata/raw/chem1_chrI_2.fastq      | 1       |


>include: define si la muestra se procesa (1) o se omite (0)

> Cada columna tiene que estar separada por tabulador (\t)

## Archivo de configuración

Archivo de configuración ubicado en `confi/congi.yaml`. Ejemplo

```bash
samples_file: "config/samples.tsv"
genome_fasta: "testdata/genome/genome.fna"
gtf: "testdata/genome/genes.gtf"
```

## Ejecutar el pipeline localmente

```bash
snakemake --use-conda --cores 8

```