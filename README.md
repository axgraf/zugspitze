# Zugspitze - Metagenomik-Projekt 🧬❄️  


## 🧪 Projektüberblick

Dieses Projekt stellt eine vollständige Analysepipeline für metagenomische ONT-Daten bereit, die im Rahmen von Umwelt- und Wildtiermonitoring rund um die Zugspitze erhoben wurden.

Ziel ist es, mit Hilfe des `epi2me-labs/wf-metagenomics`-Workflows gezielt das Vorkommen einzelner Spezies wie dem **Schneehuhn (*Lagopus muta*)**, dem **Birkhuhn (*Lyrurus tetrix*)** oder dem **Schneehasen (*Lepus timidus*)** nachzuweisen und zu quantifizieren.

Je nach Fragestellung stehen zwei Analysepfade zur Verfügung:

- **Minimap2-Workflow**:  
  Speziesspezifische Analyse mit eigenen Referenzgenomen. Ermöglicht gezielte Identifikation und Coverage-Analyse einzelner Zielorganismen.

- **Kraken2-Workflow**:  
  Breite taxonomische Klassifikation gegen eine umfassende Referenzdatenbank (z. B. PlusPFP-16), insbesondere für **Bakterien, Viren, Pilze, Protozoen, Pflanzen** und **Wirbeltiere** geeignet.

Zur vertieften Qualitätskontrolle steht ein **eigenes Python-Skript `genome_coverage.py`** zur Verfügung, das ausschließlich im Zusammenhang mit dem **Minimap2-Workflow** verwendet werden kann.  
Es erzeugt interaktive HTML-Berichte mit:

- detaillierten **Coverage-Statistiken**
- interaktiven **Chromosomenplots und Heatmaps**
- Navigation durch Regionen und visuelle Bewertung der Abdeckung

Die Analysepipeline eignet sich sowohl für gezielte Artensuche als auch für explorative metagenomische Studien in komplexen Umweltproben.

---

## 📄 Dokumentation

Die folgenden Links verweisen auf die offiziellen **EPI2ME-Workflows** von Oxford Nanopore Technologies, auf denen dieses Projekt basiert:

- Workflow-Dokumentation:  
  [https://epi2me.nanoporetech.com/epi2me-docs/workflows/wf-metagenomics/](https://epi2me.nanoporetech.com/epi2me-docs/workflows/wf-metagenomics/)

- Anleitung zur Nutzung eigener Referenzdatenbanken:  
  [https://epi2me.nanoporetech.com/how-to-meta-offline/](https://epi2me.nanoporetech.com/how-to-meta-offline/)

---

## 🧬 Referenzgenome

| Organismus | Taxonomy ID | Assembly | Level |
|-----------|-------------|----------|-------|
| *Lepus timidus* (Schneehase) | 62621 | [GCA_040893245.2](https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_040893245.2/) | Chromosome |
| *Lyrurus tetrix* (Birkhuhn) | 1233216 | [GCA_043882375.1](https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_043882375.1/) | Scaffold |
| *Lagopus muta* (Schneehuhn) | 64668 | [GCF_023343835.1](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_023343835.1/) | Chromosome |

> ⚠️ **Hinweis:** Die Referenzgenomen sollten regelmäßig auf Updates geprüft werden.  
> Besonders beim Birkhuhn (*Lyrurus tetrix*) liegt derzeit nur ein Scaffold-Level-Assembly vor – d. h. es sind noch keine vollständigen Chromosomen, sondern nur Contigs verfügbar.

---

## ⚙️ Installation

```bash
# Conda Umgebung erstellen
mamba create -n zugspitze_metagenome \
  -c conda-forge -c bioconda \
  nextflow seqkit kraken2 pysam pandas plotly mosdepth \
  --yes

mamba activate zugspitze_metagenome

# Docker installieren und aktivieren
sudo apt install docker.io
sudo systemctl enable --now docker
sudo usermod -aG docker $USER
newgrp docker  # damit Gruppenänderung sofort wirksam wird

# Testlauf zur Installation des Workflows
nextflow run epi2me-labs/wf-metagenomics --help
```
---

## 🧬 Vorbereitung der Referenzgenome

> ⚠️ Dieser Schritt ist **nur für den Minimap2-Workflow** erforderlich.  
> Beim [Kraken2-Workflow](#-kraken2-datenbank-setup) werden vorgefertigte Datenbanken verwendet (siehe entsprechender Abschnitt).


Lade die Referenzgenomen für Schneehuhn (*Lagopus muta*), Birkhuhn (*Lyrurus tetrix*) und Schneehase (*Lepus timidus*) herunter und erstelle eine kombinierte FASTA-Datei sowie einen Minimap2-Index.

```bash
mkdir -p referenceGenome/{Lepus_timidus,Lagopus_muta,Lyrurus_tetrix,combined_genomes,taxdump}

# Genomdateien herunterladen
wget -P referenceGenome/Lepus_timidus https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/040/893/245/GCA_040893245.2_mLepTim1.1_pri/GCA_040893245.2_mLepTim1.1_pri_genomic.fna.gz
wget -P referenceGenome/Lagopus_muta https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/023/343/835/GCF_023343835.1_bLagMut1_primary/GCF_023343835.1_bLagMut1_primary_genomic.fna.gz
wget -P referenceGenome/Lyrurus_tetrix https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/043/882/375/GCA_043882375.1_ASM4388237v1/GCA_043882375.1_ASM4388237v1_genomic.fna.gz

# FASTA-Dateien zusammenführen
cat referenceGenome/Lepus_timidus/*.gz \
    referenceGenome/Lagopus_muta/*.gz \
    referenceGenome/Lyrurus_tetrix/*.gz \
    > referenceGenome/combined_genomes/LeTim_LagMut_LyrTet.genome.fasta.gz

# Minimap2-Index erzeugen
minimap2 -t 8 -x map-ont -d referenceGenome/combined_genomes/LeTim_LagMut_LyrTet.genome.mmi referenceGenome/combined_genomes/LeTim_LagMut_LyrTet.genome.fasta.gz

# Taxonomie-Datenbank herunterladen
wget -P referenceGenome/taxdump https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.tar.gz
```

---

## 📄 ref2taxid-Datei erzeugen

Für die Klassifizierung mit Minimap2 benötigt der Workflow eine Datei, die jede FASTA-Sequenz einem Taxonomie-Identifier zuordnet. Diese Datei wird manuell erzeugt:

```bash
cd referenceGenome/combined_genomes

# Lepus timidus (TaxID: 62621)
zgrep "^>" ../Lepus_timidus/*.fna.gz | cut -d' ' -f1 | sed 's/^>//' | awk '{print $0 "\t62621"}' > lepus.tsv

# Lagopus muta (TaxID: 64668)
zgrep "^>" ../Lagopus_muta/*.fna.gz | cut -d' ' -f1 | sed 's/^>//' | awk '{print $0 "\t64668"}' > lagopus.tsv

# Lyrurus tetrix (TaxID: 1233216)
zgrep "^>" ../Lyrurus_tetrix/*.fna.gz | cut -d' ' -f1 | sed 's/^>//' | awk '{print $0 "\t1233216"}' > lyrurus.tsv

# Alle Einträge zusammenführen
cat lepus.tsv lagopus.tsv lyrurus.tsv > ref2taxid.targloci.tsv
```

---

## 🧭 Mapping-Datei für die Coverage-Analyse

> 🐍 Diese Datei wird ausschließlich für das Python-Skript [`genome_coverage.py`](./genome_coverage.py) benötigt,  
> das eine vertiefte Coverage-Analyse und Visualisierung einzelner Spezies ermöglicht.


Für die spätere Coverage-Visualisierung wird eine Mapping-Datei benötigt, die jede FASTA-Sequenz einem der referenzierten Organismen und einem benannten Chromosomen- oder Contignamen zuordnet.

**Beispielzeile:**  
`NC_064433.1    Lagopus muta    chromosome_1`

Diese Datei wird über ein Shell-Skript erzeugt. Lege dazu das Skript `generate_mapping.sh` im Verzeichnis `referenceGenome/` an.

### Beispielinhalt für `generate_mapping.sh`

```bash
#!/bin/bash
mkdir -p mappings

zgrep "^>" Lepus_timidus/*.fna.gz | cut -d' ' -f1 | sed 's/^>//' | awk -F'.' '{print $0 "\tLepus timidus\t"$1}' > mappings/lepus.tsv
zgrep "^>" Lagopus_muta/*.fna.gz | cut -d' ' -f1 | sed 's/^>//' | awk -F'.' '{print $0 "\tLagopus muta\t"$1}' > mappings/lagopus.tsv
zgrep "^>" Lyrurus_tetrix/*.fna.gz | cut -d' ' -f1 | sed 's/^>//' | awk -F'.' '{print $0 "\tLyrurus tetrix\t"$1}' > mappings/lyrurus.tsv

cat mappings/lepus.tsv mappings/lagopus.tsv mappings/lyrurus.tsv > mappings/combined_mapping.tsv
```

### Ausführung des Skripts

Führe das Mapping-Skript im Verzeichnis `referenceGenome` aus:

```bash
cd referenceGenome
chmod +x generate_mapping.sh
./generate_mapping.sh
```

## 🚀 Metagenome-Workflow von Epi2me mit Minimap2

Sobald alle Referenzdateien (Genom, Index, Taxonomie, `ref2taxid`, Mapping) vorbereitet sind, kann die Metagenomanalyse mit Minimap2 gestartet werden.

> ⚠️ Hinweis:  
> Mit dem Minimap2-Workflow können **nur genau die Spezies erkannt werden**, deren Genome zuvor im Schritt  
> [🧬 Vorbereitung der Referenzgenome](#-vorbereitung-der-referenzgenome) heruntergeladen und eingebunden wurden.  
> In diesem Projekt sind das:
> - *Lagopus muta* (Schneehuhn)
> - *Lepus timidus* (Schneehase)
> - *Lyrurus tetrix* (Birkhuhn)

```bash
nextflow run epi2me-labs/wf-metagenomics \
  --fastq /pfad/zur/fastq_pass/ \
  --classifier minimap2 \
  --reference referenceGenome/combined_genomes/LeTim_LagMut_LyrTet.genome.mmi \
  --ref2taxid referenceGenome/combined_genomes/ref2taxid.targloci.tsv \
  --taxonomy referenceGenome/taxdump/new_taxdump.tar.gz \
  --out_dir zugspitze_output \
  --keep_bam
```

### Parameterbeschreibung

- ``--fastq``  
  Pfad zum Verzeichnis mit den demultiplexierten ONT-Reads, z.B. `fastq_pass/`

- ``--classifier minimap2``  
  Nutzt Minimap2 als Klassifikator (statt z.B. Kraken2)

- ``--reference``  
  Minimap2-Indexdatei (`.mmi`) der kombinierten Genome

- ``--ref2taxid``  
  TSV-Datei, die jede Sequenz-ID (z.B. `NC_...`) einer NCBI Taxonomy-ID zuordnet

- ``--taxonomy``  
  Komprimiertes Taxonomie-Archiv (`new_taxdump.tar.gz`) von NCBI

- ``--out_dir``  
  Zielordner für alle Ausgabedateien (Berichte, Klassifikation, BAMs)

- ``--keep_bam``  
  Behalte die ausgerichteten BAM-Dateien (notwendig für Coverage-Analyse)

Das folgende Beispiel zeigt einen **einfachen manuellen Aufruf** mit Nextflow.  
▶️ **Vollautomatisierte Ausführung siehe hier:** [🧪 Automatisierte Analyse mit Minimap2](#-automatisierte-analyse-mit-minimap2)

---

## 🧪 Automatisierte Analyse mit Minimap2

Für die vollständige automatisierte Durchführung des Metagenomik-Workflows inklusive Coverage-Analyse kann 
das Bash-Skript [`run_metagenome_minimap2.sh`](./scripts/run_metagenome_minimap2.sh) verwendet werden.

### Aufruf:

```bash
./run_metagenome_minimap2.sh <FASTQ_INPUT_FOLDER> <OUTPUT_DIR> <SPECIES> [BIN_SIZE]
```

### Beispiel

```bash
./run_metagenome_minimap2.sh \
  20250626_1607_2F_PBE33297_9035ed7a/fastq_pass/ \
  zugspitze_analysis \
  "Lagopus muta" \
  1000
```

### Argumente

- `FASTQ_INPUT_FOLDER`  
  Pfad zum Ordner mit ONT-Fastq-Dateien, typischerweise der `fastq_pass/`-Ordner eines ONT-Runs.

- `OUTPUT_DIR`  
  Name oder Pfad des Ausgabeordners, in dem die Ergebnisse gespeichert werden.

- `SPECIES`  
  **Erforderlich.** Eine der folgenden Spezies (exakt in Anführungszeichen angeben):  
  - `"Lagopus muta"`  
  - `"Lepus timidus"`  
  - `"Lyrurus tetrix"`

- `BIN_SIZE` *(optional)*  
  Auflösung der Coverage-Binning-Fenster in Basen (z. B. `1000` für 1 kb).  
  Wird kein Wert angegeben, wird `1000` als Standard verwendet.

---

## 📦 Kraken2-Datenbank-Setup

Für die Analyse mit Kraken2 wird eine vorgefertigte Klassifikationsdatenbank benötigt. Die offiziellen Kraken2-Datenbanken sind vorkonfiguriert und können direkt verwendet werden.

> 🔗 Quelle: [https://benlangmead.github.io/aws-indexes/k2](https://benlangmead.github.io/aws-indexes/k2)

### Auswahl der passenden Datenbankgröße

Die Wahl der Datenbank hängt vom verfügbaren Arbeitsspeicher ab:

- **PlusPFP-16** (`~14.9GB`)  
  Empfohlen bei ≥16–32GB RAM  
  Enthält:
  - RefSeq Archaea, Bacteria, Viruses, Plasmids
  - Human, UniVec_Core
  - **zusätzlich**: RefSeq Protozoa, Pilze (Fungi), Pflanzen (Plants)

- **PlusPFP-8** (`~7.5GB`)  
  Empfohlen bei 8GB RAM  
  Gleicher Aufbau wie PlusPFP-16, aber reduzierter Umfang

> 💻 *Hinweis:* Auf einem Laptop mit 32GB RAM wurde im Projekt `PlusPFP-16` erfolgreich verwendet.

---

### Download von PlusPFP-16

```bash
mkdir -p referenceGenome/kraken2_db/k2_pluspf_16gb
cd referenceGenome/kraken2_db/k2_pluspf_16gb

wget https://genome-idx.s3.amazonaws.com/kraken/k2_pluspf_16gb_20250402.tar.gz
tar -xf k2_pluspf_16gb_20250402.tar.gz
```

---

## 🚀 Analyse mit Kraken2

Nach dem Setup der Kraken2-Datenbank und des NCBI-Taxonomiearchivs kann der Workflow mit Kraken2 als Klassifikator gestartet werden.

### Beispiel: Manuelle Ausführung

```bash
nextflow run epi2me-labs/wf-metagenomics \
  --fastq /data/20250626_DNA_Test_Schneehuhn/20250626_DNA_Test_Schneehuhn/20250626_1607_2F_PBE33297_9035ed7a/fastq_pass/ \
  --database referenceGenome/kraken2_db/k2_pluspf_16gb/ \
  --taxonomy referenceGenome/kraken2_db/taxonomy/ \
  --out_dir zugspitze_kraken2
```

### Parameterbeschreibung

- `--fastq`  
  Pfad zum Verzeichnis mit den demultiplexierten ONT-Fastq-Dateien (z.B. `fastq_pass/`)

- `--database`  
  Pfad zur entpackten Kraken2-Datenbank (z.B. `referenceGenome/kraken2_db/k2_pluspf_16gb/`)

- `--taxonomy`  
  Pfad zum Verzeichnis mit den entpackten NCBI-Taxonomiedaten (aus `new_taxdump.tar.gz`)

- `--out_dir`  
  Zielordner für die Analyseergebnisse (z.B. `zugspitze_kraken2`)

---

## 🧪 Automatisierter Workflow mit Kraken2

Für eine automatisierte Analyse mit Kraken2 kann das Skript [`run_metagenome_kraken2.sh`](./scripts/run_metagenome_kraken2.sh) verwendet werden.  
Es ruft den Nextflow-Workflow mit den im Skript definierten Datenbank- und Taxonomie-Pfaden auf.

### Aufruf

```bash
./run_metagenome_kraken2.sh <FASTQ_INPUT_FOLDER> <OUTPUT_DIR>
```

### Argumente

- `FASTQ_INPUT_FOLDER`  
  Pfad zum Ordner mit ONT-Fastq-Dateien, typischerweise `fastq_pass/`

- `OUTPUT_DIR`  
  Verzeichnis, in dem die Ergebnisse gespeichert werden (wird bei Bedarf automatisch angelegt)

### Hinweise

- Die im Skript verwendeten Pfade zur Kraken2-Datenbank und zur Taxonomie müssen existieren und korrekt gesetzt sein:

  ```bash
  KRAKEN2_DB="referenceGenome/kraken2_db/k2_pluspf_16gb/"
  TAXONOMY="referenceGenome/kraken2_db/taxonomy/"
  ```
  
---

## 📂 Beispielausgaben

Im Ordner [`examples/`](./examples/) befinden sich exemplarische HTML-Berichte, die im Rahmen dieses Projekts erzeugt wurden. Sie zeigen typische Resultate nach erfolgreicher Ausführung der jeweiligen Workflows und Analyse-Skripte.

### 🔬 Minimap2-Workflow

| Dateiname | Beschreibung |
|----------|--------------|
| [`wf-metagenomics-report_minimap2.html`](./examples/wf-metagenomics-report_minimap2.html) | Ergebnisübersicht des EPI2ME-Workflows mit Minimap2 (klassifiziert gegen eigene Referenzgenome) |
| [`coverage_report_Lagopus_muta_10kb.html`](./examples/coverage_report_Lagopus_muta_10kb.html) | Detaillierter Coverage-Report für *Lagopus muta* mit 10 kb Binning, generiert durch `genome_coverage.py` |

### 🌍 Kraken2-Workflow

| Dateiname | Beschreibung |
|----------|--------------|
| [`wf-metagenomics-report_kraken2.html`](./examples/wf-metagenomics-report_kraken2.html) | Ergebnisübersicht des EPI2ME-Workflows mit Kraken2-Datenbank (breite Klassifikation über viele Taxa) |