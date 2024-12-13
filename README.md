```markdown
# Microarray Data Analysis Pipeline

This document outlines the steps for processing and analyzing microarray data using R. The pipeline includes data preprocessing, quality control, and differential expression analysis.

## GEO Dataset

This analysis uses the GEO dataset GSE62817.

## Important:
You will need to create a separate metadata file (`metadata.csv`) for this dataset, which should contain the sample information such as sample IDs, condition (Metastatic, Non-metastatic, Normal), and replicate number.

## Requirements

Ensure the following R libraries are installed:

- `affy`
- `dplyr`
- `tidyr`
- `ggplot2`
- `limma`
- `affyPLM`
- `mouse4302cdf`
- `mouse4302.db`

To install any missing packages, use the command:
```R
install.packages("package_name")
# For Bioconductor packages:
BiocManager::install("package_name")
```

## Pipeline Steps

### 1. Set Up
- Define the directory containing CEL files and dataset ID.
- List and read CEL files using the appropriate CDF file for the chip.

### 2. Load Libraries
```R
# Load required libraries for different steps of the analysis
library(affy)            # For processing Affymetrix data
library(dplyr)           # For data manipulation
library("tidyr")         # Data tidying
library(ggplot2)         # For data visualization
library(limma)           # For linear modeling and differential expression
library(affyPLM)         # For probe-level modeling
library(mouse4302cdf)    # CDF file for the Mouse4302 chip
library(mouse4302.db)    # Annotation package for the Mouse4302 chip
```

### 3. Read and Explore Data
- Read CEL files and assign CDF file for the chip.
- Explore the Affy object to retrieve sample names, probe names, and annotations.

```R
# Read the CEL files and assign the correct CDF file for the mouse4302 chip
affydata = ReadAffy(celfile.path = celfiledirectory, cdfname = "mouse4302cdf")

# Explore the structure of the affydata object
class(affydata)  # Check the class of the object
sampleNames(affydata)  # Retrieve sample names
featureNames(affydata)  # Retrieve probe names
```

### 4. Preprocess Data
- Preprocess microarray data using the `threestep` method with background correction, normalization, and summarization.
- Extract expression values.

```R
# Preprocess the microarray data with background correction, normalization, and summarization
eset = threestep(
  affydata,
  background.method = "IdealMM",  # Background correction method
  normalize.method = "quantile",  # Normalization method
  summary.method = "average.log" # Summarization method
)

# Extract expression values from the preprocessed data
data = exprs(eset)
```

### 5. Quality Control
- Generate histograms and boxplots for raw and preprocessed data.

```R
# Generate histograms for raw and preprocessed data
hist(affydata, main="Histogram")  # Raw data histogram
hist(eset, main="Histogram")      # Preprocessed data histogram

# Generate boxplot for preprocessed data to visualize distribution
boxplot(eset, main="Boxplot", las=2)
```

### 6. Probe-to-Gene Mapping
- Map probe IDs to gene symbols using the annotation package `mouse4302.db`.
- Merge the expression data with gene symbols.

```R
# Map probe IDs to gene symbols using the mouse4302 annotation package
mapper = mouse4302SYMBOL
map.df = as.data.frame(mapper)

# Merge expression data with gene symbols
data2 = merge(data, map.df, by="probe_id", all.x = TRUE)

# Rearrange columns to place gene symbols first and remove probe_id column
data2 <- data2 %>% dplyr::select(symbol, everything(), -probe_id)
```

### 7. Aggregate Expression Data
- Aggregate expression data by gene symbols and clean column names.

```R
# Aggregate expression data by gene symbols (mean expression per gene)
exp_agg = aggregate(data2[-1], by=list(data2$symbol), FUN=mean)

# Clean up the row and column names for the aggregated data
rownames(exp_agg) = exp_agg$Group.1
names(exp_agg) = sub("_.*", "", names(exp_agg))  # Clean column names
exp = exp_agg[-1]  # Remove the extra group column
```

### 8. Load Metadata
- Read metadata and ensure sample IDs match with expression data.

```R
# Read metadata from the CSV file
metadata = read.csv("metadata.csv")

# Separate Condition and Replica columns from the Condition column
metadata <- metadata %>% separate(Condition, into=c("Condition", "Replica"), sep="_", remove=FALSE)

# Check that the sample IDs in the metadata match the expression data column names
all(colnames(exp) == metadata$Sample_id)
```

### 9. Differential Expression Analysis
- Define the experimental groups and create a design matrix.
- Fit a linear model, define contrasts, and extract differentially expressed genes (DEGs).

```R
# Define experimental groups and create a design matrix
groups = factor(metadata$Condition, levels = c("Metastatic", "Non_metastatic", "Normal"))
design = model.matrix(~0 + groups)

# Fit the linear model to the expression data
fit = lmFit(exp, design)

# Define contrasts to compare conditions
contrast.matrix = makeContrasts(
  metastaticVsCTR = Metastatic - Normal,
  NonmetastaticVsCTR = Non_metastatic - Normal,
  levels = design
)

# Apply contrasts and run empirical Bayes analysis to find DEGs
fit2 = contrasts.fit(fit, contrast.matrix)
fit2 = eBayes(fit2)

# Extract and save the DEGs for metastatic vs control condition
degs_metastaic = topTable(fit2, coef="metastaticVsCTR", number=Inf, adjust="BH")
write.csv(degs_metastaic, "DEGs_metastatic_vs_control.csv", row.names=FALSE)
```

## Outputs
- Histograms and boxplots for raw and preprocessed data.
- Gene expression matrix mapped to symbols.
- Differentially expressed genes saved as `DEGs_metastatic_vs_control.csv`.

## Notes
- Ensure the CEL files and metadata file are correctly placed in the specified directory.
```
