```R
# Microarray Data Analysis Pipeline
This document outlines the steps for processing and analyzing microarray data using R. The pipeline includes data preprocessing, quality control, and differential expression analysis.

## Metadata Information

The metadata file for this dataset (GEO: GSE62817) contains the following columns:

- `Sample_id`: Corresponds to the sample IDs (e.g., GSM1533823, GSM1533824, etc.).
- `Condition`: Indicates the condition for each sample (Metastatic, Non-metastatic, or Normal).
- `Replica`: Indicates the replicate number (rep1, rep2, etc.).

### Example of `metadata.csv`:

```csv
Sample_id,Condition,Replica
GSM1533823,Metastatic,rep1
GSM1533824,Metastatic,rep2
GSM1533825,Metastatic,rep3
...
```

Ensure the sample IDs in the metadata file match those in the expression data. The metadata file can be read and processed with:

```R
metadata = read.csv("metadata.csv")
metadata <- metadata %>% separate(Condition, into=c("Condition", "Replica"), sep="_", remove=FALSE)
```
# Microarray Data Analysis Pipeline

# This document outlines the steps for processing and analyzing microarray data using R.
# The pipeline includes data preprocessing, quality control, and differential expression analysis.

## Requirements
# Ensure the following R libraries are installed:

# - `affy`: For processing Affymetrix data
# - `dplyr`: For data manipulation
# - `tidyr`: For data tidying
# - `ggplot2`: For data visualization
# - `limma`: For linear modeling and differential expression
# - `affyPLM`: For probe-level modeling
# - `mouse4302cdf`: CDF file for the Mouse4302 chip
# - `mouse4302.db`: Annotation package for the Mouse4302 chip

# To install any missing packages, use the command:
```R
install.packages("package_name")
# For Bioconductor packages:
BiocManager::install("package_name")
```

## Pipeline Steps

### 1. Set Up
# Define the directory containing CEL files and dataset ID.
# List and read CEL files using the appropriate CDF file for the chip.

### 2. Load Libraries
```R
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
# Read CEL files and assign CDF file for the chip.
# Explore the Affy object to retrieve sample names, probe names, and annotations.

```R
affydata = ReadAffy(celfile.path = celfiledirectory, cdfname = "mouse4302cdf")  # Reading the CEL files
class(affydata)  # Check the class of the data
sampleNames(affydata)  # Retrieve sample names
featureNames(affydata)  # Retrieve feature names
```

### 4. Preprocess Data
# Preprocess microarray data using the `threestep` method with background correction, normalization, and summarization.
# Extract expression values.

```R
eset = threestep(
  affydata,
  background.method = "IdealMM",  # Background correction method
  normalize.method = "quantile",  # Normalization method
  summary.method = "average.log"  # Summarization method
)
data = exprs(eset)  # Extract expression values from the processed data
```

### 5. Quality Control
# Generate histograms and boxplots for raw and preprocessed data.

```R
hist(affydata, main="Histogram")  # Histogram of raw data
hist(eset, main="Histogram")  # Histogram of preprocessed data
boxplot(eset, main="Boxplot", las=2)  # Boxplot of preprocessed data
```

### 6. Probe-to-Gene Mapping
# Map probe IDs to gene symbols using the annotation package `mouse4302.db`.
# Merge the expression data with gene symbols.

```R
mapper = mouse4302SYMBOL  # Mapping probes to gene symbols
map.df = as.data.frame(mapper)  # Convert mapping to a data frame
data2 = merge(data, map.df, by="probe_id", all.x = TRUE)  # Merge expression data with gene symbols
data2 <- data2 %>% dplyr::select(symbol, everything(), -probe_id)  # Select relevant columns
```

### 7. Aggregate Expression Data
# Aggregate expression data by gene symbols and clean column names.

```R
exp_agg = aggregate(data2[-1], by=list(data2$symbol), FUN=mean)  # Aggregate expression by gene symbol
rownames(exp_agg) = exp_agg$Group.1  # Set row names to gene symbols
names(exp_agg) = sub("_.*", "", names(exp_agg))  # Clean up column names
exp = exp_agg[-1]  # Remove unnecessary column (Group.1)
```

### 8. Load Metadata
# Read metadata and ensure sample IDs match with expression data.

```R
metadata = read.csv("metadata.csv")  # Read metadata file
metadata <- metadata %>% separate(Condition, into=c("Condition", "Replica"), sep="_", remove=FALSE)  # Separate Condition and Replica
all(colnames(exp) == metadata$Sample_id)  # Ensure sample IDs match between metadata and expression data
```

### 9. Differential Expression Analysis
# Define the experimental groups and create a design matrix.
# Fit a linear model, define contrasts, and extract differentially expressed genes (DEGs).

```R
groups = factor(metadata$Condition, levels = c("Metastatic", "Non_metastatic", "Normal"))  # Define experimental groups
design = model.matrix(~0 + groups)  # Create design matrix
fit = lmFit(exp, design)  # Fit linear model
contrast.matrix = makeContrasts(
  metastaticVsCTR = Metastatic - Normal,  # Define contrasts for metastatic vs control
  NonmetastaticVsCTR = Non_metastatic - Normal,  # Define contrasts for non-metastatic vs control
  levels = design  # Set levels for the contrasts
)
fit2 = contrasts.fit(fit, contrast.matrix)  # Apply contrasts to the fit
fit2 = eBayes(fit2)  # Apply empirical Bayes moderation
degs_metastaic = topTable(fit2, coef="metastaticVsCTR", number=Inf, adjust="BH")  # Extract DEGs for metastatic vs control
write.csv(degs_metastaic, "DEGs_metastatic_vs_control.csv", row.names=FALSE)  # Save DEGs to a CSV file
```

## Outputs
# - Histograms and boxplots for raw and preprocessed data.
# - Gene expression matrix mapped to symbols.
# - Differentially expressed genes saved as `DEGs_metastatic_vs_control.csv`.

## Notes
# - Ensure the CEL files and metadata file are correctly placed in the specified directory.
# - Use `sessionInfo()`


The added comments provide a detailed explanation of each step in the pipeline, making the code easier to understand for someone following along.
