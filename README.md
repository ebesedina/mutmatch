# mutmatch

`mutmatch` is an R package designed to estimate somatic selection for a gene using various methods including neighbors baseline of mutation rates or low-impact mutations estimated with CADD score.

## Installation

### Setting Up Conda Environment

Before installing the package, it's recommended to set up an isolated conda environment to ensure compatibility and avoid conflicts with other packages.

1. If you don't have `conda` or `mamba` installed, download and install [Miniconda](https://docs.conda.io/en/latest/miniconda.html) or [Anaconda](https://www.anaconda.com/products/distribution).

2. Create a new conda environment using the provided `environment.yml`:

```bash
conda env create -f environment.yml
```

3. Activate the newly created conda environment:

```bash
conda activate mutmatch
```

### Installing mutmatch R package

You can install `mutmatch` from [GitHub](https://github.com/ebesedina/mutmatch) using:

```R
# If the devtools package is not installed, uncomment the following line:
# install.packages("devtools")
devtools::install_github("ebesedina/mutmatch")
```

## Features

- Estimate somatic selection forces using a gene's neighboring baseline mutation rates.
- Use low-impact mutations estimated with CADD score as an alternative baseline.
- Comprehensive filtering options for genomic regions.
- Flexible statistical modeling for selection.

## Usage

### Estimating Selection Using Neighboring Genes

The following example demonstrates how to estimate selection using neighboring genes as a neutral mutational rate baseline for the `KRAS` gene:

```R
# Load the mutmatch library
library(mutmatch)

# Execute the function to get selection estimates using neighboring genes
selection_estimates <- get_selection_estimates_neighbors(
  hgnc = "KRAS",
  mutationsPath = system.file("extdata", "example_mutations.csv.gz", package = "mutmatch"),
  annotationGenePath = system.file("extdata", "example_gene_annotation.csv.gz", package = "mutmatch"),
  annotationGenomeWidePath = system.file("extdata", "example_genomewide_annotation.csv", package = "mutmatch"),
  neighborsWindow = "0.5Mb",
  outlierNeighborsThreshold = 0.2
)
```

### Estimating Selection Using Low-CADD Regions

#### Step 1: Download CADD Scores

First, download the Combined Annotation-Dependent Depletion (CADD) scores. Note that these scores are not included in the package due to their large file size.

##### For UNIX-like Systems (Linux, macOS)

Replace `"your/destination/path/CADD_GRCh37-v1.4.bw"` with the path where you want to save the file.

```R
# Load the mutmatch library
library(mutmatch)

# Specify the path to store the CADD file
caddScoresPath = "your/destination/path/CADD_GRCh37-v1.4.bw"

# Download the CADD file
download_cadd_file(caddScoresPath = caddScoresPath)
```

##### For Windows Systems

For Windows machines, download the CADD scores directly from [this link](https://krishna.gs.washington.edu/download/CADD/bigWig/CADD_GRCh37-v1.4.bw).

#### Step 2: Run the Command

After obtaining the CADD scores, execute the following command to estimate selection using low-CADD regions as a neutral mutational rate baseline.

```R
# Execute the function to get selection estimates using low-CADD regions
selection_estimates <- get_selection_estimates_cadd(
  hgnc = "KRAS",
  caddScoresPath = caddScoresPath,
  mutationsPath = system.file("extdata", "example_mutations.csv.gz", package = "mutmatch"),
  annotationGenePath = system.file("extdata", "example_gene_annotation.csv.gz", package = "mutmatch"),
  annotationGenomeWidePath = system.file("extdata", "example_genomewide_annotation.csv", package = "mutmatch")
)
```

## Example Data

The package includes example datasets for convenient testing and demonstration:

- `example_mutations.csv.gz`: Example mutations dataset.
- `example_gene_annotation.csv.gz`: Example gene annotation dataset.
- `example_genomewide_annotation.csv`: Example genome-wide annotation dataset.

These files can be found under `inst/extdata/` and accessed using the `system.file()` function. 
For example:

```R
genomewide_annotation = data.table::fread(system.file("extdata", 
"example_genomewide_annotation.csv", package = "mutmatch"))
```

## License

This package is released under the MIT License. For more details, see the [`LICENSE.md`](LICENSE.md) file.


