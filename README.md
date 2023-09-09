# mutmatch

`mutmatch` is an R package designed to estimate somatic selection for a gene using various methods including neighbors baseline of mutation rates or low-impact mutations estimated with CADD score.

## Installation

You can install the development version of `mutmatch` from [GitHub](https://github.com/ebesedina/mutmatch) using:

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

Here is a simple example:

```R
library(mutmatch)

selection_estimates <- get_selection_estimates_neighbors(
hgnc = "KRAS",
mutationsPath = system.file("extdata",
"example_mutations.csv.gz", package = "mutmatch"),
annotationGenePath = system.file("extdata", "example_gene_annotation.csv.gz",
package = "mutmatch"),
annotationGenomeWidePath = system.file("extdata",
"example_genomewide_annotation.csv", package = "mutmatch"),
neighborsWindow = "0.5Mb",
outlierNeighborsThreshold = 0.2
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


