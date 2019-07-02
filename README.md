# EZscRNA

This R package is meant to be a high-level wrapper around popular, well-maintained scRNA R packages (Seurat, Slingshot, Monocle), created in an effort to streamline the end-to-end process for new users while still maintaining the ability to customize the analyses fully for more advanced users. It also attempts to include convenience functions for integrating multimodal data that other tools don't seem to address (mainly VDJ data) and cell type inference using any supplied reference dataset.

It was written with best practices as defined by the creators of Seurat, Monocle, and scanpy in mind, so few changes should be necessary for basic users.

It was developed for more convenient maintenance of our workflows, but if anyone else finds it useful, all the better.

**Under active development, things are/will be broken.**

(Better) documentation to come.

## Installation
Currently, EZscRNA can be installed with devtools:
```R
devtools::install_github("j-andrews7/EZscRNA")
```

I'll get it on CRAN once I have a stable version.

## Quick Start

To be written.

## Usage
To be written.

## Issues

Bug reports, issues, and suggestions are welcome and can be submitted to the [issue tracker](https://github.com/j-andrews7/EZscRNA/issues).


## Contributing

Pull requests fixing bugs or adding features are welcome. Please follow the [Google R style guide](https://google.github.io/styleguide/Rguide.xml) and properly document your code in [roxygen format](http://r-pkgs.had.co.nz/man.html).