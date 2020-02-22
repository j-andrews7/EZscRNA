# EZscRNA

[![Travis build status](https://travis-ci.org/j-andrews7/EZscRNA.svg?branch=master)](https://travis-ci.org/j-andrews7/EZscRNA)
[![AppVeyor build status](https://ci.appveyor.com/api/projects/status/github/j-andrews7/EZscRNA?branch=master&svg=true)](https://ci.appveyor.com/project/j-andrews7/EZscRNA)

This R package is meant to be a high-level wrapper around well-established scRNA R packages (e.g. Seurat, scran, slingshot), created in an effort to streamline the end-to-end process for new users while still maintaining the ability to customize the analyses fully for more advanced users. It also attempts to include convenience functions for integrating multimodal data that other tools don't seem to address (mainly VDJ data) and cell type inference using any supplied reference dataset.

Few changes should be necessary for basic users.

It was developed for more convenient maintenance of our workflows, but if anyone else finds it useful, all the better.

**Under active development, things are/will be broken.**

(Better) documentation to come.

## Installation
Currently, EZscRNA can be installed with devtools:
```R
devtools::install_github("j-andrews7/EZscRNA")
```

It will be added to CRAN once stable and more well-documented.

## Quick Start

To be written.

## Usage
To be written.

## Issues

Bug reports, issues, and suggestions are welcome and can be submitted to the [issue tracker](https://github.com/j-andrews7/EZscRNA/issues).


## Contributing

Pull requests fixing bugs or adding features are welcome. Please follow the [Google R style guide](https://google.github.io/styleguide/Rguide.xml) and properly document your code in [roxygen format](http://r-pkgs.had.co.nz/man.html).

## License
This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program.  If not, see <https://www.gnu.org/licenses/>.