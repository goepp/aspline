This is the first submission of this package to CRAN.
It is being resubmitted with minor changes to: 
- bad LICENSE file specification in DESCRIPTION: changed "GPL-3 + file LICENSE" to "GPL-3"
- and the redundant presence of the GPL3 LICENSE file): added "^LICENSE$" to .Rbuildignore
- an URL missing the 'https' header: changed it in "README.Rmd" and "README.md"

## Test environments
* ubuntu 20.04, locally, R 4.0.2
* ubuntu 20.04 (on github actions), R 4.1.1 
* macOS 10.15 (on github actions), R 4.1.1
* windows (on winbuilder), R 4.1.1

## R CMD check results
There were no ERRORs, no WARNINGs and one NOTE (on windows only):

* checking CRAN incoming feasibility ... NOTE
