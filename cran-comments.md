This is the second submission of this package to CRAN. 
The package was removed from CRAN on 2021-10-21 after build error occured with Solaris.

It is being resubmitted with minor changes: 
- 
- In vignettes/helmet_data_example.Rmd, I removed the blocks of code saving plots to disk.

## Test environments
Tested on the 4 environments used by `rhub::check_for_cran()`.
* Debian Linux, R-devel, GCC ASAN/UBSAN
* Fedora Linux, R-devel, clang, gfortran
* Ubuntu Linux 20.04.1 LTS, R-release, GCC
* Windows Server 2022, R-devel, 64 bit

## R CMD check results
Testing returned up to two notes, which I reckon can be ignored.

# Debian Linux, R-devel, GCC ASAN/UBSAN
No ERRORs, no WARNINGs and no NOTES.

# Fedora Linux, R-devel, clang, gfortran
One NOTE:
* checking CRAN incoming feasibility ... NOTE
Maintainer: 'Vivien Goepp <vivien.goepp@gmail.com>'
New submission
Package was archived on CRAN

# Ubuntu Linux 20.04.1 LTS, R-release, GCC
One NOTE:
* checking CRAN incoming feasibility ... NOTE
Maintainer: 'Vivien Goepp <vivien.goepp@gmail.com>'

New submission
Package was archived on CRAN

# Windows Server 2022, R-devel, 64 bit
Two NOTES:
* checking CRAN incoming feasibility ... NOTE
Maintainer: 'Vivien Goepp <vivien.goepp@gmail.com>'
New submission
Package was archived on CRAN

* checking for detritus in the temp directory ... NOTE
Found the following files/directories:
  'lastMiKTeXException'
  
