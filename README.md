This is the repository storing the text and reproducible data processing for 'Decomposing cerebral blood flow MRI into functional and structural components:  A non-local approach based on prediction'.  To successfully run the code, you need to have the following dependencies: 
* [R](www.r-cran.org)
* [ANTsR](https://github.com/stnava/ANTsR.git)
* [PatchAnalysis](https://github.com/bkandel/PatchAnalysis.git) 

Put PatchAnalysis in a place where a system call from R will be able to find it.  Building PatchAnalysis requires [CMake](http://www.cmake.org/cmake/resources/software.html) and [ITK](http://www.itk.org/ITK/resources/software.html).

To compile, run `./Rcompile` in this directory.    
