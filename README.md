# detwin


[![license](https://img.shields.io/badge/License-GPL--v3-blue)](https://www.gnu.org/licenses/gpl-3.0.en.html)
![language](https://img.shields.io/badge/Language-C-blue)


**detwin** is a third-party patch for **CrystFEL software**, which helps solve twinning problems in SFX crystal indexing.

[Method paper](http://journals.iucr.org/m/issues/2014/06/00/it5003/index.html) ( Liu, H., & Spence, J. C. (2014). The indexing ambiguity in serial femtosecond crystallography (SFX) resolved using an expectation maximization algorithm. IUCrJ, 1(6), 393-401. )

[Software paper](https://www.mdpi.com/2073-4352/10/7/588)(Shi, Y. , & Liu, H. . (2020). Em-detwin: a program for resolving indexing ambiguity in serial crystallography using the expectation-maximization algorithm. Crystals, 10(7), 588.)

### - How to install it ?

**First**, download CrystFEL package [here](http://www.desy.de/~twhite/crystfel/download.html). Unzip it.

**Then**, download detwin package from this github repo, unzip it, and go into the folder.

**Next**, use command `./replaceFile.sh <crystfel-folder-path>` to merge detwin source files into CrystFEL package. For example:

```
Macintosh:~ detwin$ ./replaceFile.sh ../crystfel-0.9.0
```

**Finally**, install CrystFEL !

### - How to use it ?

**I** . Use unix *man* command for help information.

**II**. Refer to command line help information : `detwin -h`.

```
Syntax: detwin [options]

read from stream file and convert the data to miller indexed FEL Bragg intensities.
output all crystals HKL, even numbered crystals HKL, odd numbered crystals HKL

  -h, --help                Display this help message.
      --version             Print CrystFEL version number and exit.
  -i, --input=<filename>    Specify input filename ("-" for stdin).
  -o, --output=<filename>   Specify output filename for merged intensities
                             Default: processed.hkl).
      --stat=<filename>     Specify output filename for merging statistics.
  -y, --symmetry=<sym>      Merge according to point group <sym>.
  -k, --spacegroupNum=<k>   Specify space group number (143~199).
  -m, --max-niter=<t>       Number of iterations for de-twinning, default is 30

  -s  --start-after=<n>     Skip <n> crystals at the start of the stream.
  -f  --stop-after=<n>      Stop after merging <n> crystals.
  -g, --histogram=<h,k,l>   Calculate the histogram of measurements for this
                             reflection.
  -z, --hist-parameters     Set the range for the histogram and the number of
          =<min,max,nbins>   bins. 

      --scale               Scale each pattern for best fit with the current
                             model. Operated before detwining.
      --winner-takes-all    If set, among all possible twinning modes, only the one 
                             with higest cc to reference model will merge.
                             Default: all twinning modes merge in a weighted way
      --highres=<n>         Reject reflections with resolution (A) higher than n 
                            while calculating CC in detwinning, default is Inf.
      --lowres=<n>          Reject reflections with resolution (A) lower than n 
                            while calculating CC in detwinning, default is 0.
      --no-polarisation     Disable polarisation correction.
      --cc-only             Only calculate and display CC between -i and -o hkl file
      --min-measurements=<n> Require at least <n> measurements before a
                             reflection appears in the output.  Default: 2
      --min-snr=<n>         Require individual intensity measurements to
                             have I > n * sigma(I).  Default: -infinity.
      --min-cc=<n>          Reject frames with CC less than n. Default: infinity.
      --max-adu=<n>         Maximum peak value.  Default: infinity.
      --min-res=<n>         Merge only crystals which diffract above <n> A.
      --push-res=<n>        Integrate higher than apparent resolution cutoff.
      --write-assignments   Write reindexed results of the crystals in original
              =<filename>   file to filename.
      (only packages with version >= 0.8.0 have the following options)
      --add-operators=<op>  Add twinning operators manually, such as '--add-operators=
                            -h,-k,-l/k,h,l'. DO NOT add space between characters !
      --write-stream=<fn>   Write reindexed stream file to fn.
```

### - Supported CrystFEL version

[![license](https://img.shields.io/badge/built-crystfel--0.6.3-green)](http://www.desy.de/~twhite/crystfel/crystfel-0.6.3.tar.gz)
[![license](https://img.shields.io/badge/built-crystfel--0.7.0-yellow)](http://www.desy.de/~twhite/crystfel/crystfel-0.7.0.tar.gz)
[![license](https://img.shields.io/badge/built-crystfel--0.8.0-red)](http://www.desy.de/~twhite/crystfel/crystfel-0.8.0.tar.gz)
[![license](https://img.shields.io/badge/built-crystfel--0.9.0-blue)](http://www.desy.de/~twhite/crystfel/crystfel-0.9.0.tar.gz)




