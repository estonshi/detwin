# detwin


[![license](https://img.shields.io/badge/License-GPL--v3-blue)](https://www.gnu.org/licenses/gpl-3.0.en.html)
![language](https://img.shields.io/badge/Language-C-blue)


**detwin** is a third-party patch for **crystfel software**, which helps solving twinning problems in crystal indexing.

[Reference paper](http://journals.iucr.org/m/issues/2014/06/00/it5003/index.html) ( Liu, H., & Spence, J. C. (2014). The indexing ambiguity in serial femtosecond crystallography (SFX) resolved using an expectation maximization algorithm. IUCrJ, 1(6), 393-401. )

### - How to install it ?

**First**, download crystfel package [here](http://www.desy.de/~twhite/crystfel/download.html). Unzip it.

**Then**, download detwin from this github repo, unzip and go into the folder.

**Next**, type `./replaceFile.sh` in terminal. It will ask for your crystfel download path, provide it ! (ABSOLUTE PATH)

**Finally**, install crystfel !

### - How to use it ?

**I** . Using unix *man* command for help information.

**II**. Type `$CRYSTFEL_PATH/bin/detwin -h` and you will see:

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
```

### - Supported

[![license](https://img.shields.io/badge/built-crystfel--0.6.3-green)](http://www.desy.de/~twhite/crystfel/crystfel-0.6.3.tar.gz)
[![license](https://img.shields.io/badge/built-crystfel--0.7.0-yellow)](http://www.desy.de/~twhite/crystfel/crystfel-0.7.0.tar.gz)
[![license](https://img.shields.io/badge/built-crystfel--0.8.0-red)](http://www.desy.de/~twhite/crystfel/crystfel-0.8.0.tar.gz)




