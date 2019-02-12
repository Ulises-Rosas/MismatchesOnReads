# MismatchesOnReads

Number of mismatches between reads and adapters can be graphically assessed by using regression analysis. This specific parameter from your dataset can also be coupled to Trimmomatic software. 


## Usage

Scripts presented in this repository run together in a single directory and with `*.fastq` files:

```Bash
python  curves_adapters.py *.fastq
Rscript ouput_ploter.R
```
The assessment of other target sequences is also allowed by including the `-adapt` argument in the python script.

## Output

Different files are created from the assessment of conventional indexes and adapters commonly used in an Illumina sequecing event by default

```Bash
$ ls -t .
TruSeq_Universal_Reversed.png  TruSeq_Index.png           TruSeq_Universal_Reversed.csv  template.sh.e883561  curves_adapters.py
TruSeq_Universal.png           TruSeq_Index_Reversed.csv  TruSeq_Universal.csv           template.sh.o883561
TruSeq_Index_Reversed.png      TruSeq_Index.csv           ouput_ploter.R                 template.sh
```
The following image is one of the prior results:

![](https://github.com/Ulises-Rosas/MismatchesOnReads/blob/master/TruSeq_Index.png)
