# MismatchesOnReads

Number of mismatches between reads and adapters can be graphically assessed by using regression analysis. This specific parameter from your dataset can also be coupled to Trimmomatic software. 


## Requirements

+ R (libraries: ggplot2, RColorBrewer)
+ python3 (modules: pandas & scipy)

## Usage

Scripts presented in this repository run together in a single directory and with `*.fastq` files:

```Bash
python3  curves_adapters.py *.fastq
Rscript ouput_ploter.R
```
The assessment of other target sequences is also allowed by including the `-adapt` argument in the python script.

## Output

Different files are created from the assessment of conventional indexes and adapters commonly used in an Illumina sequecing event by default

```Bash
$ ls -t .
TruSeq_Universal_Reversed.png  TruSeq_Index.png           TruSeq_Universal_Reversed.csv    curves_adapters.py
TruSeq_Universal.png           TruSeq_Index_Reversed.csv  TruSeq_Universal.csv           
TruSeq_Index_Reversed.png      TruSeq_Index.csv           ouput_ploter.R                 
```
The following image is one of the prior results:

![](https://github.com/Ulises-Rosas/MismatchesOnReads/blob/master/TruSeq_Index.png)
