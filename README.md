# MismatchesOnReads

Number of mismatches between reads and adapters can be graphically assessed by using regression analysis. This specific parameter from your dataset can also be coupled to Trimmomatic software. 


## Mode of use

Scripts presented in this repository run together in a shell. The following shell is a template for running both in a cluster:

```Bash
#!/bin/bash -l
# FILENAME: template.sh

#PBS -l nodes=1:ppn=20
#PBS -l walltime=00:10:00
#PBS -q [specify a queue here]

cd $PWD/fastqs

#load both python and r modules:

module load python/3.4.1
module load r/3.4.3

python  curves_adapters.py *.fastq

Rscript ouput_ploter.R
```

## Output

```Bash
$ ls .
checking_adapters.pbs          curves_adapters.py  TruSeq_Index.png           TruSeq_Universal.csv           TruSeq_Universal_Reversed.png
checking_adapters.pbs.e883561  ouput_ploter.R      TruSeq_Index_Reversed.csv  TruSeq_Universal.png
checking_adapters.pbs.o883561  TruSeq_Index.csv    TruSeq_Index_Reversed.png  TruSeq_Universal_Reversed.csv
```
