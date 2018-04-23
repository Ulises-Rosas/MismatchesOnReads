import re
import pandas as pd
import argparse
from scipy import stats

adapters = {'TruSeq_Universal': 'ATGATACGGCGACCACCGAG',
            'TruSeq_Universal_Reversed': 'GATCGGAAGAGCGTCGTGTA',
            'TruSeq_Index': 'GATCGGAAGAGCACACGTCT',
            'TruSeq_Index_Reversed': 'AAGCAGAAGACGGCATACGA'}

parser = argparse.ArgumentParser(description='Profile of matches between reads and \
                                 adapters')

parser.add_argument('fastq', metavar='FastQ',
                    nargs = "+",
                    type=argparse.FileType('r'),
                    help='FastQ file whose reads will be used to render input \
                    at the profiler')
parser.add_argument('-adapt', metavar='Adapters', 
                    type=argparse.FileType('r'),
                    default = adapters,
                    help='Adapters file')
args = parser.parse_args()


"""----------------------------------------------------------------------------
fas_to_dic(x) function transforms fasta files to a dictionary data type
To obtain more infomation regards the following function, please visit: 
       https://github.com/Ulises-Rosas/Fasta-to-Dictionary#fasta-to-dictionary     
-------------------------------------------------------------------------------"""

def fas_to_dic(x):       
    
    file_content = x.readlines()    
    seqs_list = [] 
    
    for i in file_content:
        seqs_list.append(i.replace("\n", "")) 
        
    keys = [] 
    values = []
    
    i = 0
    
    while(">" in seqs_list[i]):
        keys.append(seqs_list[i].replace(">", ""))
        i += 1 
        JustOneValue = []
        while((">" in seqs_list[i]) == False):
            JustOneValue.append(seqs_list[i]) 
            i += 1
            if(i == len(seqs_list)):
                i -= 1
                break
        values.append("".join(JustOneValue))
    return dict(zip(keys, values))


#print(args.fastq.readline())

if type(args.adapt) == dict:
       adapters_to_test =  args.adapt
       #print([i[1] for i in args.adapt.items()][0])
else: 
       adapters_to_test = fas_to_dic(args.adapt)
       #print([i[1] for i in fas_to_dic(args.adapt).items()][0])


def seqs_extracter(fastqs):
    
    results = []   
       
    for fastq in fastqs:
           all_reads = fastq.readlines()
           all_seqs = []
           i = 1
           
           while( i < len(all_reads)):
                  all_seqs.append(all_reads[i].replace("\n",""))
                  i +=4 
                  
           results.append([fastq.name[ fastq.name.rfind('/') + 1: ].replace(".fastq", ""), all_seqs])
           
    return results
	
	
""" 
This function retrive the number of matches given a certain and static sequence
"""

def seqs_matches(adapter_seq, all_seqs):
    matches = []    
    for i in range(len(all_seqs)):           
        if len(re.findall(r"{}".format(adapter_seq) + r"{1}", all_seqs[i])) > 0:               
            ##retrieve position (that's why you use range() function)\
            ##with non-empty elements in lists   
            matches.append(i)            
    return matches

""" ---------------------------------------------------------------------------
this function searches for matches by increasing the word size of the adapter in all samples
their input format are items from adapters and the parser argument directly retrived from comand line.
Outcome from this function is a single table which will be used in the next one function 
-------------------------------------------------------------------------------"""


def sliding_adapter(adapter_items, fastqs_read):
       
    all_seqs  = seqs_extracter(fastqs_read)
       
    all_samples_matches =  pd.DataFrame()
    
    for seq in all_seqs:
           
           nrow = max(  [ len(i[1]) for i in adapter_items   ] )
           
           ##series stored in a list
           storing_values = [
                         pd.Series([seq[0]]*nrow, name = "Sample"),
                         pd.Series([i for i in range(1, nrow + 1)], name = "Position")]
           
           for keys, values in adapter_items:
                  
                  index = range(1, len(values) + 1) 
                   
                  ##generator of series by adapters
                  matches = [ len(seqs_matches(values[:i], seq[1])) for i in index ]
                  
                  storing_values.append( pd.Series( data = matches, name = keys ) )
                 
           all_samples_matches = all_samples_matches.append( pd.concat( storing_values, axis = 1 ), ##bind columns
                                                              ignore_index = True  ##then, bind rows
                                                              )
    return all_samples_matches
	
	
""" ---------------------------------------------------------------------------
this function uses outcomes from the previous one and furnishes summaries of regressions analyses.
The output has the following structure: 
       
       Sample    Bases_from_end    Slope    Intercept    R-squared   p-value
       ...                    2      ...          ...          ...        ~1

Regression analysis starts from the last two position of the adapter, requiered for tracing a line.
Futher trimming is needed as some linear model are non-significant. Therfore, a data frame filtration in function
of the  column "p-value" will be usefull
-------------------------------------------------------------------------------"""


def sliding_regression(slinding_adapter_output, adapter_key):
       
    samples = [ i for i in set(slinding_adapter_output["Sample"]) ]

    results = pd.DataFrame()

    for one_sample in samples:
        
        sample_indexer =  slinding_adapter_output["Sample"] == one_sample   
        
        trimmed_by_sample = slinding_adapter_output[sample_indexer]
        
        slicing_index = range( 2, max(trimmed_by_sample["Position"] + 1) )

        sample = pd.Series([ one_sample]*len(slicing_index), name="Sample" )
        position = pd.Series([ i for i in slicing_index ], name= "Bases_from_end")
        slope_one = []
        intercept_one = []
        rvalue_one = []
        pvalue_one = [] 

        for i in slicing_index:
               
            indexer = trimmed_by_sample["Position"] >= min(position[-i:])
            
            trimmed_by_position = trimmed_by_sample[indexer]
            ##axis for linear regression:
            a = trimmed_by_position["Position"]
            b = trimmed_by_position[adapter_key] ##this is an unique key and is pro-
                                                 ##vided in arguments of the function
            ##linear model:
            slope, intercept, rvalue, pvalue, std = stats.linregress(a,b)
            
            slope_one.append(slope)
            intercept_one.append(intercept)
            rvalue_one.append(rvalue**2)
            pvalue_one.append(pvalue)

        tmp_df = pd.concat([ sample, position,
            pd.Series(slope_one, name="Slope"),
            pd.Series(intercept_one, name="Intercept"),
            pd.Series(rvalue_one, name= "R-squared"),
            pd.Series(pvalue_one, name="p-value")], axis= 1)
                  
        results = results.append(tmp_df, ignore_index= True)
        
    return results
       	   

""" ---------------------------------------------------------------------------
Multiple iteration are done for buiding up *.csv by adapters
-------------------------------------------------------------------------------"""

output = sliding_adapter(adapters_to_test.items(), args.fastq)


for keys in adapters_to_test.keys():
    for_csv = sliding_regression(output, keys)
    for_csv.to_csv(keys + ".csv", header = True, index = False)
