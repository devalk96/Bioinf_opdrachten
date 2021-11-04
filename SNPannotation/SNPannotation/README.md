#SNPAnnotation

## Goal
The goal of this project is to create a SNP analyser.
The user is asked to provide two files. The first file contains reference protein sequences. 
While the second file contains the protein of interest. 


## Arguments
        -r --ref:    <Str> **REQUIRED** Filepath to file containing reference proteins (fasta format)
        -t --target: <Str> **REQUIRED** Filepath to file containing target proteins (fasta format)
        -o --output: <Str>              If provided, console output will be written to file

## How to:
###1. Search Homologene for reference proteins (multi fasta format)

        - Example uses: apoptosis-inducing factor 3 isoform 1 from 15 organisms
        - Link: https://www.ncbi.nlm.nih.gov/homologene/13773
        
###2. Find your target gene and download as protein file (fasta format)
    
    
    1. Click download datasets (blue button)
    2. Select Protein sequences(FASTA)
    3. Download
    
    
    - Example used: Homo Sapiens
        - Link: https://www.ncbi.nlm.nih.gov/gene/150209

###3. Run script
*Note: all files are included in the project to run underlying commands*
#### Example input 1:
    python3 run.py -r data/reference_protein_example.txt -t data/protein_of_interest_example.faa -o data/output.txt
1. The above input will get the reference proteins from file data/homologene.txt
2. The target proteins will be parsed from data/proteins.faa
3. Output will be written to data/output.txt

#### Example input 2:
        python3 run.py -r data/homologene.txt -t data/protein.faa
1. The above input will get the reference proteins from file data/homologene.txt
2. The target proteins will be parsed from data/proteins.faa
3. As parameter -o is not provided, results will only be echoed to console        
        
## Scoring explained
Score is determined by calculating the frequency of 
the specific aminoacid at a specific column index.
Scoring is scaled from 0 to 9 whereas a 9 is the 
maximum score (deletious) and 0 is lowest score (neutral)

Score is calculated using: 9 - (P * 9)
wheras P = probability of specific letter in MSA at specific column.

A snp is seen as deletious when a score is greater then 8 (p > 0.1 of occurance) 
and possible deletious at score between 6.75 (p > 0.25 of occurance) and 8.

## Contact
Email: *s.j.bouwman@st.hanze.nl*



