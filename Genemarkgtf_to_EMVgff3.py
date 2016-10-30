import argparse
import pandas as pd
import re

# the warnings are switched off. These are raised at line 36 and 39 for
# how the columns of the df is renamed
import warnings
warnings.filterwarnings('ignore')

parser = argparse.ArgumentParser(
    description="""
        This converts a Genemark GTF file to a format compatible with Evidence modler.
        Use with x.f.food.gtf file from braker1 output.
        It requires the following modules.
        pandas, re, argparse.
         """)

parser.add_argument('GM_in_file', help="Name of Genemark_GTF_file", type=str)

args = parser.parse_args()

GM_file = args.GM_in_file


# here starts the main script reading in the file using pandas and making
# a copy of it
def add_column_9(x):
    """Function to add a 9th column to assure assocation across different contigs.
    Needs to be deleted at the end again."""
    match = re.search('(?<=gene_id ")\w+', x) #remeber that '.' are not part of the \w+ characters
    gene_id = match.group(0)
    return gene_id

GM_df = pd.read_csv(GM_file, header=None, sep="\t")
EVM_df = GM_df.iloc[:, :]
EVM_df[9] = EVM_df[8].apply(add_column_9)

set_of_genes = set(EVM_df[9].tolist())

#quick test if the number of start, stop codon and number of genes add up
if len(set_of_genes) == len(set(EVM_df[EVM_df[2] =='stop_codon'][9].tolist() + EVM_df[EVM_df[2] =='start_codon'][9].tolist())):
    print('There are as many genes as start and stop codons\n All good to proceed\n')
else:
    print('Something wrong with the number of gene vs. stop codpns vs. start codons!\n please check input gtf!\n')
    exit()

for x in set_of_genes:
    #print(EVM_df[(EVM_df[9] == x)][6])
    start_side = ''
    stop_side =''
    if EVM_df[(EVM_df[9] == x)][6].tolist()[0] == '+':
        #get the start and stop sides depending on the orintation
        start_side = int(EVM_df[(EVM_df[9] == x)&(EVM_df[2]=='start_codon')][3])
        stop_side = int(EVM_df[(EVM_df[9] == x)&(EVM_df[2]=='stop_codon')][4])
        start_index = EVM_df[(EVM_df[9] == x) & (EVM_df[2] == 'start_codon')].index
        stop_index = EVM_df[(EVM_df[9] == x) & (EVM_df[2] == 'stop_codon')].index
        #convert the start and stop codons into gene and fmRNA
        #print(start_index)
        #print(start_side, stop_side)
        EVM_df.loc[start_index, 2] = 'gene'
        EVM_df.loc[start_index, 3] = start_side
        EVM_df.loc[start_index, 4] = stop_side
        EVM_df.loc[stop_index, 2] = 'fmRNA'
        EVM_df.loc[stop_index, 3] = start_side
        EVM_df.loc[stop_index, 4] = stop_side
    if EVM_df[(EVM_df[9] == x)][6].tolist()[0] == '-':
        start_side = int(EVM_df[(EVM_df[9] == x)&(EVM_df[2]=='start_codon')][4])
        stop_side = int(EVM_df[(EVM_df[9] == x)&(EVM_df[2]=='stop_codon')][3])
        start_index = EVM_df[(EVM_df[9] == x) & (EVM_df[2] == 'start_codon')].index
        stop_index = EVM_df[(EVM_df[9] == x) & (EVM_df[2] == 'stop_codon')].index
        #convert the start and stop codons into gene and fmRNA
        #print(start_index)
        #print(start_side, stop_side)
        EVM_df.loc[start_index, 2] = 'gene'
        EVM_df.loc[start_index, 3] = stop_side
        EVM_df.loc[start_index, 4] = start_side
        EVM_df.loc[stop_index, 2] = 'fmRNA'
        EVM_df.loc[stop_index, 3] = stop_side
        EVM_df.loc[stop_index, 4] = start_side



# EVM_df.loc[:,7] ='.' #first move to make columne 7 a '.' instead of frame
EVM_df.loc[:, 1] = 'GeneMark.hmm'

# now after we generated the overall EVM dataframe we need to bring the house in order and sort correctly
# before we iterate over the each row to adapt the format

# adding row #9 ensures that genes that belong together
EVM_df.sort_values(by=[0, 9, 3, 2], ascending=[0, 0, 1, 0], inplace=True)
# stay together
EVM_df.reset_index(drop=True, inplace=True)

name = 'GeneMark.hmm'
# go over the whole data frame and change the nameing in columne 8
# according to EVM syntax
for x, y in EVM_df.iterrows():
    if y[2] == 'gene':
        gene_id = ''
        counter = 1
        gene_id = y[0] + '-' + y[9]
        EVM_df.iloc[x, 8] = 'ID=' + gene_id + ';' + name + '_' + gene_id
    if y[2] == 'fmRNA':
        EVM_df.iloc[x, 8] = 'ID=' + gene_id + '.t1;Parent=' + \
            gene_id + ';' + name + '_' + gene_id
        EVM_df.iloc[x, 2] = 'mRNA'
    if y[2] == 'exon':
        EVM_df.iloc[x, 8] = 'ID=' + gene_id + '.t1.exon' + \
            str(counter) + ';Parent=' + gene_id + '.t1'
        counter += 1
    if y[2] == 'CDS':
        EVM_df.iloc[x, 8] = 'ID=cds.' + gene_id + \
            '.t1;Parent=' + gene_id + '.t1'

outname = GM_file[:-4]+ '.EVM.gff3'

# now safe it and move on

EVM_df.iloc[:, 0:9].to_csv(outname, sep='\t', header=None, index=None)

exit()
