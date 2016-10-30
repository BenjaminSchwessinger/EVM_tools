import argparse
import pandas as pd
import re

# the warnings are switched off. These are raised at line 36 and 39 for
# how the columns of the df is renamed
import warnings
warnings.filterwarnings('ignore')

parser = argparse.ArgumentParser(
    description="""
        This converts a Coding_Quarry GFF3 file to a format compatible with Evidence modler.
        It requires the following modules.
        pandas, re, argparse.
         """)

parser.add_argument('CQ_in_file', help="Name of CQ_gff3_file", type=str)

args = parser.parse_args()

CQ_file = args.CQ_in_file


# here starts the main script reading in the file using pandas and making
# a copy of it
def add_column_9(x):
    """Function to add a 9th column to assure assocation across different contigs.
    Needs to be deleted at the end again."""
    match = re.search('(?<=ID=)\w+;', x.replace('.', '_').replace('CDS:', '')
                      )  # remeber that '.' are not part of the \w+ characters
    gene_id = match.group(0)[:-1]
    return gene_id

CQ_df = pd.read_csv(CQ_file, header=None, sep="\t")
EVM_df = CQ_df.iloc[:, :]
EVM_df[9] = EVM_df[8].apply(add_column_9)

gdf = EVM_df[EVM_df[2] == 'gene']  # add the exon part for EVM df
gdf[2] = 'fmRNA'

edf = EVM_df[EVM_df[2] == 'CDS']
edf[2] = 'exon'

EVM_df = EVM_df.append(gdf).append(edf)  # make the EVM dataframe
# EVM_df.loc[:,7] ='.' #first move to make columne 7 a '.' instead of frame
EVM_df.loc[:, 1] = 'CodingQuarry_v2'

# now after we generated the overall EVM dataframe we need to bring the house in order and sort correctly
# before we iterate over the each row to adapt the format

# adding row #9 ensures that genes that belong together
EVM_df.sort_values(by=[0, 9, 3, 2], ascending=[0, 0, 1, 0], inplace=True)
# stay together
EVM_df.reset_index(drop=True, inplace=True)

name = 'Name=CodingQuarry_v2_prediction'

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

outname = CQ_file.split('.')[0] + '.EVM.gff3'

# now safe it and move on

EVM_df.iloc[:, 0:9].to_csv(outname, sep='\t', header=None, index=None)

exit()
