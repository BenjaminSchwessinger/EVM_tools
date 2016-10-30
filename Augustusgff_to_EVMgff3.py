import argparse
import pandas as pd
import re

#the warnings are switched off. These are raised at line 36 and 39 for how the columns of the df is renamed
import warnings 
warnings.filterwarnings('ignore')

parser = argparse.ArgumentParser(
        description="""
        This converts a Augustus GFF file to a format compatible with Evidence modler.
        It requires the following modules.
        pandas, re, argparse.
         """)

parser.add_argument('Augustus_in_file', help="Name of Augustus_gff_file", type=str)

args = parser.parse_args()

Augustus_file = args.Augustus_in_file

def add_column_9(x):
    """Function to add a 9th column to assure assocation across different contigs.
    Needs to be deleted at the end again."""
    match = re.search('g[0-9]+', x) 
    gene_id = match.group(0)
    return gene_id

def fix_column_9(x):
    """Function to fix 9th column to assure assocation across different alternative transcripts.
    Needs to be deleted at the end again."""
    if re.search('g[0-9]+.t[2-9][0-9]?', x) :
        gene_id = re.search('g[0-9]+.t[2-9][0-9]?', x).group(0)
        gene_id = gene_id.replace('.t','_')
    else:
        match = re.search('g[0-9]+', x) 
        gene_id = match.group(0)
    return gene_id


Augustus_df = pd.read_csv(Augustus_file, sep='\t', header=None)

EVM_df = Augustus_df[(Augustus_df[2]=='gene')|(Augustus_df[2]=='transcript')|(Augustus_df[2]=='CDS')]
EVM_df[9] = EVM_df[8].apply(add_column_9) #this binds everything together and enables appropriate sorting later on

#add the gene part for EVM df when t2+
gdf = EVM_df[(EVM_df[8].str.contains('t[2-9][0-9]?'))&(EVM_df[2]=='transcript')] 
gdf.loc[:,2] = 'gene'

edf= EVM_df[EVM_df[2]=='CDS'] #add the exon part for EVM df
edf.loc[:,2]= 'exon'

EVM_df = EVM_df.append(gdf).append(edf)

#EVM_df.loc[:,7] ='.' #first move to make columne 7 a '.' instead of frame
EVM_df.loc[:,5] ='.' #first move to make columne 5 a '.' probability value from Augustus
EVM_df.loc[:,1] ='Augustus' # to make columne 2 a '.' instead of frame

EVM_df.replace(to_replace='transcript', value='fmRNA', inplace=True)

EVM_df[9] = EVM_df[8].apply(fix_column_9)

EVM_df.sort_values(by=[0,9,3,2],ascending=[0,0,1,0],inplace=True) #adding row #9 ensures that genes that belong together 
#stay together
EVM_df.reset_index(drop=True,inplace=True)

name = 'Name=braker1_Augustus_prediction'
#go over the whole data frame and change the nameing in columne 8 according to EVM syntax
for x,y in EVM_df.iterrows():
    if y[2] =='gene':
        counter = 1
        gene_id = y[0]+'-'+y[9]
        EVM_df.iloc[x,8] = 'ID='+gene_id+';'+name   
    if y[2] == 'fmRNA':
        EVM_df.iloc[x,8] = 'ID='+gene_id+'.t1;Parent='+gene_id+';'+name
        EVM_df.iloc[x,2] = 'mRNA'
    if y[2] == 'exon':
        EVM_df.iloc[x,8] = 'ID='+gene_id+'.t1.exon'+str(counter)+';Parent='+gene_id+'.t1'
        counter += 1
    if y[2] == 'CDS':
        EVM_df.iloc[x,8] = 'ID=cds.'+gene_id+'.t1;Parent='+gene_id+'.t1'

outname = Augustus_file.split('.')[0]+'.augustus.EVM.gff3'

#now safe it and move on

EVM_df.iloc[:,0:9].to_csv(outname, sep='\t', header=None, index=None)
