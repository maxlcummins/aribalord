import pandas as pd
import numpy as np
import regex as re
import warnings
import argparse
import glob
import os
from functools import reduce

parser = argparse.ArgumentParser(description='Process and join all genotypic and phylogenetic ARIBA output')
parser.add_argument('input_dir', help='Input directory containing csv files')
parser.add_argument('output_file', help='Output filename')
parser.add_argument('--clean', action='store_true', help='Use to clean anything after _R1 from filenames')
parser.add_argument('-chars','--characters', action='store_true', help='Pattern to trim after')
args = parser.parse_args()

pattern = '*.csv'

csv_files = glob.glob(pattern)

dflist = []

names = []

for csv in csv_files:

    df = pd.read_csv(csv)


    if args.clean:
        df['name'] = df['name'].replace("{}.*".format(args.characters), "", regex = True)
    else:
        df['name'] = df['name'].replace("_R1.*", "", regex = True)

    df.set_index('name',inplace=True)

    #strip whitespace from column heads
    df_clean = df.rename(columns=lambda x: x.strip())

    #convert string values to representative integers
    cats = ['yes', 'no', 'yes_nonunique', 'partial', 'interrupted', 'fragmented']
    values = [1, 0, 1, 0, 0, 0]
    df_test = df_clean.replace(to_replace = cats, value = values)

    #initialise empty list for assembled,ref_seq dataframes
    df_list = []

    #define number of columns, convert to a list, iterate this list and bundle into batches of two
    new_df_length = int(len(df.columns))
    split = list(range(new_df_length))
    it = iter(split)
    zip_it = zip(it, it)

    #for the two columns associated with each gene hit
    #save these columns and their index to a data frame
    #then append this dataframe to a list of dataframes
    for i in zip_it:
        df_pair = df_test.iloc[:,[i[0],i[1]]]
        df_list.append(df_pair)

    #initialise empty list
    frame = []    

    #for dataframe in df_list
    #pivot the tables based so gene hits are colnames and 'assembled' status are observations
    #remove trash column that is generated in column 0
    #change NaNs to zeroes
    #append these modified dataframes to new list of dataframes called frame
    for i in df_list :
        data_pivot = (i.pivot(values = i.columns[0], columns = i.columns[1]))
        data_pivot = data_pivot.iloc[:,1:]
        data_pivot = data_pivot.fillna(0)
        frame.append(data_pivot)

    #make list of original sample names using index of original df    
    base = df.iloc[:,0:0]

    #join new dataframes from frames to base dataframe
    data_all = base.join(frame)
    
    dflist.append(data_all)

    #combine base names for a merge later
    names.append(base)



#concatenate together the names list
namelist = pd.concat(names)

#reset the index
namelist.reset_index(inplace=True)

for csv in dflist:
    
    #print(csv.head())
    csv.reset_index(inplace=True)
    namelist = pd.merge(namelist,csv,how='outer')

full_dedupe = namelist.drop_duplicates()

namelist = pd.concat(names)

#reset the index
namelist.reset_index(inplace=True)

dflist2 = []

#regular expressions to clean column names
for df in dflist:
    
    #patterns used to recognise spreadsheets as being of serotype, custom, virulence etc
    #to apply the correct regexs to clean column names
    phylogroup = [(re.search(r"^arpA|^tspE4\.C2|^yjaA", i)) for i in df]
    virulence = [(re.search(r"^gad", i)) for i in df]
    resistance = [(re.search(r"^dfrA|^tet|^bla|^sul|^aadA", i)) for i in df]
    plasmid = [(re.search(r"^Inc", i)) for i in df]
    insertion = [(re.search(r"^IS", i)) for i in df]
    custom = [(re.search(r"^fimH", i)) for i in df]
    serotype = [(re.search(r"^wzy|^wzx|^wzt|^wzm|^fliC", i)) for i in df]
    MLST = [(re.search(r"^adk|^fumC|^gyrB", i)) for i in df]
    
    if any(phylogroup):
        df = df.rename(columns=lambda x: re.sub('(.*)',r'phylogroup_\1',x))

    elif any(custom):
        df = df.rename(columns=lambda x: re.sub('(.*)(_|\.)[0-9]+_.*',r'\1',x))
        df = df.rename(columns=lambda x: re.sub('^(.*)(NC_)?_.*',r'v_\1',x))
        df = df.rename(columns=lambda x: re.sub('ipaH','v_ipaH',x))
        df = df.rename(columns=lambda x: re.sub('^v_IS','i_IS',x))
        df = df.rename(columns=lambda x: re.sub('^v_int','i_int',x))
        df = df.rename(columns=lambda x: re.sub('^v_(mer|ter|sil|pco|czc|scs)',r'r_\1',x))
        df = df.rename(columns=lambda x: re.sub('_$','',x))
        df = df.rename(columns=lambda x: re.sub('_NC$','',x))
        df = df.rename(columns=lambda x: re.sub('v_kpsMT_II__K2_CP000468.1__APEC','v_kpsMT_II',x))
        df = df.rename(columns=lambda x: re.sub('v_kpsMT_II__K1_AE014075.1','v_kpsMT_II',x))
        df = df.rename(columns=lambda x: re.sub('v_kpsMT_III__K54','v_kpsMT_III',x))
        df = df.rename(columns=lambda x: re.sub('v_kpsMT_II__K5_CP002212.1_clone_D','v_kpsMT_II',x))

    elif any(virulence):
        df = df.rename(columns=lambda x: re.sub('(.*)(_|\.)[0-9]+_.*',r'\1',x))
        df = df.rename(columns=lambda x: re.sub('^(.*)',r'v_\1',x))

    elif any(serotype):
        df = df.rename(columns=lambda x: re.sub('^([^_]+).*(_.*)$',r'sero_\1\2',x))

    elif any(resistance):
        df = df.rename(columns=lambda x: re.sub('(.*)(_|\.)[0-9]+_.*',r'\1',x))
        df = df.rename(columns=lambda x: re.sub('^(.*)',r'r_\1',x))

    elif any(plasmid):
        df = df.rename(columns=lambda x: re.sub('(.*)(_|\.)[0-9]+_.*',r'\1',x))
        df = df.rename(columns=lambda x: re.sub('^FIA',r'IncFIA',x))
        df = df.rename(columns=lambda x: re.sub('^(.*)',r'p_\1',x))
        df = df.rename(columns=lambda x: re.sub('(p_IncF[A-Z]+)_.*',r'\1',x))

    elif any(insertion):
        df = df.rename(columns=lambda x: re.sub('([^-])_.*',r'\1',x))
        df = df.rename(columns=lambda x: re.sub('^(IS.*)',r'i_\1',x))

    elif any(MLST):
        warnings.warn('CSV file contains MLST alleles.', Warning)
    else:
        warnings.warn('CSV file contents could not be identified and therefore may have sloppy column headers.', Warning)
    
    df = df.rename(columns=lambda x: re.sub('.*name','name',x))
    dflist2.append(df)

    
#merge dataframes
simple_dedupe = reduce(lambda x, y: pd.merge(x, y, on = 'name', how='outer'), dflist2)

#set index to 'name' column
df.set_index('name',inplace=True)

#sort columns by their names
simple_dedupe = simple_dedupe.sort_index(axis=1)

#collapse columns that share a column name
simple_dedupe = simple_dedupe.groupby(simple_dedupe.columns, axis=1).sum()
full_dedupe = full_dedupe.groupby(full_dedupe.columns, axis=1).sum()

#remove columns with colsums of zero
simple_dedupe = simple_dedupe.loc[:, (simple_dedupe != 0).any(axis=0)]
full_dedupe = full_dedupe.loc[:, (full_dedupe != 0).any(axis=0)]

#Move 'name' to start of the dataframe
simple_dedupe = simple_dedupe.reindex(columns=(['name'] + list([a for a in simple_dedupe.columns if a != 'name']) ))

simple_dedupe.to_csv(args.output_file+'simple_binary.csv',index=False)

simple_dedupe.set_index('name', inplace=True)
simple_dedupe[simple_dedupe > 1] = 1

#replace all >1 hits with 1

simple_dedupe.reset_index(inplace=True)

#Write to CSV
simple_dedupe.to_csv(args.output_file+'simple_multi.csv',index=False)
full_dedupe.to_csv(args.output_file+'full.csv',index=False)