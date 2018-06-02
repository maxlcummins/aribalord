import pandas as pd
import numpy as np
import regex as re
import warnings
import glob
import os
from functools import reduce

#dont forgets**********



def geno(glob_path, clean, chars):
    """ Processes ARIBA data
    """

    #initialise empty list for assembled,ref_seq dataframes
    df_list = []
    dflist1 = []

    #initialise empty list to append dfs
    frame = []    

    #initialise empty list for unique names from all tables
    names = []

    for n, csv in enumerate(glob.glob('{}/FH*.csv'.format(glob_path))):

        print('\t {} found {}...'.format(n+1, csv))

        df = pd.read_csv(csv)

        if clean:
            df['name'] = df['name'].replace("{}.*".format(chars), "", regex = True)
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
    
        dflist1.append(data_all)

        names.append(base)

    return (dflist1, names)

def simple_clean(dflist):
#regular expressions to clean column names
    dflist2 = []
    
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
            df = df.rename(columns=lambda x: re.sub('(^arpA|^tspE4\.C2|^yjaA|^chuA)',r'phylogroup_\1',x))

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
            df = df.rename(columns=lambda x: re.sub('^([A-z]+)_.*_((O|H).*)$',r'\1_\2',x))
            df = df.rename(columns=lambda x: re.sub('^([A-z]+)_((O|H).*)$',r'sero_\1_\2',x))

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

        df = df.groupby(df.columns, axis=1).sum()

        dflist2.append(df)
    
    return dflist2


def mlst(glob_path):

    for n, csv in enumerate(glob.glob('{}/*.txt'.format(glob_path))):

        print('\t {} MLST file found: {}...'.format(n+1, csv))

        MLST = pd.read_csv(csv, delimiter='\t')

        MLST_pattern = [(re.search(r"^adk|^fumC|^gyrB|^icd|^mdh|purA|recA", i)) for i in MLST]

        if any(MLST_pattern) :
            MLST = MLST.rename(columns=lambda x: re.sub(r'(adk|^fumC|^gyrB|^icd|^mdh|^purA|^recA)',r'MLST_\1',x))
        else :
            warnings.warn('MLST text file does not contain expected values...', Warning)

        MLST.columns.values[-1] = 'name'

        pattern = "^ST"
        filter_ = MLST['ST'].str.contains(pattern)
        MLST = MLST[~filter_]

        mlst_table = MLST

        return mlst_table
    

def sero(table, simple_csv=True):
    
    pd.options.mode.chained_assignment = None

    table = table.reset_index()

    table = table.filter(regex=r'(^name|^sero_)', axis=1)

    table = table.rename(columns=lambda x: re.sub(r'^sero_','',x))
    
    #melt df
    table = pd.melt(table, id_vars='name')


    #filter non 1s
    table = table[table.value == 1]

    #split string of fliC_H4, for example, to two columns, one with fliC and one with H4
    table['simple'] = table.variable.str.split('_').str[0]
    table['type'] = table.variable.str.split('_').str[1]

    
    #remove unwanted columns
    table = table.iloc[:,[0,3,4]]

    #group by name and gene type ('simple') and concatenate together multiple hits for the same gene type ('simple')
    table['O_or_H'] = table.groupby(['name','simple'])['type'].transform(lambda x: '-'.join(x))

    #remove unwanted columns
    table = table.iloc[:,[0,1,3]]

    # #drop duplicate rows
    table.drop_duplicates(inplace=True)

    # #pivot table
    table = table.pivot(index = 'name', columns = 'simple', values = 'O_or_H')

    #create empty columns with NaNs - i believe this is redundant...
    table['O1cat'] = np.nan
    table['O2cat'] = np.nan
    table['O_type'] = np.nan
    table['H_type'] = np.nan

    #if no non-fliC H hits then set H_type to fliC, otherwise set to fliC and add an asterisk
    if 'flnA' or 'flmA' or 'fllA' or 'flkA' in table.columns:
        table['H_type'] = table['fliC']
    else:
        table['H_type'] = str(table['fliC']+"*")

    #replace null with blank
    table = table.where((pd.notnull(table)), str(''))

    #add two new columns combining the gene hits for wzm/wzt and wzx/wzy
    table['O1cat'] = table['wzm'].map(str) + '*/' + table['wzt'] + '*'
    table['O2cat'] = table['wzx'].map(str) + '*/' + table['wzy'] + '*'

    #if wzm==wzt then set O1 (O_type 1) to wzm, otherwise set O1 to O1cat
    table['O1'] = np.where(table['wzm']==table['wzt'], table['wzm'], table['O1cat'])

    #if wzx==wzy then set O2 (O_type 2) to wzx, otherwise set O2 to O2cat
    table['O2'] = np.where(table['wzx']==table['wzy'], table['wzx'], table['O2cat'])

    #replace */* with *
    table = table.replace("\*/\*$", "*", regex = True)

    #combined O1 and O2 with a slash between them
    table['O_type'] = table['O1'].map(str) + '/' + table['O2']

    #replace separating strings that exist in absence of gene hits
    table = table.replace("^\*-", "", regex = True)
    table = table.replace("^/", "", regex = True)
    table = table.replace("/$", "", regex = True)

    #replace non-hits with ONT
    table['O_type'] = table['O_type'].replace("^$", "ONT", regex = True)

    #combine O and H type into O:H format in column OH_type
    table['OH_type'] = table['O_type'].map(str) + ':' + table['H_type']

    table = table.reset_index()

    #make simple table of EcOH data
    EcOH = table.loc[:,['name','O_type','H_type', 'OH_type']]

    #create a more informative table
    seroframe = table
    
    if simple_csv:
        return EcOH
    else:
        return table


def phylog(table):

    phylogroup = table.filter(regex=r'(^name|^phylogroup_)', axis=1)

    phylogroup = phylogroup.rename(columns=lambda x: re.sub(r'^phylogroup_','',x))

    phylogroup.set_index('name')

    B2_or_D = phylogroup.loc[phylogroup['chuA'] == 1]

    A_or_B1 = phylogroup.loc[phylogroup['chuA'] == 0]

    B2 = B2_or_D.loc[phylogroup['yjaA'] == 1]

    D = B2_or_D.loc[phylogroup['yjaA'] == 0]

    B1 = A_or_B1.loc[phylogroup['tspE4.C2'] == 1]

    A = A_or_B1.loc[phylogroup['tspE4.C2'] == 0]

    listdfs = []

    listdfs = listdfs.append([A,B1,B2,D])


    A['phylogroup'] = 'A'
    B1['phylogroup'] = 'B1'
    B2['phylogroup'] = 'B2'
    D['phylogroup'] = 'D'

    combined = pd.DataFrame()

    combined = pd.DataFrame(combined.append([A,B1,B2,D]))

    summary = combined.loc[:,['name','phylogroup']]

    return summary



if __name__ == '__main__':
    import argparse

parser = argparse.ArgumentParser(description='Process and join all genotypic and phylogenetic ARIBA output')
parser.add_argument('input_dir', help='Input directory containing csv files')
parser.add_argument('output_file', help='Output filename')
parser.add_argument('--clean', action='store_true', help='Use to clean anything after _R1 from filenames')
parser.add_argument('-chars','--characters', action='store_true', help='Pattern to trim after')
args = parser.parse_args()

if args.clean:
    print('\tRunning geno on CSV files in {} \n\t trimming characters in column \'name\' after \'{}\'...'.format(args.input_dir, args.characters))
    dflist1, names = geno(args.input_dir, args.clean, args.characters)
else:
    print('\tRunning geno on CSV files in {}... \n\tArgument \'--trim\' not used.'.format(args.input_dir))
    dflist1, names = geno(args.input_dir, args.clean, args.characters)

catnames = pd.concat(names)

catnames.reset_index(inplace=True)

names = catnames.drop_duplicates(keep='first')

print('\tCleaning headers of processed CSV files, adding identifiers prefixes')
dflist2 = simple_clean(dflist1)

print('\tGenerating MLST table')
mlst_table = mlst(args.input_dir)


mlst_table = mlst_table.set_index('name')

mlst_simple = mlst_table.iloc[:,0].to_frame()

full = mlst_table.join(dflist1)
simple = mlst_simple.join(dflist2)

EcOH = sero(simple, simple_csv=True)

simple = simple.reset_index()

phylogroup = phylog(simple)

simple = pd.merge(simple,EcOH,how='outer',on='name')

simple = pd.merge(simple,phylogroup,how='outer',on='name')

simple = simple.loc[:, ~simple.columns.str.startswith('sero_')]
simple = simple.loc[:, ~simple.columns.str.startswith('phylogroup_')]

simple = simple.set_index(['name','ST','phylogroup','O_type','H_type','OH_type'])

simple = simple.reindex_axis(sorted(simple.columns), axis=1)

simple = simple.reset_index()

simple = simple.loc[:, (simple != 0).any(axis=0)]

#Write to CSV
print('\tWriting simplified ARIBA table to {}simple.csv')
simple.to_csv(args.output_file+'simple.csv')
print('\tWriting full ARIBA table to {}simple.csv')
full.to_csv(args.output_file+'full.csv')
