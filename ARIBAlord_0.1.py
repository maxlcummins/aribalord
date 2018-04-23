import pandas as pd
import argparse

parser = argparse.ArgumentParser(description='Summarise gene hits from ARIBA output')
parser.add_argument('-f','--file',type=str, help='Path to file')
parser.add_argument('-o','--output',type=str, help='Output name')
args = parser.parse_args()

#read in csv file
df = pd.read_csv(args.file, index_col=0)

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
#   data_pivot = data_pivot.drop(data_pivot.columns[0])
    data_pivot = data_pivot.iloc[:,1:]
    data_pivot = data_pivot.fillna(0)
    frame.append(data_pivot)

#make list of original sample names using index of original df    
base = df.iloc[:,0:0]    

#join new dataframes from frames to base dataframe
data_all = base.join(frame)    

writefile = args.output+'.csv' 
   
#to write the file to disk unhash the line below and choose a filename
data_all.to_csv(writefile)