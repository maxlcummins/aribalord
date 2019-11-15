# ARIBAlord
ARIBA allows for a fast, efficient and scalable approach to the genotyping of short-read sequence data. The inbuilt ARIBA summary function provides in-depth insight into the quantity and quality of genes detected within short-read sequence data, however the cluster ID's generated through the use of CD-hit can complicate the en-masse generation of easily interpretable reports detailing the specific reference sequences identified within a given sample's sequence data.

Enter ARIBAlord; a script used for the processing of ARIBA data. It streamlines the processing of genotypic and phylogenetic data produced by ARIBA and allows serotype and phylogroup classification of E. coli isolates.

## NOTE: THIS SOFTWARE IS OPTIMISED FOR E. coli AND MAY NOT WORK FOR OTHER GENERA!

## Getting Started

No information will be given on the installation or usage of ARIBA.
For usage of ARIBA see the [ARIBA](https://github.com/sanger-pathogens/ariba) github.

### Prerequisites

  * [Python3][python] version >= 3.3.?

There are also several Python packages upon which ARIBAlord was built:
  * pandas
  * numpy
  * argparse
  * regex
  
 If you're using conda, try this:
 
 ```
 conda create -n aribalord pandas regex argparse
 ```
  
## Data Preprocessing

Following read-mapping by ARIBA to generate phylogenetic and genotypic data there is one of two pre-processing steps required for the generation of summary files, depending on the type of data to be processed.

### General pre-processing
Generation of a summary file using the ARIBA function 'summary' is required for the processing of all ARIBA output by ARIBAlord, with the only exception being a MLST summary file, which will soon be described.

To generate the summary files in a format that ARIBAlord can handle, run the following command:

```
ARIBA summary --cluster_cols assembled,ref_seq <prefix> *<ARIBA_output_directory/report.tsv>
```

### MLST pre-processing
Following the instructions of ARIBA to undertake multi-locus sequence typing, two steps are required to prepare MLST outputs for processing by ARIBAlord.

The first command appends each row of each MLST report with the name of the file so that following concatenation each row has the appropriate ID. The output is written to a new file (appended with '.new') to maintain the integrity of the original files.

The following command should do the trick:
```
for f in *<ARIBA_output_directory>/mlst_report.tsv; do paste $f <(yes $f | head -n $(cat $f | wc -l)) > $f.new; done
```
Next we concatenate these files using cat:
```
cat *<ARIBA_output_directory>/mlst_report.tsv.new > <output.tsv>
```

The concatenated output is now ready for processing with ARIBAlord.

Note: This concatenated output also concatenates the header from each file. For your viewing pleasure, you can delete the extra header lines as follows, however this step is not required if you plan to run the file through the ARIBAlord.

```
grep -v ^ST <output.tsv> > <noheaders_output.tsv>
```

## ARIBAlord Usage
ARIBAlord arguments can be viewed using the following command:

```
python ARIBAlord.py -h
```

Basic usage is as follows.
```
python ARIBAlord.py input_directory output_prefix
```
Note that no CSVs can be present in the input directory other than those to be processed by ARIBAlord. Ensure they are formatted correctly, with two columns for each gene hit, gene_assembled and gene_ref_seq, as described above. Also make sure that the MLST file to be processed, if present, ends in MLST.tsv, and that Phylogroup and Serotype files end in Phylogroup.csv and EcOH.csv, respectively. Note that ARIBAlord can process any combination of MLST, Phylogroup and Serotype data, and will work fine in absence of all files if you aren't interested in generating such data.

## Known Issues/Upcoming Changes
* Untested effect of having non-fliC H alleles
* fimH typing to be introduced
* Support for MLST of species other than E. coli


## Authors

* **Max Cummins** - *Initial work* - [maxlcummins](https://github.com/maxlcummins)
* **Cameron Reid** - *Initial work* - [CJReid](https://github.com/CJReid)

## License
                    GNU GENERAL PUBLIC LICENSE
                       Version 3, 29 June 2007

 Copyright (C) 2007 Free Software Foundation, Inc. <http://fsf.org/>
 Everyone is permitted to copy and distribute verbatim copies
 of this license document, but changing it is not allowed.

See LICENSE file for more details.

## Acknowledgments
* Michael Liu
