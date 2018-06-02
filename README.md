# ARIBAlord
ARIBA allows for a fast, efficient and scalable approach to the genotyping of short-read sequence data. The inbuilt ARIBA summary function provides in-depth insight into the quantity and quality of genes detected within short-read sequence data, however the cluster ID's generated through the use of CD-hit convolute the en-masse generation of easily interpretable reports detailing the specific reference sequences identified for a given sample.

Enter ARIBAlord; a series of scripts used for the processing of ARIBA data. It streamlines the processing of genotypic and phylogenetic data (multi-locus sequence type, serotype and phylogroup) and streamlines the conglomeration of such data into phylogenetic, genotypic or combined tables.



## Getting Started

No information will be given on the installation or usage of ARIBA.
For usage of ARIBA see the [ARIBA](https://github.com/sanger-pathogens/ariba) github.

### Prerequisites

  * [Python3][python] version >= 3.3.?

There are also several Python packages upon which ARIBAlord was built:
  * pandas
  * numpy
  * argparse
  
## Data Preprocessing

Following read-mapping by ARIBA to generate phylogenetic and genotypic data there is one of two pre-processing steps required for the generation of summary files, depending on the type of data to be processed.

### General pre-processing
Generation of a summary file using the ARIBA function 'summary' is required for the processing of all ARIBA output by ARIBAlord, with the only exception being a MLST summary file, which will soon be described.

To generate the summary files in a format that ARIBAlord can handle, run the following command:

```
ARIBA summary --cluster_cols assembled,ref_seq <prefix> <*_ARIBA_output_prefix/report.tsv>
```

### MLST pre-processing
Following the instructions of ARIBA to undertake multi-locus sequence typing, two steps are required to prepare MLST outputs for processing by ARIBAlord.

The first command appends each row of each MLST report with the name of the file so that following concatenation each row has the appropriate ID. The output is written to a new file (appended with '.new') to maintain the integrity of the original files.

The following command should do the trick:
```
for f in *<ARIBA_output_prefix>/mlst_report.tsv; do paste $f <(yes $f | head -n $(cat $f | wc -l)) > $f.new; done
```
Next we concatenate these files using cat:
```
for f in *<ARIBA_output_prefix>/mlst_report.tsv.new; do cat $f > <output.tsv>
```

Then we need to strip unnecessary characters from the sample names to facilitate downstream processing. In this example, the read filename is appended by \_R1.fastq.gz. Feel free to edit the regular expression in the following script if your filename doesn't conform to this structure. Not that the following command replaces the string of interest in place, therefore you may want to ensure you get your string substitution correct by running the replacement command on a copy of your ARIBA summary files or else you might have to regenerate them.
```
perl -p -i -e 's/_R1\.fastq\.gz.*\/.*report\.tsv//g' <output.tsv>
```
The concatenated output is now ready for processing with ARIBAlord.

Note: This concatenated output also concatenates the header from each file. For your viewing please, you can delete the extra header lines as follows, however this step is not required if you plan to run the file through the ARIBAlord subcommands MLST, phylojoin or genophy.

```
grep -v ^ST <output.tsv> > <noheaders_output.tsv>
```

## ARIBAlord Usage
ARIBAlord subcommands can be viewed using the following command:

```
python ARIBAlord.py -h
```
The arguments for the various ARIBAlord subcommands can be viewed similarly, for example:
```
python genophy -h
```

## Authors

* **Max Cummins** - *Initial work* - [maxlcummins](https://github.com/maxlcummins)
* **Cameron Reid** - *Initial work* - [CJReid](https://github.com/CJReid)

See also the list of [contributors](https://github.com/your/project/contributors) who participated in this project.

## License

This project is licensed under an MIT license:

MIT License

Copyright (c) 2018 maxlcummins

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

## Acknowledgments
* Michael Liu
