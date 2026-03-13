
## About
This is the repo for the Personal Genomics (CSE 284) project of Michael Iter, Dan Musachio, and Jonah Silverman.

We implemented a method to detect Identical By Descent (IBD) segments using unphased genotypes between two individuals.

To accomplish this we are created a Python based implementation of the tool TRUFFLE[^1]. We are calling our tool pyTRUFFLE.
TRUFFLE leverages long runs of Identical by State (IBS) SNPs to classify regions as IBD. Our tool is open source and easy to run.

[^1]: Dimitromanolakis A, Paterson AD, Sun L. Fast and Accurate Shared Segment Detection and Relatedness Estimation in Un-phased Genetic Data via TRUFFLE. Am J Hum Genet. 2019 Jul 3;105(1):78-88. doi: 10.1016/j.ajhg.2019.05.007. Epub 2019 Jun 6. PMID: 31178127; PMCID: PMC6612710.

# Data

To access the full datasets we ran, go to 'data_access_link.txt'. There will be a link to google drive containg full VCFs.
Further, for each access we included Chromesome 20 of a sibling pair (HG00581, HG00635) in this repo directly `data/HG00581_HG00635.vcf.gz`.

We also included full resutls of pyTRUFFLE and TRUFFLE on the full genomes of this sibling pair in `HG00581_vs_HG00635` and `HG00581_vs_HG00635_TRUFFLE` respectively.

# Usage

## Installation
Our tools runs in native python. All necessary packages can be installed using conda and the provided requirements.yaml
```
git clone https://github.com/Jonahs11/IBD_B_IBS_B_IDS.git

conda env create -f requirements.yaml
conda activate IBD_env
```
## Running
Data paths and parameter choices are organized using a config.json file. An example of one of these files can be found in the config_files directory. Many choices have defaults set, but the vcf path of interest and samp1/2 ids **must** be set in this file prior to running.


### Config Parameters

| Parameter | Description | Default |
|-----------|-------------|---------|
| `vcf` | Path to the VCF file containing samples to compare | — |
| `sample1` | ID of first sample | — |
| `sample2` | ID of second sample | — |
| `A` | Error cost | 1 |
| `B` | Error budget | 3 |
| `min_af` | Minimum minor allele fraction | 0.1 |
| `min_sig_length_ibd1` | Minimum contiguous length not breaking error budget to classify as IBD1. Based on analysis from the TRUFFLE paper. | 5 MB |
| `min_sig_length_ibd2` | Minimum contiguous length not breaking error budget to classify as IBD2. Based on analysis from the TRUFFLE paper. | 2 MB |
| `out_dir` | Directory to output results | — |


Once these fields are set, pass in the config file path using the --config flag.

```
python run_pyTRUFFLE.py --config <Path to config>
```

## Outputs

The `out_dir` will contain an `ibd_segments.tsv` file with IBD classifications for the genomic regions in the VCF file, structured as follows:

| Column | Description |
|--------|-------------|
| `chrom` | Chromosome |
| `start` | Start position |
| `end` | End position |
| `ibd_status` | IBD classification |



# TRUFFLE Usage

The TRUFFLE binary can be installed from the docs found at https://adimitromanolakis.github.io/truffle-website/.

Once installed here is the command we used to run truffle:
```
./truffle --vcf $file --segments --nofiltering #OPTIONALLY: --cpu 4
```
Note that we removed filtering because the default filtering was giving us weird results (e.g. IBD 2 = 0.04 for known full siblings)

### Outputs
- `truffle.ibd` holds the IBD 0,1,2 predicted percentages per pair in your VCF
- `truffle.segments` holds the predicted IBD segments in every chromosome you have data for


# Vignette

For a demo we will run and compare pyTRUFFLE (our implementation) against TRUFFLE on chromosome 20, for a sibling pair from the 1000Genomes dataset (HG00581 and HG00635). The vcf for these two samples are found in `./data/HG00581_HG00635.vcf.gz` of the repo.


First clone and create the conda environment if not done so already.
```
git clone https://github.com/Jonahs11/IBD_B_IBS_B_IDS.git
conda env create -f requirements.yaml
conda activate IBD_env

```

A config file is already created in `config_files/run1.json`

```
python run_pyTRUFFLE.py --config ./config_files/run1.json
```

The file `ibd_segments.tsv` should not exist in the `HG00581_vs_HG00635_chr20` directory.

The notebook `plot_ibd_results.ipynb` details how to plot this result as well as the provided full results.







 
