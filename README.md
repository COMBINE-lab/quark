Quark
=====

semi-reference-based short read compression

<p align="center">
<img src="qimage.001.png">
</p>


## Assumption

The read files are in gzipped format i.e. they should be like .. 1.fastq.gz and 2.fastq.gz

The software is tested on paired end and single end data on bash compatible shell (redirection might not work with *fish* kind of ad on), ~~*single end support will be added to the "quark.sh" script soon*.~~

## Dependency

Quark depends on `plzip` for downstream compression. More information about Plzip and installation guide can be found [here](http://www.nongnu.org/lzip/plzip.html). 



## Compile

```{r, engine='bash', encode and decode}
$git clone www.github.com/COMBINE-lab/quark.git
$cd quark
$mkdir build
$cd build
$cmake ..
$make
$cd ..
```

##Running Quark

To see the options

```{r, engine='bash', encode and decode}
$./quark.sh -h

```

### To build the index with kmer size k

```{r, engine='bash', encode and decode}
snakemake -s quark.snake make_index --config out="<output dir>" fasta="<fasta file>" kmer=<#k>
```


### To Encode

#### Single End

```{r, engine='bash', encode and decode}
snakemake -s quark.snake encode --config out="<output dir>" index="<index dir>" r="<mate>" p=<#threads> lib="single" quality=0
```


#### Paired end

```{r, engine='bash', encode and decode}
snakemake -s quark.snake encode --config out="<output dir>" index="<index dir>" m1="<mate1>" m2="<mate2>" p=<#threads> lib="paired" quality=0
```

### To Decode

```{r, engine='bash', encode and decode}
snakemake -s quark.snake decode --config in="<in dir>" out="<out dir>" lib="paired/single" quality=0
``` 



## To check the encoded and decoded sequences are same !! (it is lossless) 


```{r, engine='bash', encode and decode}
$./check_pair.sh <original left end> <original right end> <quark left end> <quark right end>

```
## Link to the preprint

[Quark enables semi-reference-based compression of RNA-seq data](http://dx.doi.org/10.1101/085878) by Hirak Sarkar, Rob Patro

