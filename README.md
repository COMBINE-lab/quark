#Quark 

semi-reference-based short read compression

##Assumption
The read files are in gzipped format i.e. they should be of the form 1.fastq.gz and 2.fastq.gz

The software is tested on paired end data, single end support would be added soon!!


##Compile
```{r, engine='bash', encode and decode}
git clone www.github.com/COMBINE-lab/quark.git
cd quark
mkdir build
cd build
cmake ..
make
cd ..
```

To see the options

```{r, engine='bash', encode and decode}
./mainscript.sh -h

```

To build the index


```{r, engine='bash', encode and decode}
./mainscript.sh index -t <transcript fasta> -o <out dir> -k <k mer length>

```

To Encode

```{r, engine='bash', encode and decode}
./mainscript.sh -1 <left_end> -2 <right_end> -i <index> -p <threads> -o <out dir>

```

To Decode

```{r, engine='bash', encode and decode}
./mainscript.sh -d decode -l [P/S] -i <input dir> -p <threads> -o <out dir>

```

To check the encoded and decoded sequences 


```{r, engine='bash', encode and decode}
./check_pair.sh <original left end> <original right end> <quark left end> <quark right end>

```
