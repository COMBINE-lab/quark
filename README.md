#quark

semi-reference-based short read compression


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
