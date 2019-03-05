#!/bin/bash

usage() { echo "Usage $0 [-1 left end] [-2 right end]/[-r single end read] [-i path to sailfish index] [-p #threads] [-o output Dir]" 1>&2;
          echo "Usage $0  -l [P/S] [-i input directory] [-p #threads] [-o output Dir]" 1>&2;
          echo "Usage $0 index -t <transcript file> -o <out dir> -k <kmer length>"
          echo "-h for usage"
          exit 1;
        }


SCRIPT=$(readlink -f "$0")
SCRIPTPATH=$(dirname "$SCRIPT")

quark=$SCRIPTPATH/build/src/quark
decoder=$SCRIPTPATH/build/src/decoder
mince=$SCRIPTPATH/Mince-Binaries-0.6.1/mince_linux
echo "$SCRIPTPATH"

makeindex() {
    echo "$quark index -t $1 -o $2 -k $3";
    $quark index -t $1 -o $2 -k $3;
    exit 1;
}

echo "$#"
decode=false
quality=false

[[ $1 == "-h" || "$#" == 0 ]] && usage

[[ $1 == "index" ]] && makeindex $3 $5 $7

if [[ $1 == "-d" ]]; then
    decode=true
fi

if [ "$decode" = false ];then
    if [ "$#" -eq 12 ];then
        while getopts ":1:2:i:p:o:q:" e; do
            case "${e}" in
                1)
                    left=${OPTARG}
                    ;;
                2)
                    right=${OPTARG}
                    ;;
                i)
                    ind=${OPTARG}
                    ;;
                p)
                    th=${OPTARG}
                    ;;
                o)
                    out=${OPTARG}
                    ;;
                q)
                    quality=${OPTARG}
                    ;;
            esac
        done
    elif [ "$#" -eq 10 ];then
        while getopts ":r:i:p:o:q:" e; do
            case "${e}" in
                r)
                    read=${OPTARG}
                    ;;
                i)
                    ind=${OPTARG}
                    ;;
                p)
                    th=${OPTARG}
                    ;;
                o)
                    out=${OPTARG}
                    ;;
                q)
                    quality=${OPTARG}
                    ;;
            esac
        done
    fi
else
    while getopts ":d:l:i:p:o:q:" e; do
        case "${e}" in
            d)
                dec=${OPTARG}
                ;;
            l)
                libtype=${OPTARG}
                ;;
            i)
                inputdir=${OPTARG}
                ;;
            p)
                th=${OPTARG}
                ;;
            o)
                out=${OPTARG}
                ;;
            q)
                quality=${OPTARG}
                ;;
        esac
    done
fi



pair=false
if [ "$#" -eq 12 ]; then
    pair=true
fi


shift $((OPTIND-1))
echo "$#"

decodepair=false
if [[ $libtype == "P" ]];then
    decodepair=true
fi

if [ ! -d $out ];then
    mkdir -p $out;
fi
echo "quality $quality"

if [[ $quality == "1" ]];then
    quality=true
fi

if [ "$decode" = false ];then
    if [ "$pair" = true ];then
        #encoding
        name=`echo $left | awk -F"/" '{print $NF}' | cut -d \. -f 1 | cut -d \_ -f 1`
        if [ "$quality" = true ];then
            echo "Encoding with quality score"
            echo "$quark quant -i $ind -l IU -1 <(gunzip -c $left) -2 <(gunzip -c $right) -p $th -o $out --quality"
            $quark quant -i $ind -l IU -1 <(gunzip -c $left) -2 <(gunzip -c $right) -p $th -o $out --quality
        else
            echo "Encoding without quality score"
            $quark quant -i $ind -l IU -1 <(gunzip -c $left) -2 <(gunzip -c $right) -p $th -o $out
        fi
        cd $out/aux
        echo $PWD
        $mince -e -l IU -1 unmapped_1.fastq -2 unmapped_2.fastq -p $th -o m_
        cd -
        cp $out/aux/*.lz $out/
        rm -r $out/aux
        rm -r $out/logs
        rm $out/cmd_info.json
    else
        if [ "$quality" = true ];then
            $quark quant -i $ind -l U -r <(gunzip -c $read) -p $th -o $out --quality
        else
            $quark quant -i $ind -l U -r <(gunzip -c $read) -p $th -o $out
        fi
        cd $out/aux
        $mince -e -l U -r unmapped.fastq -o m_
        cd -
        cp $out/aux/*.lz $out/
        rm -r $out/aux
        rm -r $out/logs
        rm $out/cmd_info.json
    fi
else
    echo "decompressing $inputdir/islands.txt.lz ... "
    plzip -k -d $inputdir/islands.txt.lz
    echo "entering $inputdir"
    cd $inputdir
    $mince -d -i m_ -o um_
    cd -
    if [ "$decodepair" = true ];then
        if [ "$quality" = true ];then
            echo "decoding with quality score"
            plzip -k -d $inputdir/quality_1.quark.lz
            plzip -k -d $inputdir/quality_2.quark.lz
            $decoder $inputdir $out P Q
        else
            echo "decoding without quality score"
            $decoder $inputdir $out P N
        fi
    else
        if [ "$quality" = true ];then
            plzip -k -d $inputdir/quality.quark.lz
            $decoder $inputdir $out S Q
        else
            $decoder $inputdir $out S N
        fi
    fi
    # rm $inputdir/um_*
    # rm $inputdir/islands.txt
fi
