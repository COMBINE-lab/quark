usage() { echo "Usage $0 [-1 left end] [-2 right end]/[-r single end read] [-i path to sailfish index] [-p #threads] [-o output Dir]" 1>&2; 
          echo "Usage $0 -d [to decode] -l [P/S] [-i input directory] [-p #threads] [-o output Dir]" 1>&2; 
          echo "Usage $0 index -t <transcript file> -o <out dir> -k <kmer length>"
          echo "-h for usage"
          exit 1; 
        }


sailfish=$PWD/build/src/sailfish
decoder=$PWD/build/src/decoder
mince=$PWD/Mince-Binaries-0.6.1/mince_linux

makeindex() {
    echo "$sailfish index -t $1 -o $2 -k $3";
    $sailfish index -t $1 -o $2 -k $3;
    exit 1;
}

echo "$#"
decode=false

[[ $1 == "-h" || "$#" == 0 ]] && usage 

[[ $1 == "index" ]] && makeindex $3 $5 $7

if [[ $1 == "-d" ]]; then
    decode=true
fi

if [ "$decode" = false ];then
    if [ "$#" -eq 10 ];then
        while getopts ":1:2:i:p:o:" e; do
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
            esac
        done
    elif [ "$#" ];then
        while getopts ":r:i:p:o:" e; do
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
            esac
        done
    fi
else
    while getopts ":d:l:i:p:o:" e; do
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
        esac
    done
fi



pair=false
if [ "$#" -eq 10 ]; then
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

if [ "$decode" = false ];then
    if [ "$pair" = true ];then
        #encoding
        name=`echo $left | awk -F"/" '{print $NF}' | cut -d \. -f 1 | cut -d \_ -f 1`
        $sailfish quant -i $ind -l IU -1 <(gunzip -c $left) -2 <(gunzip -c $right) -p $th -o $out
        cd $out/aux
        echo $PWD
        $mince -e -l IU -1 unmapped_1.fastq -2 unmapped_2.fastq -p $th -o m_
        cd -
        cp $out/aux/*.lz $out/
        rm -r $out/aux
        rm -r $out/logs
        rm $out/cmd_info.json
    else
        $sailfish quant -i $ind -l U -r <(gunzip -c $read) -p $th -o $out
        $mince -e -l U -r $out/aux/unmapped.fastq -o $out/aux/m_
    fi
else
    if [ "$decodepair" = true ];then
        plzip -k -d $inputdir/islands.txt.lz
        cd $inputdir
        $mince -d -i m_ -o um_
        cd -
        $decoder $inputdir $out P
        rm $inputdir/um_*
    else
        $mince -d -i $out/aux/m_ -o $out/aux/um_
        $decoder $inputdir $out S
    fi
fi

