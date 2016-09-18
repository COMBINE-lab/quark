usage() { echo "Usage $0 [-1 left end] [-2 right end] [-i path to sailfish index] [-p #threads] [-o output Dir]" 1>&2; exit 1;}

echo "$#"
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
else
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


#echo "$#"
sailfish=~/Projects/quark/quickQuark/build/src/sailfish
leon=~/Projects/quark/quickQuark/leon/leon
scalce=~/Projects/quark/scalce/scalce
mince=/home/rob/mince/bin/mince


echo "mince binary $mince"
echo "scalce binary $scalce"
echo "leon binary $leon"
echo "sailfish binary $sailfish"
echo "left: $left"
echo "right: $right"
echo "sailfish index $ind"
echo "threads #$th"
echo "output dir $out"

pair=false
if [ "$#" -eq 10  ]; then
    pair=true
fi

shift $((OPTIND-1))
echo "$#"
#echo "$sailfish quant -i $ind -l IU -1 $left -2 $right -p $th -o $out"

if [ "$pair" = true ];then
    #sailfish module
    name=`echo $left | awk -F"/" '{print $NF}' | cut -d \. -f 1 | cut -d \_ -f 1`
    echo "Running Quark module"
    $sailfish quant -i $ind -l IU -1 <(gunzip -c $left) -2 <(gunzip -c $right) -p $th -o $out/quark
    plzip -k -f -n 20 $out/quark/aux/reads.quark
    plzip -k -f -n 20 $out/quark/aux/islands.quark
    $scalce $out/quark/aux/unmapped_1.fastq -r -o res -T 20 -p 100 -n library
    n1=$(stat -c%s res_1.scalcen)
    n2=$(stat -c%s res_2.scalcen)
    q1=$(stat -c%s res_1.scalceq)
    q2=$(stat -c%s res_2.scalceq)
    r1=$(stat -c%s res_1.scalcer)
    r2=$(stat -c%s res_2.scalcer)
    rq=$(stat -c%s $out/quark/aux/reads.quark.lz)
    iq=$(stat -c%s $out/quark/aux/islands.quark.lz)
    quarkSize=$((n1 + n2 + q1 + q2 + r1 + r2 + rq +iq))
    rm res_*
    #mince module
    echo "Running mince module"
    echo "$mince -e -l IU -1 <(gunzip -c $left) -2 <(gunzip -c $right) -p $th -o mince_"
    $mince -e -l IU -1 <(gunzip -c $left) -2 <(gunzip -c $right) -p $th -o mince_
    minceSize=$(stat -c%s *.lz | awk '{ sum += $1 }; END {print sum}')
    rm mince_*
    #scalce
    echo "Running scalce Module"
    echo "$scalce $left -r -o res -T 20 -p 100 -n library"
    $scalce $left -r -o res -T 20 -p 100 -n library
    nn=$(stat -c%s *.scalcen | awk '{ sum += $1  }; END {print sum}')
    qq=$(stat -c%s *.scalceq | awk '{ sum += $1  }; END {print sum}')
    rr=$(stat -c%s *.scalcer | awk '{ sum += $1  }; END {print sum}')
    scalceSize=$((nn + qq + rr))
    echo "$scalceSize"
    rm res_*
    #leon module
    echo "Running leon Module"
    echo "$leon -file $left -c -nb-cores 5 -seq-only"
    echo "$leon -file $right -c -nb-cores 5 -seq-only"
    $leon -file $left -c -nb-cores 5 -seq-only
    $leon -file $right -c -nb-cores 5 -seq-only
    leonSize=$(stat -c%s $out/*.leon | awk '{ sum += $1  }; END {print sum}')
    echo "$leonSize"
    rm *.leon
    echo "Q: $quarkSize M: $minceSize S: $scalceSize L: $leonSize"
    echo "$name Q: $quarkSize M: $minceSize S: $scalceSize L: $leonSize">>result.txt
else
    name=`echo $left | awk -F"/" '{print $NF}' | cut -d \. -f 1`
    echo "Running Quark module"
    $sailfish quant -i $ind -l U -r <(gunzip -c $read) -p $th -o $out/quark
    plzip -k -f -n 20 $out/quark/aux/reads.quark
    plzip -k -f -n 20 $out/quark/aux/islands.quark
    $scalce $out/quark/aux/unmapped.fastq -o res -T $th -p 100 -n library
    n1=$(stat -c%s *.scalcen)
    q1=$(stat -c%s *.scalceq)
    r1=$(stat -c%s *.scalcer)
    rq=$(stat -c%s $out/quark/aux/reads.quark.lz)
    iq=$(stat -c%s $out/quark/aux/islands.quark.lz)
    quarkSize=$((n1 + q1 + r1 + rq +iq))
    rm res_*
    #mince module
    echo "Running mince module"
    echo "$mince -e -l U -r <(gunzip -c $read) -p $th -o mince_"
    $mince -e -l U -r <(gunzip -c $read) -p $th -o mince_
    minceSize=$(stat -c%s *.lz | awk '{ sum += $1 }; END {print sum}')
    rm mince_*
    #scalce
    echo "Running scalce Module"
    echo "$scalce $read -o res -T $th -p 100 -n library"
    $scalce $read -o res -T $th -p 100 -n library
    nn=$(stat -c%s *.scalcen | awk '{ sum += $1  }; END {print sum}')
    qq=$(stat -c%s *.scalceq | awk '{ sum += $1  }; END {print sum}')
    rr=$(stat -c%s *.scalcer | awk '{ sum += $1  }; END {print sum}')
    scalceSize=$((nn + qq + rr))
    echo "$scalceSize"
    rm res_*
    #leon module
    echo "Running leon Module"
    echo "$leon -file $read -c -nb-cores 5 -seq-only"
    $leon -file $read -c -nb-cores 5 -seq-only
    leonSize=$(stat -c%s $out/*.leon | awk '{ sum += $1  }; END {print sum}')
    echo "$leonSize"
    rm *.leon
    echo "Q: $quarkSize M: $minceSize S: $scalceSize L: $leonSize"
    echo "$name Q: $quarkSize M: $minceSize S: $scalceSize L: $leonSize">>result.txt

fi
