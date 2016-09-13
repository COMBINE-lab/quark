usage() { echo "Usage $0 [-1 left end] [-2 right end] [-s path to sailfish binary] [-i path to sailfish index] [-p #threads] [-o output Dir]" 1>&2; exit 1;}

while getopts ":1:2:s:i:p:o:" e; do
    case "${e}" in
        1)
            left=${OPTARG}
            ;;
        2)
            right=${OPTARG}
            ;;
        s)
            sailfish=${OPTARG}
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
        *)
            usage
            ;;
    esac
done

#echo "$#"
echo "{left read $left}"
echo "{right read $left}"
echo "sailfish binary $sailfish"
echo "sailfish index $ind"
echo "threads #$th"
echo "output dir $out"

if [ "$#" -ne 12  ]; then
    exit 1
fi

shift $((OPTIND-1))

#echo "$sailfish quant -i $ind -l IU -1 $left -2 $right -p $th -o $out"
$sailfish quant -i $ind -l IU -1 $left -2 $right -p $th -o $out
./decoder $out/aux/islands.quark $out/aux/reads.quark $out
