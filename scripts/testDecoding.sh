while [[ $# > 1 ]]
do
    key="$1"

    case $key in
        --dec1)
            DECODED1="$2"
            shift # past argument
            ;;
        --dec2)
            DECODED2="$2"
            shift # past argument
            ;;
        --unmap1)
            UNMAP1="$2"
            shift # past argument
            ;;
        --unmap2)
            UNMAP2="$2"
            shift
            ;;
        --orig1)
            ORIG1="$2"
            shift
            ;;
        --orig2)
            ORIG2="$2"
            shift
            ;;
        *)
            # unknown option
            ;;
    esac
    shift # past argument or value
done

echo "comparing read sets"
cmp <(sort <(paste -d'' <(cat <(awk '{ if(NR % 4 == 2) { print $0; }}' ${UNMAP1}) ${DECODED1}) <(cat <(awk '{ if(NR % 4 == 2) { print $0; }}' ${UNMAP2}) ${DECODED2}))) <(sort <(paste -d'' <(awk '{ if(NR % 4 == 2) { print $0; }}' ${ORIG1})  <(awk '{ if(NR % 4 == 2) { print $0; }}' ${ORIG2})))
out=$?
if [ $out -eq 0 ]; then
    echo "Decoding was successful."
fi
