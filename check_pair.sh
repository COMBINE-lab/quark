ORI1="$1"
ORI2="$2"

Q1="$3"
Q2="$4"

paste <(gunzip -c $ORI1) <(gunzip -c $ORI2) \
    | gawk 'BEGIN {n=0} {if(n++ % 4 == 1) print $0}' \
    | sort > original.seq

paste  $Q1 $Q2 \
    | gawk 'BEGIN {n=0} {if(n++ % 4 == 1) print $0}' \
    | sort > quark.seq

cmp original.seq quark.seq
if [ ${?} -eq 0 ] ; then
    echo "Looks ok; deleting files"
    rm original.seq
    rm quark.seq
fi
