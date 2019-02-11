ind=$1
out=$2

cat /dev/null > $out.preatoms
cat /dev/null > $out.postatoms

a=`awk 'NR==1{print NF}' $ind`
b=`awk 'NR==2{print NF}' $ind`
c=`awk 'NR==3{print NF}' $ind`

i=$(($a+$b+$c))
for j in `seq $(($a+1)) $i`; do
    echo $j >> $out.preatoms
done

i=$(($a+$b))
for j in `seq $i`; do
    echo $j >> $out.postatoms
done
