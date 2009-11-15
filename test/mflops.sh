#!/bin/sh

queryid=
#queryhead="SELECT input_value.value AS p, POW(2,input_value.value) AS n"
queryhead="SELECT input_value.value AS p"
querybody=
queryconnect="input_value.input_value_id=t$1.input_value_id"

for i in $*
do 
  queryid=${queryid}$i
  s="SELECT input_value_id,value FROM gfft.output_value WHERE result_id=$i"
  if [ $i -ne $1 ]; then
    querybody="${querybody}, ";
    queryconnect="${queryconnect} AND t$1.input_value_id=t$i.input_value_id";
  fi
  queryhead="${queryhead},POW(2,input_value.value)*5*input_value.value/t$i.value/1000000 AS CPUtime$i"
#  queryhead="${queryhead},t$i.value AS CPUtime$i"
  querybody="${querybody}(${s}) AS t$i"
done

#fname="/home/myrnyy/db/gfft/compare${queryid}.txt"
#query="${queryhead} INTO OUTFILE '${fname}' FROM gfft.input_value, ${querybody} WHERE ${queryconnect};"
query="${queryhead} FROM gfft.input_value, ${querybody} WHERE ${queryconnect} ORDER BY gfft.input_value.input_value_id;"
#echo $query > tmpsql.txt
mysql -u myrnyy -s --column-names=FALSE -e "${query}" > tmpquery.txt

