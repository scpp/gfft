#!/bin/sh

## Example: ./gfftimport.sh res.txt

mysql_run="mysql -u myrnyy -s --column-names=FALSE -e"

config_fields="compiler_id software_id type_id value_type_id decimation_id parall_id array_type_id dim"
result_fields="hardware_id system_id output_type_id config_id"

## copy input file into import.txt
cp $1 import.txt

## import data from import.txt into temporary table gfft.import
mysqlimport -u myrnyy -d gfft `pwd`/import.txt

## construct field lists to build queries with gfft.config 
sel=""
fields=""
on_cond=""
i=0
for name in $config_fields
do
  i=`expr $i + 1`
  if [ $i -ge 2 ]; then
    sel="$sel,"
    fields="$fields,"
    on_cond="$on_cond and"
  fi
  sel="$sel tg.$name"
  fields="$fields $name"
  on_cond="$on_cond config.$name=tg.$name"
done

## select unique records from gfft.import to insert into gfft.config
import_select4config="(SELECT "\*" FROM gfft.import GROUP BY $fields) AS tg"
import_select4result="(SELECT "\*" FROM gfft.import GROUP BY $fields ,output_type_id ) AS tg"

## select the records that are not already in gfft.config
config_insert="SELECT $sel FROM gfft.config RIGHT JOIN $import_select4config ON ($on_cond) WHERE ISNULL(config_id)"
#$mysql_run "${result_query}"

## records will be inserted into gfft.config, if "insert" is given as the second parameter,
## otherwise the number of records will be displayed
if test "$2" = "insert"
then
$mysql_run "INSERT INTO gfft.config ($fields) ${config_insert}"
#echo "INSERT INTO gfft.config ($sel) ${config_insert}"
else
echo -n "Config: "
$mysql_run "SELECT COUNT(*) FROM (${config_insert}) AS t"
#$mysql_run "${config_insert}"
fi

config_query4result="SELECT gfft.config."\*",hardware_id,system_id,output_type_id FROM gfft.config RIGHT JOIN $import_select4result ON ($on_cond)"
#$mysql_run "${config_query}"
#echo "$config_query"

config_query4output="SELECT gfft.config."\*",hardware_id,system_id,output_type_id,input_value_id,value FROM gfft.config RIGHT JOIN  gfft.import AS tg ON ($on_cond)"
#$mysql_run "${config_query}"
#echo

##########################################################

## construct field lists to build queries with gfft.result
sel=""
fields=""
on_cond=""
i=0
for name in $result_fields
do
  i=`expr $i + 1`
  if [ $i -ge 2 ]; then
    sel="$sel,"
    fields="$fields,"
    on_cond="$on_cond and"
  fi
  sel="$sel tr.$name"
  fields="$fields $name"
  on_cond="$on_cond result.$name=tr.$name"
done

result_query="SELECT $sel FROM gfft.result RIGHT JOIN ($config_query4result) AS tr ON ($on_cond) WHERE ISNULL(result_id)"
#echo "$result_query"
#$mysql_run "${result_query}"

## records will be inserted into gfft.result, if "insert" is given as the second parameter,
## otherwise the number of records will be displayed
if test "$2" = "insert"
then
$mysql_run "INSERT INTO gfft.result ($fields) ${result_query} AND (tr.config_id>0)"
#echo "INSERT INTO gfft.result ($fields) ${result_query} WHERE config_id>0"
else
echo -n "Result: "
$mysql_run "SELECT COUNT(*) FROM (${result_query}) AS t"
#$mysql_run "${result_query}"
fi

## query of records to be inserted into gfft.output_value
output_query="SELECT input_value_id,result_id,value FROM gfft.result RIGHT JOIN ($config_query4output) AS tr ON ($on_cond)"

## check, if the records are already in gfft.output_value
output_query_check="SELECT COUNT(ttt.input_value_id) FROM gfft.output_value RIGHT JOIN ($output_query) AS ttt ON (output_value.input_value_id=ttt.input_value_id and output_value.result_id=ttt.result_id and output_value.value=ttt.value) WHERE ISNULL(output_value.result_id)"
#echo "$output_query"


if test `$mysql_run "${output_query_check}"` -ne 0
then
if test "$2" = "insert"
then
$mysql_run "INSERT INTO gfft.output_value ${output_query}"
#echo "INSERT INTO gfft.output_value ${output_query}"
else
echo -n "Output: "
$mysql_run "SELECT COUNT(*) FROM (${output_query}) AS t"
fi
else
echo "The file seems to be already added to the database!"
fi
