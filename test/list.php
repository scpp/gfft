<html>
  <head>
    <title>GFFT Database</title>
   <meta http-equiv="Content-Type" content="text/html; charset=iso-8859-1">
  </head>
<body bgcolor="#F4F2DB">
<H2>Data analysis</H2>
<?php

function makecombobox($db_link,$table,$idname,$where) 
{
        $fid = "{$idname}_id";
	$query = "SELECT {$table}.{$fid},name FROM $table";
#	if (strlen($where)>0) {
 #          $query = "$query INNER JOIN (SELECT * FROM result_config_output GROUP BY $fid HAVING $where) AS t1 ON {$table}.{$fid}=t1.{$fid} ORDER BY $fid ;";
  #      }
#echo $query,"<br>";
	echo "<td>{$idname}<br>";
	$result = mysql_query($query, $db_link);
#echo $result,"<br>";
	echo "<select name=\"${idname}[]\" size=\"10\" multiple=\"multiple\">\n";
	echo "<option value=\"-1\">all</option>\n";
	while ($row = mysql_fetch_row($result)) {
		echo "<option value=\"$row[0]\"";
		if (!empty($_POST[$idname])) {
		  for ($i=0; $i<count($_POST[$idname]); $i++) {
		      $k = $_POST[$idname][$i];
			if ($row[0]==$k) echo " selected=\"selected\"";
		  }
		}
		echo ">$row[1]</option>\n";
	}
	echo '</select>';
	echo "</td>";
}

	$idnames = array('software','hardware','system','compiler','type','parall','decimation','array_type','value_type','output_type');
	$relnames = array('config','result','result','config','config','config','config','config','config','result');
	$tnames = array('gfft.software_view','common.hardware','common.system','gfft.compiler_view','gfft.fft_type','gfft.fft_parall','gfft.fft_decimation','gfft.fft_array_type','gfft.fft_value_type','gfft.output_type');
	$tsnames = array('software_view','hardware','system','compiler_view','fft_type','fft_parall','fft_decimation','fft_array_type','fft_value_type','output_type');
	$qwhere = '';
# build string $qwhere
if ($_POST) {
	for ($i=0; $i<count($idnames); $i++) {
		$value = $idnames[$i];
		if (!empty($_POST[$value])) {
		  $qwhere = "$qwhere (";
		  for ($j=0; $j<count($_POST[$value]); $j++) {
		      $k = $_POST[$value][$j];
		      $qwhere = "$qwhere {$value}_id=$k OR";
		  }
		  $qwhere = substr($qwhere,0,-3).") AND";
		}
	}
	if (strlen($qwhere)>0) {
		$qwhere = substr($qwhere,0,-4);        // remove last " AND"
        }
}

	$db_link = mysql_connect('localhost', 'myrnyy', '');
	//echo "$query";
	echo "<form method=\"POST\" action=\"list.php\">\n";
	#echo "<fieldset><legend>Config</legend>\n";
	echo "<table><tr>";
	for ($i=0; $i<count($idnames); $i++) 
		makecombobox($db_link,$tnames[$i],$idnames[$i],$qwhere);
	echo "</tr></table>\n";
	echo "<table width=100%>\n";
	echo "<tr><td valign=top>\n";
	echo "Legend<br>\n";
	echo "<select name=\"legend[]\" size=\"10\" multiple=\"multiple\">\n";
	// This variable saves first non-set (=='all') config entry
	$def_legend_id = -1;
	$k = 0;
	for ($i=0; $i<count($idnames); $i++) {
		if (($_POST && empty($_POST[$idnames[$i]])) || (!$_POST)) {
			echo "<option value=\"$i\"";
			if (!empty($_POST['legend']) && $k<count($_POST['legend'])) {
				if ($_POST['legend'][$k]==$i) {
					$k++;
					echo " selected=\"selected\"";
				}
			}
			echo ">$idnames[$i]</option>\n";
		}
		if ($_POST && empty($_POST[$idnames[$i]]) && $def_legend_id==-1) 
			$def_legend_id = $i;
	}
	if ($def_legend_id==-1) $def_legend_id = 0;
	echo '</select><br><br>';
	echo '<input type="submit" value="Query">';
	echo '</td>';

if ($_POST) {
	$qfrom = "FROM gfft.result INNER JOIN gfft.config ON result.config_id=config.config_id";
	$query = "SELECT result_id $qfrom";
#echo "where: ",$qwhere,"<br><br>";
	$legendfields = '';
	$legendtables = '';
	$legendon = '';
	if (!empty($_POST['legend'])) {
	for ($i=0; $i<count($_POST['legend']); $i++) {
		$k = $_POST['legend'][$i];
		$value = $idnames[$k];
		$legendfields = "$legendfields,',',$tsnames[$k].name";
		$legendtables = "$legendtables,$tnames[$k]";
		$legendon = "$legendon AND $tsnames[$k].{$value}_id=$relnames[$k].{$value}_id";
	}
	$legendfields = "concat(".substr($legendfields,5).")";
	} else {    // default legend query variables
		$value = $idnames[$def_legend_id];
		$legendfields = "$tsnames[$def_legend_id].name";
		$legendtables = ",$tnames[$def_legend_id]";
		$legendon = " AND $tsnames[$def_legend_id].{$value}_id=$relnames[$def_legend_id].{$value}_id";
	}
	$legendquery = "SELECT $legendfields FROM gfft.result INNER JOIN (gfft.config$legendtables) ON (result.config_id=config.config_id $legendon)";
	if (strlen($qwhere)>0) {
		$query = "$query WHERE $qwhere ORDER BY result_id";
		$legendquery = "$legendquery WHERE $qwhere ORDER BY result_id";
	}
#echo $query,"<br><br>";
	$result = mysql_query($query, $db_link);
	$nres = mysql_num_rows($result);

	if ($nres>0 && $nres<=10) { 
		$resstring = '';
		while ($row = mysql_fetch_row($result)) {
			$resstring = "$resstring$row[0] ";
		}

		$plot = "set terminal png notransparent font arial 12\nset output 'tmpplot.png'\nset grid\nset xrange [1:]\nset xlabel 'Power of two'\nset ylabel 'MFLOPS'\nplot";
		$result = mysql_query($legendquery, $db_link);
		$i = 1;
		while ($row = mysql_fetch_row($result)) {
			$i++;
			$plot = "$plot 'tmpquery.txt' using 1:$i title '$row[0]' with lines lw 2 lc $i";
			if ($i<=$nres) $plot = "$plot,";
		}
	
#echo "Legend: ",$legendquery,"<br><br>";
	
		$handle = fopen("tmpplot.plt","w");
		fwrite($handle,$plot);
		fclose($handle);
echo $resstring,"<br>";
		system("/home/vladimir/db/gfft/mflops.sh $resstring");
		system("gnuplot tmpplot.plt");
		echo '<td halign="center"><img src="tmpplot.png"></td>';
	}
	echo '</tr>';
	echo '</table>';
	echo '</form>';
	echo '<hr>';
	echo "Number of results: $nres<br>";
}
	mysql_close($db_link);

?>
  </body>
</html>
