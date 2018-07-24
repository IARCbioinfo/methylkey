<html>
<head>
</head>
<body>
</br></br>
<HR>
<h1> Differentially methylated Probes : {{ nbdmps }}</h1>
process : {{ process_id }}</br>
Model : {{ formula1 }}</br>
Method : {{ method }}</br>
Contrast matrix : {{ cmtx }}</br>
fdr cutoff : {{ fdr }}</br>
</br>
<a href="{{ out }}/{{ process_id }}/toptable.txt">TopTable ({{ nbdmps }})</a><br/>
<a href="{{ out }}/{{ process_id }}/annotated.txt">Complete annotation</a><br/>
Continue to R : <a href="{{ out }}/meth.rdata">RDATA</a><br/>
</br>

<table>
<tr><td><h3>QQplot</h3> Lambda = {{ lambda }}</td>
<td><h3>Top DMRs</h3></td>
<td><h3>Circus plot</h3></td>
</tr>
<tr>
<td><a href="{{ out }}/{{ process_id }}/qqplot.jpg"  ><img src="{{ out }}/{{ process_id }}/qqplot.jpg" height=420 width=420 /></a></td>
<td><a href="{{ out }}/{{ process_id }}/heatmap.jpg" ><img src="{{ out }}/{{ process_id }}/heatmap.jpg" height=420 width=420 /></a></td>
<td><a href="{{ out }}/{{ process_id }}/circleplot.jpg" ><img src="{{ out }}/{{ process_id }}/circleplot.jpg" height=420 width=420 /></a></td>
</tr>

</table>
</br>
Tip 1 : Dots on the red line of the qqplot are not significant while dots above it represent significantly differentilly methylated probes.</br>
Tip 2 : Lambda value reflect the qqplot's shape. It should be close to 1.</br>
Tip 3 : Lambda superior to 3, the qqplot is inflated you get a lot of significant dmps, try to improve your model by adding covariates.</br>
Tip 4 : Lambda inferior to 1, something wrong in the data or the model.</br></br>

</br>

<table>
<tr><td><h3>Gene structure annotation</h3></td>
<td><h3>CpG Island annotation</h3></td>
<td><h3>Histone Mark</h3></td>
</tr>
<tr>
<td><a href="{{ out }}/{{ process_id }}/cpgi_barplot.jpg"  ><img src="{{ out }}/{{ process_id }}/cpgi_barplot.jpg" height=420 width=420 /></a></td>
<td><a href="{{ out }}/{{ process_id }}/genes_barplot.jpg" ><img src="{{ out }}/{{ process_id }}/genes_barplot.jpg" height=420 width=420 /></a></td>
<td><a href="{{ out }}/{{ process_id }}/chromHMM_barplot.jpg" ><img src="{{ out }}/{{ process_id }}/chromHMM_barplot.jpg" height=420 width=420 /></a></td>
</tr>
</table>
</br>


</br></br>

</body>
</html>
