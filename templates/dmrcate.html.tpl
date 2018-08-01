<html>
<head>
</head>
<body>
</br></br>
<HR>
<h1> Differentially methylated Regions : {{ nbdmrs }}</h1>
process : {{ process_id }}</br>
Model : {{ formula1 }}</br>
Method : {{ method }}</br>
Contrast matrix : {{ cmtx }}</br>
fdr cutoff : {{ fdr }}</br>
pvalue cutoff : {{ pcutoff }}</br>
</br>
<a href="{{ process_id }}/toptable.txt">TopTable ({{ nbdmrs }})</a><br/>
<a href="{{ process_id }}/annotated.txt">Complete annotation</a><br/>
Continue to R : <a href="meth.rdata">RDATA</a><br/>
</br>

<table>
<tr><td><h3>Gene structure annotation</h3><br/></td>
<td><h3>CpG Island annotation</h3></td>
</tr>
<tr>
<td><a href="{{ process_id }}/cpgi_barplot.jpg"  ><img src="{{ process_id }}/cpgi_barplot.jpg" height=420 width=420 /></a></td>
<td><a href="{{ process_id }}/genes_barplot.jpg" ><img src="{{ process_id }}/genes_barplot.jpg" height=420 width=420 /></a></td>
</tr>
</table>
</br>

</body>
</html>
