process_id : {{ process_id }}<br/>
Samples : {{ samples }}<br/>
Platform : {{ platform }}<br/>
Genome : {{ genome }}</br>
Pipeline : {{ pipeline }}<br/>
Normalisation : {{ normalize }}<br/>
Missing values : {{ missing }}<br/>
Excluded because NA > 20% (<a href="removed.txt">{{ filteredNAprobes }}</a>)<br/>
Cell mixture composition : {{ cell }} <br/>
<a href="report.html_files/pdata.txt">pdata</a><br/>
<a href="report.html_files/betas.txt">betas</a><br/>
<a href="report.html_files/deltabetas.txt">DeltaBetas</a><br/>
Continue to R : <a href="meth.rdata">RDATA</a><br/>

</br></br>
<HR>
<h1>Betas distribution by group</h1>

<select id="toogle" onchange="setPicture(this);">
{{ options }}
</select>

<table class="vignette">
<tr align="center">
<td><a href="{{ group }}/densityPlot2.jpg"><img class="gpt" src="{{ group }}/densityPlot2.jpg" width="500" height="500"/></a></td>
<td><a href="{{ group }}/densityBeanPlot2.jpg"><img class="gpt" src="{{ group }}/densityBeanPlot2.jpg" width="500" height="500"/></a></td>
</tr>
<tr>
<td><a href="{{ group }}/mdsPlot2.jpg"><img class="gpt"  src="{{ group }}/mdsPlot2.jpg" width="500" height="500"/></a></td>
<td></td>
</tr>
</table>

