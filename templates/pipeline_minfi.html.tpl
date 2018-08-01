process_id : {{ process_id }}<br/>
Samples : {{ samples }}<br/>
Platform : {{ platform }}<br/>
Pipeline : {{ pipeline }}<br/>
Normalisation : {{ normalize }}<br/>
Probes : {{ nbprobes2 }}/{{ nbprobes1 }} ({{ nbprobes1-nbprobes2 }} filtered) <br/>
Excluded from list = {{ filteredFromList }} ( {{ filters }} )</br>
Missing values : {{ missing }}<br/>
Excluded because NA > 20% (<a href="removed.txt">{{ filteredNAprobes }}</a>)<br/>
Cell mixture composition : {{ cell }} <br/>
<a href="pdata.txt">pdata</a><br/>
<a href="betas.txt">betas</a><br/>
<a href="deltabetas.txt">DeltaBetas</a><br/>
Continue to R : <a href="meth.rdata">RDATA</a><br/>

</br></br>
<HR>
<h1>Signal intensities distribution</h1>

<table class="vignette">
<tr>
<td><h3>Boxplot of Red (methylated) and Green (unmethylated) Color Channels :</h3></br>
Tip 1: The median should be aligned, a box very low compared to other is suspicious and should be remove from the analysis.</br>
Tip 2: The 8th row of each sentrix is often lower compared to others. This is corrected by the normalization
</td>
<td><h3>Ratio of median intensities</h3></br>
Tip 1: All samples above the line (black) are good quality, while all sample above (red) are bad quality.
Tip 2: Fresh Frozen tissus are often bad quality
</td>
</tr>
<tr>
<td><a href="boxplot_colour1.jpg"><img src="boxplot_colour1.jpg" width="500" height="500" /></a></td>
<td><a href="plotQC.jpg"><img src="plotQC.jpg" width="500" height="500" /></a></td>
</tr>
<tr>
<td><h2>Boxplot of Red (methylated) Color Channels :</h2></br>
</td>
<td><h2>Boxplot of Green (unmethylated) Color Channels :</h2></br>
</td>
</tr>
<tr>
<td><a href="methylated1.jpg"><img src="methylated1.jpg" width="500" height="500" /></a></td>
<td><a href="unmethylated1.jpg"><img src="unmethylated1.jpg" width="500" height="500" /></a></td>
</tr>
</table>

</br></br>
<HR>
<h1>Betas distribution by group</h1>

<select id="toogle" onchange="setPicture(this);">
{{ options }}
</select>

<table class="vignette">
<tr align="center">
<td><h3>Raw</h3></td><td><h3>Filtered & Normalized</h3></td>
</tr>
<tr>
<td><a href="{{ group }}/densityPlot1.jpg"><img class="gpt" src="{{ group }}/densityPlot1.jpg" width="500" height="500"/></a></td>
<td><a href="{{ group }}/densityPlot2.jpg"><img class="gpt"  src="{{ group }}/densityPlot2.jpg" width="500" height="500"/></a></td>
</tr>
<tr>
<td><a href="{{ group }}/densityBeanPlot1.jpg"><img class="gpt" src="{{ group }}/densityBeanPlot1.jpg" width="500" height="500"/></a></td>
<td><a href="{{ group }}/densityBeanPlot2.jpg"><img class="gpt"  src="{{ group }}/densityBeanPlot2.jpg" width="500" height="500"/></a></td>
</tr>
<tr>
<td><a href="{{ group }}/mdsPlot1.jpg"><img class="gpt"  src="{{ group }}/mdsPlot1.jpg" width="500" height="500"/></a></td>
<td><a href="{{ group }}/mdsPlot2.jpg"><img class="gpt"  src="{{ group }}/mdsPlot2.jpg" width="500" height="500"/></a></td>
</tr>
</table>

