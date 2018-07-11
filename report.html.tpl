<html>
<head>
<style>
.vignette table {
border-collapse: collapse;
border: 1px solid black;
}
.vignette td{border : 1px solid black; width:100px;}
</style>
</head>
<body>
process_id : {{ process_id }}
Samples : {{ samples }}<br/>
Platform : {{ platform }}<br/>
Pipeline : {{ pipeline }}<br/>
Normalisation : {{ normalize }}<br/>
Probes : {{ nbprobes3 }}/{{ nbprobes1 }} ({{ nbprobes1-nbprobes3 }} filtered) <br/>
Excluded from list = {{ nbprobes1-nbprobes2 }} ( SNP_EPIC.csv,Sex_EPIC.csv )</br>
Missing values : {{ missing }}<br/>
Excluded because NA > 20% (<a href="removed.txt">{{ nbprobes2-nbprobes3 }}</a>)<br/>
Cell mixture composition : {{ cell }} <br/>
<a href="report.html_files/pdata.txt">pdata</a><br/>
<a href="report.html_files/betas.txt">betas</a><br/>
<a href="report.html_files/deltabetas.txt">DeltaBetas</a><br/>

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
Only with minfi !</br>
Tip 1: All samples above the line (black) are good quality, while all sample above (red) are bad quality.
Tip 2: Fresh Frozen tissus are often bad quality
</td>
</tr>
<tr>
<td><a href="{{ out }}/boxplot_colour1.jpg"><img src="{{ out }}/boxplot_colour1.jpg" width="500" height="500" /></a></td>
<td><a href="{{ out }}/plotQC.jpg"><img src="{{ out }}/plotQC.jpg" width="500" height="500" /></a></td>
</tr>
<tr>
<td><h2>Boxplot of Red (methylated) Color Channels :</h2></br>
</td>
<td><h2>Boxplot of Green (unmethylated) Color Channels :</h2></br>
</td>
</tr>
<tr>
<td><a href="{{ out }}/methylated1.jpg"><img src="{{ out }}/methylated1.jpg" width="500" height="500" /></a></td>
<td><a href="{{ out }}/unmethylated1.jpg"><img src="{{ out }}/unmethylated1.jpg" width="500" height="500" /></a></td>
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
<td><a href="{{ out }}/{{ group }}/densityPlot1.jpg"><img class="gpt" src="{{ out }}/{{ group }}/densityPlot1.jpg" width="500" height="500"/></a></td>
<td><a href="{{ out }}/{{ group }}/densityPlot2.jpg"><img class="gpt"  src="{{ out }}/{{ group }}/densityPlot2.jpg" width="500" height="500"/></a></td>
</tr>
<tr>
<td><a href="{{ out }}/{{ group }}/densityBeanPlot1.jpg"><img class="gpt" src="{{ out }}/{{ group }}/densityBeanPlot1.jpg" width="500" height="500"/></a></td>
<td><a href="{{ out }}/{{ group }}/densityBeanPlot2.jpg"><img class="gpt"  src="{{ out }}/{{ group }}/densityBeanPlot2.jpg" width="500" height="500"/></a></td>
</tr>
<tr>
<td><a href="{{ out }}/{{ group }}/mdsPlot1.jpg"><img class="gpt"  src="{{ out }}/{{ group }}/mdsPlot1.jpg" width="500" height="500"/></a></td>
<td><a href="{{ out }}/{{ group }}/mdsPlot2.jpg"><img class="gpt"  src="{{ out }}/{{ group }}/mdsPlot2.jpg" width="500" height="500"/></a></td>
</tr>
</table>


</br></br>
<HR>
<h1>Violin plots</h1>

<table class="vignette">
<tr align="center">
<td><a href="{{ out }}/violin1.jpg"><img class="gpt"  src="{{ out }}/violin1.jpg" width="500" height="500"/></a></td>
<td><a href="{{ out }}/violin2.jpg"><img class="gpt"  src="{{ out }}/violin2.jpg" width="500" height="500"/></a></td>
<td><a href="{{ out }}/violin3.jpg"><img class="gpt"  src="{{ out }}/violin3.jpg" width="500" height="500"/></a></td>
</tr>
</table>


</br></br>
<HR>
<h1>PCA</h1>
Indicate if variables are significantly responsible for variability of one of the Principal Component.</br>
If your variable of interest is not the more significant for PC1 and followings you should concidere to use batch correction.
<table>
<tr align="center">
<td></td>
<td><select class="toogle" onchange="showhidepq(this);"><option value=pvalue>pvalue</option><option value=qvalue>qvalue</option></select></td>
</tr>
<tr>
<td><a href="{{ out }}/pca_contributions_raw.jpg"><img class="gpt"  src="{{ out }}/pca_contributions_raw.jpg" width="500" height="500"/></a></td>
<td><div class="pvalue">{{ pvalue }}</div> <div class="qvalue">{{ qvalue }}</div></td>
<tr>
</table>


</body>
<script type="text/javascript">

/*to color table of pvalues*/
function colortables() {

    var qtables = document.getElementsByTagName("table")
    for (var qt = 0; qt < qtables.length; qt++){
	var cellules = qtables[qt].getElementsByTagName("td");
	for (var i = 0; i < cellules.length; i++) {
		if ( cellules[i].textContent < 0.05) {
		    cellules[i].style.backgroundColor = "#b0e2ff";
		}
		if ( cellules[i].textContent < 0.01) {
		    cellules[i].style.backgroundColor = "#8db6cd";
		}
		if ( cellules[i].textContent < 0.001) {
		    cellules[i].style.backgroundColor = "#87ceff";
		}
		if ( cellules[i].textContent < 0.0001) {
		    cellules[i].style.backgroundColor = "#00bfff";
		}
		if ( cellules[i].textContent < 0.00001) {
		    cellules[i].style.backgroundColor = "#3a5fcd";
		}
	}
	}

}

/*group selection*/
var selection="Sample_Group";

function setPicture(select) {
var images = document.getElementsByClassName("gpt");
for(i=0;i<images.length;i++)
{
images[i].src = images[i].src.replace(selection, select.options[select.selectedIndex].value)
}
selection=select.options[select.selectedIndex].value
}

/*toogle pvalue qvalue*/
function toggle(className, displayState){
    var elements = document.getElementsByClassName(className)

    for (var i = 0; i < elements.length; i++){
        elements[i].style.display = displayState;
    }
}

function showhidepq(item){
    
    if (item.value === "pvalue"){
	toggle('pvalue', 'block')
	toggle('qvalue', 'none');
    } else {
	toggle('pvalue', 'none'); 
	toggle('qvalue', 'block'); 
    }

    var elements = document.getElementsByClassName("toogle")
    for (var i =0; i< elements.length; i++){
	elements[i].value=item.value;
    }

}

/*load*/
toggle('pvalue', 'block'); // Shows
toggle('qvalue', 'none'); // hides
colortables();

</script>

</html>
