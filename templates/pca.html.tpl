</br></br>
<HR>
<h1> PCA {{ correction }} </h1>
process : {{ process_id }}</br>
Indicate if variables are significantly responsible for variability of one of the Principal Component.</br>
If your variable of interest is not the more significant for PC1 and followings you should concidere to use batch correction.
<table>
<tr align="center">
<td></td>
<td><select class="toogle" onchange="showhidepq(this);"><option value=pvalue>pvalue</option><option value=qvalue>qvalue</option></select></td>
</tr>
<tr>
<td><a href="pca_contributions_{{ process_id }}.jpg"><img class="gpt"  src="pca_contributions_{{ process_id }}.jpg" width="500" height="500"/></a></td>
<td><div class="pvalue">{{ pvalue }}</div> <div class="qvalue">{{ qvalue }}</div></td>
<tr>
</table>
