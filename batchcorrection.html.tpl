<html>
<head>
</head>
<body>

</br></br>
<HR>
<h1> PCA Batch correction : {{ correction }} with model {{ formula1 }} </h1>
process : {{ process_id }}</br>
Indicate if variables are significantly responsible for variability of one of the Principal Component.</br>
If your variable of interest is not the more significant for PC1 and followings you should concidere to use batch correction.
<b>{{ comment }}</b>
<table>
<tr align="center">
<td></td>
<td><select class="toogle" onchange="showhidepq(this);"><option value=pvalue>pvalue</option><option value=qvalue>qvalue</option></select></td>
</tr>
<tr>
<td><a href="{{ out }}/pca_contributions{{ process_id }}.jpg"><img class="gpt"  src="{{ out }}/pca_contributions{{ process_id }}.jpg" width="500" height="500"/></a></td>
<td><div class="pvalue">{{ pvalue }}</div> <div class="qvalue">{{ qvalue }}</div></td>
<tr>
</table>


</body>
<script>
toggle('pvalue', 'block'); // Shows
toggle('qvalue', 'none'); // hides
colortables();
</script>
</html>
