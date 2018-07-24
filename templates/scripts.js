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

