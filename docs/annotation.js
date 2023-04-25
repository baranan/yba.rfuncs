//var button = document.getElementById('add-annotation-button');

  var button = document.createElement('button');
  button.id = 'add-annotation-button';
  button.innerHTML = 'Add Annotation';
  
button.addEventListener('click', function() {
	var annotationDiv = document.createElement('div');
	annotationDiv.className = 'ybaannotation';
	var annotationId = 'annotation-' + Math.random().toString(36).substr(2, 9); // generate a random string
	annotationDiv.setAttribute('data-annotation-id', annotationId);


	annotationDiv.style.left = '50%';
	annotationDiv.style.top = (window.innerHeight / 2 - annotationDiv.offsetHeight / 2 + window.pageYOffset) + 'px';
	//annotationDiv.style.height = '200px';
	annotationDiv.style.width = '240px';

	var p = document.createElement('p'); // create a new p element
    p.innerHTML = 'Click to edit annotation'; // set the initial text of the p element
    annotationDiv.appendChild(p);

    var textarea = document.createElement('textarea');
    textarea.style.display = 'none'; // hide the textarea initially
	textarea.style.width = '100%'; // set the width to 100%
	textarea.style.boxSizing = 'border-box'; // set the box-sizing property to border-box

	textarea.addEventListener('input', function() {
		this.style.height = 'auto';
		this.style.height = this.scrollHeight + 'px';
	});


    annotationDiv.appendChild(textarea);
    document.body.appendChild(annotationDiv);

    // Add blur event listener to the textarea to replace it with the p element
    textarea.addEventListener('blur', function() {
        p.innerHTML = this.value || 'Click to add annotation'; // set the text of the p element to the textarea value, or the default text if it's empty
        textarea.style.display = 'none';
        p.style.display = 'block';
    });

    // Add click event listener to the p element to replace it with the textarea
    p.addEventListener('click', function() {
        p.style.display = 'none';
        textarea.style.display = 'block';
        textarea.focus();
    });
	
	var isDragging = false;
	var startX, startY, startLeft, startTop;

	annotationDiv.addEventListener('mousedown', function(e) {
		isDragging = true;
		startX = e.clientX;
		startY = e.clientY;
		startLeft = parseInt(window.getComputedStyle(this).getPropertyValue('left'), 10);
		startTop = parseInt(window.getComputedStyle(this).getPropertyValue('top'), 10);
	});

	document.addEventListener('mousemove', function(e) {
		if (isDragging) {
			var dx = e.clientX - startX;
			var dy = e.clientY - startY;
			annotationDiv.style.left = startLeft + dx + 'px';
			annotationDiv.style.top = startTop + dy + 'px';
		}
	});

	document.addEventListener('mouseup', function(e) {
		isDragging = false;
	});
});

document.body.appendChild(button);


//var saveButton = document.getElementById('save-annotations-button');

var saveButton = document.createElement('button');
saveButton.id = 'save-annotations-button';
saveButton.innerHTML = 'Save annotated file';

saveButton.addEventListener('click', function() {saveAnnotations();});

function saveAnnotations() {
	var annotations = []; // initialize an empty array to store the annotations
	var annotationDivs = document.getElementsByClassName('ybaannotation'); // get all the annotationDiv elements
	for (var i = 0; i < annotationDivs.length; i++) {
		var annotationId = annotationDivs[i].getAttribute('data-annotation-id'); // get the unique identifier
		var annotationHtml = annotationDivs[i].outerHTML; // get the outerHTML
		// add the text from the textarea to the annotation HTML
		/*annotationHtml = '<div class="ybaannotation" style="left:'+annotationDivs[i].style.left + 
		' top:'+annotationDivs[i].style.top + ' width:'+annotationDivs[i].style.width + 
		'><p>' + annotationText + '</p></div>';*/
		annotations.push({ id: annotationId, html: annotationHtml}); // add the annotation to the array
	}
	
  // Get the existing HTML content of the page
  var pageHtml = document.documentElement.outerHTML;

  // Append the outerHTML of the annotation divs to the page HTML
  var annotationCss = '';
  for (var j = 0; j < annotations.length; j++) {
	pageHtml += annotations[j].html;
	annotationCss += '.ybaannotation[data-annotation-id="' + annotations[j].id + '"] {' +
	annotationDivs[j].style.cssText + '}\n';
  }

  // Create a new style element and append the CSS rules
  var style = document.createElement('style');
  style.type = 'text/css';
  style.appendChild(document.createTextNode(annotationCss));

  // Create a new Blob object with the updated HTML content
  var blob = new Blob([pageHtml], {type: "text/html;charset=utf-8"});

  // Create a link element and simulate a click on it to download the file
  var link = document.createElement("a");
  link.download = "annotations.html";
  link.href = URL.createObjectURL(blob);
  link.click();
}

document.body.appendChild(saveButton);
