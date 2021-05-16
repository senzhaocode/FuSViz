// modified from https://jsfiddle.net/rkkghrvL/7/

var element = document.getElementById('overviewmain');

var drawLines = function(event) {
  var divwidth = document.getElementById('overviewmain').offsetWidth;
  var rulerline = document.getElementById("dashline");
  var x = 0;
  var y = 0;
  var straightLine = element.querySelector('.straightLine');
  var hrLine = element.querySelector('.hrLine');
 
  if ( rulerline.checked == true ) {
      x = event.pageX; x = x - (divwidth + 80)/2;
      y = event.pageY;

      var slTrans = 'translate(' + x + 'px, 0px)';
      var hrTrans = 'translate(0px, ' + y + 'px)';
      if(!straightLine) {
         straightLine = document.createElement('div');
         straightLine.classList.add('straightLine');
         straightLine.style.height = "100%";
         straightLine.style.width = '0.5px';
         element.appendChild(straightLine);
      }
      straightLine.style.transform = slTrans;
  
      if(!hrLine) {
         hrLine = document.createElement('div');
         hrLine.style.height = "0.5px";
         hrLine.classList.add('hrLine');
         hrLine.style.width = '100%';
         element.appendChild(hrLine);
      }
      hrLine.style.transform = hrTrans;
  }
}

element.addEventListener('mousemove', function(event) {
   drawLines(event);
});

element.addEventListener('mousedown', function(event) {
   drawLines(event);   
});

element.addEventListener('mouseup', function(event) {
   drawLines(event);
});

element.addEventListener('mouseout', function(event) {
   drawLines(event);
});
