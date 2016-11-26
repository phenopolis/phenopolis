Graphviz *dot* in your browser
==============================
Bower component `graphviz-d3-renderer` renders [Graphviz](http://graphviz.org) source in the browser with [d3.js](https://github.com/mbostock/d3). Check it out on [Graphviz fiddling website](http://graphviz.it/).

[![Build Status](https://travis-ci.org/mstefaniuk/graph-viz-d3-js.svg?branch=master)](https://travis-ci.org/mstefaniuk/graph-viz-d3-js)
[![Coverage Status](https://coveralls.io/repos/mstefaniuk/graph-viz-d3-js/badge.svg?branch=master)](https://coveralls.io/r/mstefaniuk/graph-viz-d3-js?branch=master)
![Forks](https://img.shields.io/github/forks/mstefaniuk/graph-viz-d3-js.svg)
![Stars](https://img.shields.io/github/stars/mstefaniuk/graph-viz-d3-js.svg)

Contents
--------
* `dot` parser with lax mode to verify Graphviz input
* `dot` mode for ACE editor
* stage data renderer with `d3.js`

Design
------
DOT parser is written in [PEG.js](https://github.com/dmajda/pegjs) has lax mode to parse source to the end with all errors. [Graphviz](http://graphviz.org) is embedded in browser using [viz.js](https://github.com/mdaines/viz.js).
Instead of using SVG directly it uses `xdot` format and parses it. Data structure of the output is drawn using
[d3.js](https://github.com/mbostock/d3) with animations during rendering.

Usage
-----
To include it in your project simply use `Bower`:
```
bower install graphviz-d3-renderer --save
```
Note that it needs `require.js` to work. Before loading proper paths should be defined for renderer and its dependecies (`d3.js` and `worker` for `require.js` plugin):
```javascript
requirejs.config({
	//By default load any module IDs from js/lib
	baseUrl: 'js',
	//except, if the module ID starts with "app",
	//load it from the js/app directory. paths
	//config is relative to the baseUrl, and
	//never includes a ".js" extension since
	//the paths config could be for a directory.
	paths: {
		d3: '/bower_components/d3/d3',
		"dot-checker": '/bower_components/graphviz-d3-renderer/dist/dot-checker',
		"layout-worker": '/bower_components/graphviz-d3-renderer/dist/layout-worker',
		worker: '/bower_components/requirejs-web-workers/src/worker',
		renderer: '/bower_components/graphviz-d3-renderer/dist/renderer'
	}
});
```
Then you can inject it into you app:
```javascript
require(["renderer"],
  function (renderer) {

  dotSource = 'digraph xyz ...';
  // initialize svg stage
  renderer.init("#graph");

  // update stage with new dot source
  renderer.render(dotSource);
});
```
Now you can even zoom / drag your graph

```javascript
require(["renderer"],
	function (renderer) {
		dotSource = 'digraph xyz ...';
		// initialize svg stage. Have to get a return value from renderer.init 
		//   to properly reset the image.
	    zoomFunc = renderer.init({element:"#graph", extend:[0.1, 10]});

	    // update stage with new dot source
	    renderer.render(dotSource);
	    
	    // for saving the image, 
	    $('#copy-button').on('click', function(){
		    $('#copy-div').html(renderer.getImage({reset:true, zoomFunc:zoomFunc}));
	    });	  
	    
	    // if do not need to reset the image before saving the image
	    $('#copy-button').on('click', function(){
		    $('#copy-div').html(renderer.getImage());
	    });	
});  
```
Even more powerful, you can add customised labelAttributer to change the position of the texts, and a callback function to add some fancy effects! (like a tooltip?)

```javascript
require(["renderer"],
	function (renderer) {
		$.get('my.dot', function(data){
			
			// initialize svg stage
			zoomFunc = renderer.init({element:"#graph", extend:[0.1, 10]});
			// update stage with new dot source
			renderer.render({source:data, labelAttributer: myLabelAttributer, callBack:callback});
			$('#copy-button').on('click', function(){
				$('#copy-div').html(renderer.getImage({reset:true, zoomFunc:zoomFunc}));
			});
			//check_element = $('#graph > svg > g > g');
			//check_load(check_element, callback);
		});	  
});
//Hide the tooltip when the mouse moves away
function removeTooltip() {

	//Fade out the circle to normal opacity
	d3.select(this).select('path')
		.attr('style', 'stroke:rgba(255,255,254,1);fill:rgba(176,224,230,1)');

	//Hide tooltip
	$('.popover').each(function() {
		$(this).remove();
	}); 

}//function removeTooltip

//Show the tooltip on the hovered over slice
function showTooltip(d) {

	//Define and show the tooltip
	$(this).popover({
		placement: 'auto top',
		container: '#graph',
		trigger: 'manual',
		position: 'fixed',
		html : true,
		title: function() { return d.labels[d.labels.length - 1].text; },
		content: function() {
			var content = '';
			$.each(d.labels, function(i, e){
				if(i){
					if (i < d.labels.length - 1){
						content = content + "<br/><span style='font-size: 11px;'>"
							+  e.text + "</span>";
					}
				}else{
					content = "<span style='font-size: 11px;'>"
					+  e.text + "</span>"
				}	
			}); 
			return content; 
		}
	});
	$(this).popover('show');
	
	//highlight
	d3.select(this).select('path')
		.attr("style", 'stroke:rgba(255,255,254,1);fill:#58FAF4');
			
}//function showTooltip

// shrink the texts, and add a tooltip
callback = function (){
	svg = d3.select('#graph').select('svg');
	all_g = svg.select('g').select('g').selectAll('g.node');
	all_g.selectAll('text')
		.style('font-size', '7px');
	all_g
		.on('mouseover', showTooltip)
		.on('mouseout', removeTooltip);
	$('.node').on('click',function() {
		// get the HPO id
		var this_node = $(this);
		var id = this_node.children().last().html();
		window.location.href = "http://localhost:3000/" + id;
	});
};
// custom labelAttributer
var myLabelAttributer = function() {
	  this
		.attr("x", function (d) {
		  return d.x;
		})
		.attr('y', function(d,i){
				var pData = d3.select(this.parentNode).datum().labels;
				//console.log(pData.length);
				return -d.y - 7*(i - pData.length/2 + 1);
		})
		.text(function (d) {
		  return d.text;
		});

	  this.each(function(d) {
		var self = d3.select(this);
		d.style.map(function(e) {
		  switch (e.key) {
			case "stroke":
			  return {key: "color", value: e.value};
			case "font-size":
			  return {key: e.key, value: e.value + "px"};
			default:
			  return e;
		  }
		}).forEach(function(e) {
		  self.style(e.key, e.value);
		});
	  });
};
```
And in CSS, you may add,

```CSS
.node:hover {
	cursor:pointer;
}
#graph path {
	-webkit-transition: fill 0.8s; /* Safari */
	transition: fill 0.8s;
}
.popover {
	text-align: center;
```
Roadmap
-------
* Test suite using Graphviz gallery examples (50% done)
* Improve animations with path tweening and concatenation of arrow heads with arrow arcs
* Custom `viz.js` compile with `xdot` output only to optimize size

License
-------
Currently project is available on LGPL so you can use it unmodified in free or commercial projects. If you add improvements
to it you must share them.
