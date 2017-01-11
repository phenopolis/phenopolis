//modified from main.js, to display table in /individual/<pid>

//tooltip
$('.tip').tooltip({html:true});
//popover
$('.pop').popover({html:true,title:'loading...', trigger: 'hover manual'});

//check if an object is in an array
function include(arr,obj) {
    return (arr.indexOf(obj) != -1);
}

// Sending post via js
function post(path, params, method) {
    method = method || "post"; // Set method to post by default if not specified.

    // The rest of this code assumes you are not using a library.
    // It can be made less wordy if you use one.
    var form = document.createElement("form");
    form.setAttribute("method", method);
    form.setAttribute("action", path);

    for(var key in params) {
        if(params.hasOwnProperty(key)) {
            var hiddenField = document.createElement("input");
            hiddenField.setAttribute("type", "hidden");
            hiddenField.setAttribute("name", key);
            hiddenField.setAttribute("value", params[key]);

            form.appendChild(hiddenField);
         }
    }

    document.body.appendChild(form);
    form.submit();
}

// display the result from /tools/batch_pubmed
function displayPubmedTable(data) {
	//max score to match with color
	var maxi = 50;
	//initialise colour
	var color = d3.scale.quantize().range([
		"rgb(255,255,255)",
		"rgb(255,255,230)",
        "rgb(255,255,204)",
        "rgb(255,237,160)",
        "rgb(254,217,118)",
        "rgb(254,178,76)",
        "rgb(253,141,60)",
        "rgb(252,78,42)",
        "rgb(227,26,28)"]);
    color.domain([0, maxi]); //max score is maxi

	var titles = data[0];
	
	var col_width = new Array(data[0].length + 1).join('80 ').split(' ').map(parseFloat); // define width of columns.
	col_width[0] = col_width[6] = 100;
	col_width[2] = 180;
	var table = d3.select('#table-content').append('table')
		.attr('class', 'table table-bordered sortable')
		//.style('width', 'auto')
		.style('table-layout', 'fixed');
		
	table.append('colgroup')
		.selectAll('col')
		.data(col_width)
		.enter()
		.append('col')
		.attr('style', function(d){ return 'width: ' + d + 'px;';});
	
	// add an empty element to the beginning of each row. this is to add a button to delete a given row
	titles.unshift('Quality');
	// create table header
    var th = table.append('thead').append('tr')
        .selectAll('th')
        .data(titles).enter()
        .append('th')
        .attr('class', 'list-title')
        .style('overflow', 'hidden')
        .style('text-overlow', 'ellipsis')
        .text(function(d) {
            d = d.replace(/&#46;/g, '.');
            return d
        });
	$('th').each(function(){
		$(this).prepend('<div class="some-handle"></div>');
	});
	
	var results = data[1];
	results.sort(sortPubmed);
	
	// create table body
    var tr = table.append('tbody')
        .selectAll('tr')
        .data(results).enter()
        .append('tr')
        .selectAll('td')
        .data(function(row, i) {
            return titles.map(function(c) {
                // compute cell values for this specific row
                if(c == 0){
                	// 1st cell, return nothing
                	return '';
                }
                
                var cell = row[c];
                if(c == 'ref(pubmedID)'){
                	if (cell.results == 'masked'){
                		//masked by the user-provided mask-gene list
                		return {results: 'masked', known: cell.known};
                	}
                	var results = [];
                	//sort results according to score.
                	cell.results.sort(function(a,b){
                		return (b.score - a.score);
                	});
                	cell.results.map(function(r){
                		// title and abstract can't have single quote/</>. replace it with &#39;/&#60;/&#62;
                		var sq = new RegExp("'", 'g');
                		var lt = new RegExp('<', 'g');
                		var gt = new RegExp('>', 'g');
                		if (r.title){
                			r.title = r.title.replace(sq, "&#39;");
                			r.title = r.title.replace(lt, "&#60;");
                			r.title = r.title.replace(gt, "&#62;");
                		} else {
                			r.title = 'NA';
                		}
                		if (r.abstract){
                			r.abstract = r.abstract.replace(sq, "&#39;");
                			r.abstract = r.abstract.replace(lt, "&#60;");
                			r.abstract = r.abstract.replace(gt, "&#62;");
                		} else {
                			r.absract = 'NA';
                		}
                		
                		
                		results.push("<span class='pop pointer' data-toggle='popover' title='<a href="
                						+ '"http://www.ncbi.nlm.nih.gov/pubmed/'
                						+ r.id
                						+ '" target="_blank" >' 
                						+ r.title
                						+ "</a>' data-content='"
                						+ r.abstract
                						+ "'>"
                						+ r.id
                						+ "</span>");
                	});
                	return {results: results.join(', '), known: cell.known};
                } else {
                	// need to get rid of the special characters again...
                	if (cell && cell.typeof == 'string'){
						var sq = new RegExp("'", 'g');
						var dq = new RegExp('"', 'g');
						var lt = new RegExp('<', 'g');
						var gt = new RegExp('>', 'g');
						cell = cell.replace(sq, "&#39;");
						cell = cell.replace(dq, "&#34;");
						cell = cell.replace(lt, "&#60;");
						cell = cell.replace(gt, "&#62;");
                	}
                	// return
                	return cell;
                }
            });
        }).enter()
        .append('td')
        .style('overflow', 'hidden')
        .style('text-overflow', 'ellipsis')
        .style('color', function(d, i){
        	// add color to known genes
        	if (titles[i] == 'HUGO'){
        		var pData = d3.select(this.parentNode).datum();
        		if (pData['ref(pubmedID)'].known == 1){
        			return 'green';
        		} else {
        			return 'black';
        		}
        	}
        })
        .style('background-color', function(d, i){
        	if (titles[i] == 'pubmed_score'){
				
				if (d > maxi){
					//color top maxi
					d = maxi;
				}
				return color(d);
        	} else {
        		return '';
        	}
        })
        .attr('title', function(d, i){ 
        		if (titles[i] == 'Samples' || titles[i] == 'HUGO' || titles[i] == 'Func' || titles[i] == 'ExonicFunc' || titles[i] == 'Description' || titles[i] == 'AAChange' || titles[i] == 'Gene' || titles[i] == 'Start' || titles[i] == 'End'){
        			return d;
        		} else {
        			return '';
        		}
        	})
        .html(function(d, i) {
        	if (i == 0){
        		// check filter first. Pass = green button, Fail = red, and dubious = yellow
        		var pData = d3.select(this.parentNode).datum();
        		var filter = ('FILTER' in pData) ? pData.FILTER : 'PASS';
        		var col = '';
        		if (filter == 'PASS'){
        			col = 'success'; //green
        		} else if (filter == 'FAIL') {
        			col = 'danger'; //red
        		} else {
        			col = 'warning'; //yellow
        		}
        		
        		return '<button type="button" class="btn btn-' + col + ' del-row">Del row</button> ';
        	} else if (titles[i] == 'ref(pubmedID)'){
        		return d.results;
            } else if (titles[i] == 'HUGO'){
                // add hyperlink if retnet
                var pData = d3.select(this.parentNode).datum();
                if ('disease' in pData['ref(pubmedID)']){
                    var return_string = "<span class='pop pointer' data-toggle='popover' title='<a href="
                            + '"https://sph.uth.edu/Retnet/disease.htm" target="_blank" >' 
                            + 'A RetNet gene'
                            + "</a>' data-content='";
                    if (pData['ref(pubmedID)'].omim){
                        // add omim
                        return_string += '<p style="color:black;margin-bottom:0"><b>OMIM:</b> ';
                        pData['ref(pubmedID)'].omim.map(function(o){
                            return_string += '<a href="http://omim.org/entry/' + o + '">'
                                + o + '</a> ';
                        });
                    }
                    return_string += '</p><hr/><p style="color:black;margin-top:0"><b>Disease:</b><br />'
                        + pData['ref(pubmedID)'].disease
                        + "</p>'>" + '<span class="retnet-mode"><b>'
                        + pData['ref(pubmedID)'].mode.toUpperCase()
                        + '</b></span>' + d + '</span>';
                    return return_string;
                } else {
                    return d;
                }
        	} else {
        		return d;
        	}
        });
	
	$('#table-content').show();
	$('#tablehead-div').show();
	
	//populate old-name field for renaming
	$('#old-name').val(data[3]);
	
	//popover
	$("[data-toggle='popover']").popover({html:true, trigger: 'hover manual'});
	//sortable
	$("table").tablesorter();
	// delete row
	$(".del-row").on('click', function(event) {
		$(this).parent().parent().remove();
	});
	// export results to excel
	$('#export').on('click', function(){
		exportTableToCSV.apply(this, [$('#table-content'), 'export.csv']);
	});
	
}

function sortPubmed(a, b) {
	return (b['ref(pubmedID)'].total_score - a['ref(pubmedID)'].total_score);
}

function exportTableToCSV($table, filename) {
	var $rows = $table.find('tr:has(th), tr:has(td)'),
		// Temporary delimiter characters unlikely to be typed by keyboard
		// This is to avoid accidentally splitting the actual contents
		tmpColDelim = String.fromCharCode(11), // vertical tab character
		tmpRowDelim = String.fromCharCode(0), // null character

		// actual delimiter characters for CSV format
		colDelim = '","',
		rowDelim = '"\n"',

		// Grab text from table into CSV formatted string
		csv = '"' + $rows.map(function (i, row) {
			var $row = $(row),
				$cols = $row.find('th, td');

			return $cols.map(function (j, col) {
				var $col = $(col),
					text = $col.text();

				return text.replace(/"/g, '""').replace(/\n/g,'').replace(/^ +/g,'').replace(/\w {2,}/g, ''); // escape double quotes

			}).get().join(tmpColDelim);

		}).get().join(tmpRowDelim)
			.split(tmpRowDelim).join(rowDelim)
			.split(tmpColDelim).join(colDelim) + '"',
		// Data URI
		csvData = 'data:application/csv;charset=utf-8,' + encodeURIComponent(csv);
	

	$(this)
		.attr({
		'download': filename,
			'href': csvData,
			'target': '_blank'
	});
}

// get hash value by key pattern
function getValueByFilter (obj, filter) {
	var result = [];

    Object.keys(obj).map(function(d) {
        if (!filter || filter.test(d)) {
            result[result.length] = obj[d];
        }
    });
    return result;
}
