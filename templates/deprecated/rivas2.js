
function variant_colors(d) {
    if (d.category == 'lof_variant') {
        return '#cd2932';
    } else if (d.category == 'missense_variant') {
        return '#FF9966';
    } else if (d.category == 'synonymous_variant') {
        return '#228B22';
    }
    else{ 
         return "#d3d3d3";
    }
}

function change_coverage_chart(data, new_data, variant_data, _transcript, scale_type, metric, skip_utrs, container) {
    var coords = skip_utrs ? 'pos_coding_noutr' : 'pos_coding';
    var metricor = 'odds_ratio'; 
    var maxn = d3.max(data, function(d) { return d[metric] }) + 2;
    var maxnor = d3.max(data, function(d){ return d[metricor]}) + .3;
    var minor =  d3.min(data, function(d) { return d[metricor] }) - .3;
    var minor = minor < 0 ? 0 : minor; 
    var max_cov = (metric == '-log10pvalue') ? maxn : maxnor;
    var min_cov = (metric == '-log10pvalue') ? 0 : minor;
    var coding_coordinate_params = get_coding_coordinate_params(_transcript, skip_utrs);
    var chart_width;
    if (scale_type == 'overview') {
        chart_width = gene_chart_width;
    } else {
        chart_width = coding_coordinate_params.size*2;
    }

    var exon_x_scale = d3.scale.linear()
        .domain([0, coding_coordinate_params.size])
        .range([0, chart_width]);

    var svg = d3.select(container).select('#inner_svg')
        .attr("width", chart_width + gene_chart_margin.left + gene_chart_margin.right)
        .attr("height", gene_chart_height + gene_chart_margin.top + gene_chart_margin.bottom)
        .select('#inner_graph');

    var y = d3.scale.linear()
        .domain([min_cov, max_cov])
        .range([gene_chart_height, 0]);

    var area = d3.svg.area()
        .x( function(d) {
            return exon_x_scale(d[coords]);
        }).y0( function(d) {
            return gene_chart_height;
        }).y1( function(d) {
            return (metric in d) ? y(d[metric]) : gene_chart_height;
        });

    svg.selectAll("circle")
        .data(data)
        .transition()
        .duration(500)
        .attr("class", "dot")
<!-- stuff i added -->
        .attr("data-toggle", "tooltip")
        .attr('filter_status', function(d) {
            return d.flag;
        })
        .attr('category', function(d) {
            return d.category;
        })
        .attr("class", function(d) {
            return "track_variant " + d.category;
        })
 <!-- stuff i added -->

      .attr("cx", function(d) {
         if (d['pos_coding_noutr'] == undefined) {
                return -100;
         } else {
                return exon_x_scale(d[coords]);
         }
        })
        .attr("cy", function(d){
               return (metric in d) ? y(d[metric]) : gene_chart_height;
        })
        .attr('r', 4 )
        .style("fill", variant_colors);
;

    // plot exons
    var svg_outer = d3.select(container).select('#track_svg')
        .attr("width", chart_width + gene_chart_margin_lower.left + gene_chart_margin_lower.right)
        .attr("height", lower_gene_chart_height).select('#track');

    svg_outer.select('#boundary_line')
        .attr("x2", exon_x_scale(coding_coordinate_params.size));

    // plot exon rounded rects
    svg_outer.selectAll("rect")
        .data(_transcript.exons)
        .transition()
        .duration(500)
        .attr("x", function(d, i) { return exon_x_scale(get_coding_coordinate(_transcript, d.start, skip_utrs)); })
        .attr("width", function(d, i) {
            if (get_coding_coordinate(_transcript, d.start, skip_utrs) == -100) {
                return exon_x_scale(175);
            }
            return exon_x_scale(d.stop-d.start+1);
        })
        .attr("height", function(d, i) {
            if (d.feature_type == 'CDS') {
                return lower_gene_chart_height;
            } else {
                return lower_gene_chart_height/2;
            }
        });

    // plot variants
    svg_outer.selectAll("a")
        .data(variant_data)
        .transition()
        .duration(500)
        .selectAll('ellipse')
        .attr("cx", function(d) {
            if (d[coords] == undefined) {
                return -100;
            } else {
                return exon_x_scale(d[coords]);
            }
        });

    var yAxis = d3.svg.axis()
        .scale(y)
        .orient("left");

    svg.select(".y.axis")
        .transition()
        .duration(200)
        .call(yAxis);
}

function create_new_data(data, coords) {
    var data_object = {};
    $.each(data, function(i, d) {
        data_object[d[coords]] = d;
    });
    var new_data = [];
    for (var i = d3.min(data, function(d) { return d[coords] }); i < d3.max(data, function(d) { return d[coords] }); i++) {
        var x = {'has_coverage': false};
        x[coords] = i;
        //Check the previous base to see if this is the beginning of an exon

        if (i in data_object && !(i-1 in data_object)) {
            new_data.push(x);
        }
        //Check the previous base to see if this is the end of an exon
        if (!(i in data_object) && i-1 in data_object) {
            x[coords] = i-1;
            new_data.push(x);
        }
        if (i in data_object) {
            new_data.push(data_object[i]);
        } else {
            new_data.push(x);
        }
    }
    return new_data;
}

function coverage_sum(key) {
    var total = 0;
    $.map(window.coverage_stats, function(entry) {
        total += entry[key];
    });
    return (total/window.coverage_stats.length).toPrecision(4);
}

function refresh_links() {
    $("#coverage_plot_download").attr('href', set_plot_image('gene_plot_container', 0));
    $("#exon_plot_download").attr('href', set_plot_image('gene_plot_container', 1));
}


