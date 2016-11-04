// Some functions to display HPO terms
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
            container: $(this).parent().parent().parent().parent(),
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

// shrink the fucking texts
callback = function (svg){
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
            window.location.href = "https://uclex.cs.ucl.ac.uk/hpo/" + id;
            });
}
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
