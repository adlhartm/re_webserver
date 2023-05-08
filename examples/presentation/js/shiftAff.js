var thickness = 2;
var thicknessHighlight = 4.0;

var seq_c = $("#coding_seq").val(),
    type_p = $("#prot_selector > .btn.active").text().trim(),
    wndw_s = 21,
    seq_1 = translate(seq_c, 0),
    seq_2 = translate(seq_c, 1),
    offset = 0;

var window_s = $('#shift-window').slider()
                                .on("slide", function() {
                                    wndw_s = Number(this.value);
                                    updateDrawingShift();
                                })
                                .on("slideStop", function() {
                                    wndw_s = Number(this.value);
                                    updateDrawingShift();
                                })
                                .data('slider');

var offset_slider = $('#shift').slider()
                                .on("slide", function() {
                                    offset = Number(this.value);
                                    updateDrawingShift();
                                })
                                .on("slideStop", function() {
                                    offset = Number(this.value);
                                    updateDrawingShift();
                                })
                                .data('slider');

$('.btn-draw-p').on("click", function () {
    type_p = $(this).text().trim();
    updateDrawingShift();
});
$('#btn-more').on("click", function() {
    offset_slider.setValue(offset_slider.getValue()+1);
    offset +=1;
    updateDrawingShift();
});
$('#btn-less').on("click", function() {
    offset_slider.setValue(offset_slider.getValue()-1);
    offset -=1;
    updateDrawingShift();
});

var profile_1 = protein_calculation(seq_1,type_p,wndw_s),
    profile_2 = protein_calculation(seq_2,type_p,wndw_s);

var x1_min = 0,
    x1_max = seq_1.length,
    x2_min = 0,
    x2_max = seq_2.length,
    y1_min = d3.min(profile_1, function(d) {return d.y}),
    y1_max = d3.max(profile_1, function(d) {return d.y}),
    y2_min = d3.min(profile_2, function(d) {return d.y}),
    y2_max = d3.max(profile_2, function(d) {return d.y}),
    padding = (d3.max([y1_max, y2_max]) - d3.min([y1_min, y2_min]))*0.1;

var svg = d3.select(".canvas"),
    margin = {top: 30, right: 70, bottom: 30, left: 70},
    width = +svg.attr("width") - margin.left - margin.right,
    height = +svg.attr("height") - margin.top - margin.bottom,
    g = svg.append("g").attr("transform", "translate(" + margin.left +"," + margin.top +")");

var text = g.append('text')
            .classed("cor", true)
            .attr('font-size', '15px')
            .attr('font-weight', 'bold')
            .attr("transform", "translate(" + width*0.047 + "," + height*0.05 + ")")
            .text('R =')

var rmstext = g.append('text')
            .classed("rmsd", true)
            .attr('font-size', '15px')
            .attr('font-weight', 'bold')
            .attr("transform", "translate(" + width*0.01 + "," + height*0.12 + ")")
            .text('RMSD =')

var mask = svg.append("defs").append("clipPath")
                .attr("id", "clip")
              .append("rect")
                .attr("width", width)
                .attr("height", height)

var x = d3.scaleLinear()
    .domain([d3.min([x1_min, x2_min]), d3.max([x1_max, x2_max])])
    .rangeRound([0,width]);

var y = d3.scaleLinear()
    .domain([d3.min([y1_min, y2_min])-padding, d3.max([y1_max, y2_max])+padding])
    .rangeRound([height, 0]);

var xAxis = d3.axisBottom(x);
var yAxis = d3.axisLeft(y);

var gX = g.append('g')
           .classed("axisBlue", true)
           .classed("x", true)
           .classed("axis", true)
           .attr("transform", "translate(0," + height + ")")
           .call(xAxis)
var gY = g.append("g")
           .classed("axisBlue", true)
           .classed("yAxis", true)
           .classed("axis", true)
           .call(yAxis)

var protLine = d3.line()
    .x(function(d) {return x(d.x); })
    .y(function(d) {return y(d.y); });

var rnaLine  = d3.line()
    .x(function(d) {return x(d.x); })
    .y(function(d) {return y(d.y); });

var rL = g.append("path")
        .classed("rnaline", true)
        .datum(profile_2)
        .attr("clip-path", "url(#clip)")
        .attr("fill", "none")
        .attr("stroke", "tomato")
        .attr("stroke-linejoin", "round")
        .attr("stroke-linecap", "round")
        .attr("stroke-width", thickness)
        .attr("d", rnaLine)

var pL = g.append("path")
        .classed("protline", true)
        .datum(profile_1)
        .attr("clip-path", "url(#clip)")
        .attr("fill", "none")
        .attr("stroke", "steelblue")
        .attr("stroke-linejoin", "round")
        .attr("stroke-linecap", "round")
        .attr("stroke-width", thickness)
        .attr("d", protLine)

limits = findOverlap(profile_1, profile_2);
pSlice = profile_1.slice(limits[0], limits[1]);
rSlice = profile_2.slice(limits[2], limits[3]);
cor = spearson.round(spearson.correlation.pearson(yVal(pSlice), yVal(rSlice), true),2).toString();

rmsd = spearson.round(RMSD(yVal(pSlice), yVal(rSlice)),4).toString();
rmstext.text('RMSD = ' + rmsd);
rmstext.raise();

text.text('R = ' + cor);
text.raise();
document.getElementById('corr_output').innerHTML = cor;

var zoom = d3.zoom()
    .scaleExtent([0.25, 8])
    .translateExtent([[-width, -height], [2 * width, 2 *height]])
    .on("zoom", zoomed);

var zoomRect = svg.append("rect")
    .attr("width", width)
    .attr("height", height)
    .attr("fill", "none")
    .attr("pointer-events", "all")
    .attr("transform", "translate(" + margin.left +"," + margin.top +")")
    .call(zoom);

function zoomed() {
    var xz = d3.event.transform.rescaleX(x);
    // var xy1 = d3.event.transform.rescaleY(y1);
    // var xy2 = d3.event.transform.rescaleY(y2);

    gX.call(xAxis.scale(xz));
    // gY1.call(yAxis1.scale(xy1));
    // gY2.call(yAxis2.scale(xy2));

    rnaLine
        .x(function(d) {return xz(d.x); })
        // .y(function(d) {return xy2(d.y); })
    protLine
        .x(function(d) {return xz(d.x); })
        // .y(function(d) {return xy1(d.y); })
    pL.attr("d", protLine)
    rL.attr("d", rnaLine)

}

function updateProteinSequences(shift1 = 0, shift2 = 1) {
    seq_c = $("#coding_seq").val();

    seq_1 = translate(seq_c, shift1);
    seq_2 = translate(seq_c, shift2);

    updateDrawingShift();    
}

function updateDrawingShift() {

    profile_1 = protein_calculation(seq_1,type_p,wndw_s);
    profile_2 = protein_calculation(seq_2,type_p,wndw_s);

    offset_slider.setAttribute('min', -profile_2.length+3)
    offset_slider.setAttribute('max', profile_1.length-3)
    offset_slider.setAttribute('focus', true)

    for (i=0;i<profile_2.length; i++) {
        profile_2[i].x = profile_2[i].x+offset
    }
    
    x1_max = seq_1.length;
    x2_max = seq_2.length;
    y1_min = d3.min(profile_1, function(d) {return d.y});
    y1_max = d3.max(profile_1, function(d) {return d.y});
    y2_min = d3.min(profile_2, function(d) {return d.y});
    y2_max = d3.max(profile_2, function(d) {return d.y});

    x.domain([0, d3.max([x1_max, x2_max])])
    y.domain([d3.min([y1_min, y2_min])-padding, d3.max([y1_max, y2_max])+padding])
    gX.call(xAxis.scale(x));
    gY.call(yAxis.scale(y));

    rL.datum(profile_2)
      .attr("d", rnaLine);
    pL.datum(profile_1)
      .attr("d", protLine);

    limits = findOverlap(profile_1, profile_2);
    pSlice = profile_1.slice(limits[0], limits[1]);
    rSlice = profile_2.slice(limits[2], limits[3]);
    cor = spearson.round(spearson.correlation.pearson(yVal(pSlice), yVal(rSlice), true),2).toString();

    rmsd = spearson.round(RMSD(yVal(pSlice), yVal(rSlice)),4).toString();
    rmstext.text('RMSD = ' + rmsd);
    rmstext.raise();

    text.text('R = ' + cor);
    text.raise();
    document.getElementById('corr_output').innerHTML = cor;
    // if (profile_p.length == profile_r.length) {
    //     text.
    //         text('R = ' + spearson.round(spearson.correlation.pearson(yVal(profile_p), yVal(profile_r), true),2).toString());
    //     document.getElementById('corr_output')
    //         .innerHTML = spearson.round(
    //                               spearson.correlation.pearson(yVal(profile_p), yVal(profile_r), true)
    //                               , 2);
    // } else {
    //     text.
    //         text("Err - Lengths don't match");
    //     document.getElementById('corr_output')
    //         .innerHTML = "Lengths don't match"
    // }

}




function parseENA() {
    var xhttp = new XMLHttpRequest();
    var parser = new DOMParser();
    var rnaID = $(".enaID").val();
    $(".alert_box").empty();

    try {
        

        var mrna_request = "https://www.ebi.ac.uk/ena/data/view/" + rnaID.split('.')[0] + "&display=xml"
        xhttp.open("GET", mrna_request, false);
        xhttp.setRequestHeader("Content-type", "text/xml");
        xhttp.send();
        
        mRNA_response = parser.parseFromString(xhttp.responseText, "text/xml");
        mRNA = mRNA_response.getElementsByTagName("sequence")[0].childNodes[0].nodeValue;
        mRNA = mRNA.replace(/\s/g,'');
        mRNA = mRNA.toUpperCase();
        mRNA = mRNA.replace(/T/g,'U');

        $("#coding_seq").val(mRNA);
    }
    catch (err) {
        $(".alert_box").empty();
        $(".alert_box").append("<div class='alert alert-danger alert-dismissable'><button type='button' class='close' data-dismiss='alert' aria-hidden='true'>&times;</button><strong>Error</strong> the provided ID could not be located!</div>");
    }

    updateDrawingShift();
};