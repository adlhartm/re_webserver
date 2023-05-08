var thickness = 2;
var thicknessHighlight = 4.0;

var seq_p = $("#prot_seq").val(),
    seq_r = $("#rna_seq").val(),
    type_p = $("#prot_selector > .btn.active").text().trim(),
    type_r = $("#rna_selector > .btn.active").text().trim(),
    wndw_p = 21,
    wndw_r = 63,
    offset = 0;

var window_p = $('#prot-window').slider()
                                .on("slide", function() {
                                    wndw_p = Number(this.value);
                                    wndw_r = 3*Number(this.value);
                                    updateDrawing();
                                })
                                .on("slideStop", function() {
                                    wndw_p = Number(this.value);
                                    wndw_r = 3*Number(this.value);
                                    updateDrawing();
                                })
                                .data('slider')
var window_r = $('#rna-window').slider()
                                .on("slide", function() {
                                    wndw_r = Number(this.value);
                                    updateDrawing();
                                })
                                .on("slideStop", function() {
                                    wndw_r = Number(this.value);
                                    updateDrawing();
                                })
                                .data('slider')
var offset_slider = $('#shift').slider()
                                .on("slide", function() {
                                    offset = Number(this.value);
                                    updateDrawing();
                                })
                                .on("slideStop", function() {
                                    offset = Number(this.value);
                                    updateDrawing();
                                })
                                .data('slider')

var profile_p = protein_calculation(seq_p,type_p,wndw_p),
    profile_r = rna_calculation(seq_r,type_r,wndw_r);

offset_slider.setAttribute('min', -profile_r.length+3);
offset_slider.setAttribute('max', profile_p.length-3);
offset_slider.setAttribute('focus', true);

var x1_min = 0,
    x1_max = seq_p.length,
    y1_min = d3.min(profile_p, function(d) {return d.y}),
    y1_max = d3.max(profile_p, function(d) {return d.y}),
    y2_min = d3.min(profile_r, function(d) {return d.y}),
    y2_max = d3.max(profile_r, function(d) {return d.y}),
    padding1 = (y1_max - y1_min)*0.1,
    padding2 = (y2_min - y2_max)*0.1;

var svg = d3.select(".canvas"),
    margin = {top: 30, right: 70, bottom: 30, left: 70},
    width = +svg.attr("width") - margin.left - margin.right,
    height = +svg.attr("height") - margin.top - margin.bottom,
    g = svg.append("g").attr("transform", "translate(" + margin.left +"," + margin.top +")");

var text = g.append('text')
            .classed("cor", true)
            .attr('font-size', '15px')
            .attr('font-weight', 'bold')
            .attr("transform", "translate(" + width*0.01 + "," + height*0.05 + ")")
            .text('R =')

var mask = svg.append("defs").append("clipPath")
                .attr("id", "clip")
              .append("rect")
                .attr("width", width)
                .attr("height", height)

var x1 = d3.scaleLinear()
    .domain([x1_min, x1_max])
    .rangeRound([0,width]);

var y1 = d3.scaleLinear()
    .domain([y1_min - padding1, y1_max + padding1])
    .rangeRound([height, 0]);

var y2 = d3.scaleLinear()
    .domain([y2_max - padding2, y2_min + padding2])
    .rangeRound([height, 0]);

var xAxis = d3.axisBottom(x1);
var yAxis1 = d3.axisLeft(y1);
var yAxis2 = d3.axisRight(y2);

var gX1 = g.append('g')
           .classed("axisBlue", true)
           .classed("x1", true)
           .classed("axis", true)
           .attr("transform", "translate(0," + height + ")")
           .call(xAxis)
var gY1 = g.append("g")
           .classed("axisBlue", true)
           .classed("yAxis", true)
           .classed("axis", true)
           .call(yAxis1)
var gY2 = g.append("g")
           .classed("axisBlue", true)
           .classed("yAxis", true)
           .classed("axis", true)
           .attr("transform", "translate(" + width + ",0)")
           .call(yAxis2)

var protLine = d3.line()
    .x(function(d) {return x1(d.x); })
    .y(function(d) {return y1(d.y); });

var rnaLine  = d3.line()
    .x(function(d) {return x1(d.x); })
    .y(function(d) {return y2(d.y); });

var rL = g.append("path")
        .classed("rnaline", true)
        .datum(profile_r)
        .attr("clip-path", "url(#clip)")
        .attr("fill", "none")
        .attr("stroke", "tomato")
        .attr("stroke-linejoin", "round")
        .attr("stroke-linecap", "round")
        .attr("stroke-width", thickness)
        .attr("d", rnaLine)
        // .on("mouseover", mouseover)
        // .on("mouseout", mouseout);

var pL = g.append("path")
        .classed("protline", true)
        .datum(profile_p)
        .attr("clip-path", "url(#clip)")
        .attr("fill", "none")
        .attr("stroke", "steelblue")
        .attr("stroke-linejoin", "round")
        .attr("stroke-linecap", "round")
        .attr("stroke-width", thickness)
        .attr("d", protLine)
        // .on("mouseover", mouseover)
        // .on("mouseout", mouseout);

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

// g.raise()

limits = findOverlap(profile_p, profile_r);
pSlice = profile_p.slice(limits[0], limits[1]);
rSlice = profile_r.slice(limits[2], limits[3]);
cor = spearson.round(spearson.correlation.pearson(yVal(pSlice), yVal(rSlice), true),2).toString();

text.text('R = ' + cor);
text.raise();
document.getElementById('corr_output').innerHTML = cor;

// function mouseover(d, i) {
//     console.log("over");
//     d3.select(this).attr("stroke-width", thicknessHighlight);
// };

// function mouseout(d, i) {
//     console.log("out");
//     d3.select(this).attr("stroke-width", thickness);
// };

function zoomed() {
    var xz = d3.event.transform.rescaleX(x1);
    // var xy1 = d3.event.transform.rescaleY(y1);
    // var xy2 = d3.event.transform.rescaleY(y2);

    gX1.call(xAxis.scale(xz));
    // gY1.call(yAxis1.scale(xy1));
    // gY2.call(yAxis2.scale(xy2));

    rnaLine
        .x(function(d) {return xz(d.x); })
        // .y(function(d) {return xy2(d.y); })
    protLine
        .x(function(d) {return xz(d.x); })
        // .y(function(d) {return xy1(d.y); })
    d3.selectAll(".protline").attr("d", protLine)
    // pL.attr("d", protLine)
    d3.selectAll(".rnaline").attr("d", rnaLine)

}

function updateProteinSequence() {
    seq_p = $("#prot_seq").val();

    // offset_slider.setAttribute('min', -profile_r.length+3)
    // offset_slider.setAttribute('max', profile_p.length-3)
    // offset_slider.setAttribute('focus', true)

    updateDrawing();
}

function updateRNASequence() {
    seq_r = $("#rna_seq").val();

    // offset_slider.setAttribute('min', -profile_r.length+3)
    // offset_slider.setAttribute('max', profile_p.length-3)
    // offset_slider.setAttribute('focus', true)
    
    updateDrawing();
}

function updateDrawing() {

    // wndw_p = window_p.getValue();
    // wndw_r = window_r.getValue();

    // offset = offset_slider.getValue();

    profile_r = rna_calculation(seq_r,type_r,wndw_r);
    profile_p = protein_calculation(seq_p,type_p,wndw_p);

    offset_slider.setAttribute('min', -profile_r.length+3);
    offset_slider.setAttribute('max', profile_p.length-3);
    offset_slider.setAttribute('focus', true);

    for (i=0;i<profile_r.length; i++) {
        profile_r[i].x = profile_r[i].x+offset
    }

    x1_max = seq_p.length;
    y1_min = d3.min(profile_p, function(d) {return d.y});
    y1_max = d3.max(profile_p, function(d) {return d.y});
    y2_min = d3.min(profile_r, function(d) {return d.y});
    y2_max = d3.max(profile_r, function(d) {return d.y});
    padding1 = (y1_max - y1_min)*0.1;
    padding2 = (y2_min - y2_max)*0.1;

    x1.domain([0,x1_max]);
    y1.domain([y1_min - padding1, y1_max + padding1]);
    y2.domain([y2_max - padding2, y2_min + padding2]);
    gX1.call(xAxis.scale(x1));
    gY1.call(yAxis1.scale(y1));
    gY2.call(yAxis2.scale(y2));

    rL.datum(profile_r)
      .attr("d", rnaLine);
    pL.datum(profile_p)
      .attr("d", protLine);

    limits = findOverlap(profile_p, profile_r);
    pSlice = profile_p.slice(limits[0], limits[1]);
    rSlice = profile_r.slice(limits[2], limits[3]);
    cor = spearson.round(spearson.correlation.pearson(yVal(pSlice), yVal(rSlice), true),2).toString();

    text.text('R = ' + cor);
    text.raise();
    document.getElementById('corr_output').innerHTML = cor;
}

$('.btn-draw-p').on("click", function () {
    type_p = $(this).text().trim();
    updateDrawing();
});
$('.btn-draw-r').on("click", function () {
    type_r = $(this).text().trim();
    updateDrawing();
});
$('#btn-more').on("click", function() {
    offset += 1;
    offset_slider.setValue(offset);
    updateDrawing();
});
$('#btn-less').on("click", function() {
    offset -= 1;
    offset_slider.setValue(offset);
    updateDrawing();
});


function parseUniprot() {
    var xhttp = new XMLHttpRequest();
    var parser = new DOMParser();
    var uniprotID = $(".UniprotID").val();
    var protName = "";

    $(".alert_box").empty()

    if (uniprotID.length == 6) {
        xhttp.open("GET", "https://www.ebi.ac.uk/proteins/api/proteins/" + uniprotID, false);
        xhttp.setRequestHeader("Content-type", "application/json");
        xhttp.send();
        if (xhttp.status == "200") {
            
            response = JSON.parse(xhttp.responseText);
            $("#prot_seq").val(response.sequence.sequence)
            protseq = response.sequence.sequence
            
            var rnaTargets = [];

            try {
                protName = response.protein.recommendedName.fullName.value;
            } catch (err) {
                protName = "*Name not found*";
            }

            response["dbReferences"].forEach(function (obj) {
                if (obj.type == 'EMBL' && obj["properties"].hasOwnProperty("protein sequence ID")) {
                    rnaTargets.push(obj.properties["protein sequence ID"].split(".")[0])
                }
            });

            $("#rna_seq").val('');
            
            if (rnaTargets.length == 0) {
                $(".alert_box").append("<div class='alert alert-warning alert-dismissable'><button type='button' class='close' data-dismiss='alert' aria-hidden='true'>&times;</button><strong>Warning!</strong> successfully loaded data for " + protName + " but could not retrieve mRNA information!</div>");    
            } else {

                for (let i=0; i<rnaTargets.length; i++) {
                    var mrna_request = "https://www.ebi.ac.uk/ena/data/view/" + rnaTargets[i] + "&display=xml"
                    xhttp.open("GET", mrna_request, false);
                    xhttp.setRequestHeader("Content-type", "text/xml");
                    xhttp.send();

                    mRNA_response = parser.parseFromString(xhttp.responseText, "text/xml");
                    mRNA = mRNA_response.getElementsByTagName("sequence")[0].childNodes[0].nodeValue;
                    mRNA = mRNA.replace(/\s/g,'');
                    mRNA = mRNA.toUpperCase();
                    mRNA = mRNA.replace(/T/g,'U');
                    
                    if (mRNA.length == 3*protseq.length + 3) {
                        $("#rna_seq").val(mRNA);
                        console.log(rnaTargets[i]);
                        $(".alert_box").append("<div class='alert alert-success alert-dismissable'><button type='button' class='close' data-dismiss='alert' aria-hidden='true'>&times;</button><strong>Success</strong> loading data for " + protName + "</div>");
                        break;

                    } else {
                        console.log(mRNA);                      
                        console.log(rnaTargets[i]);
                    };
                    $(".alert_box").append("<div class='alert alert-warning alert-dismissable'><button type='button' class='close' data-dismiss='alert' aria-hidden='true'>&times;</button><strong>Warning!</strong> successfully loaded data for " + protName + " but could not retrieve mRNA information!</div>");

                };
            };

        } else {
            $(".alert_box").append("<div class='alert alert-danger alert-dismissable'><button type='button' class='close' data-dismiss='alert' aria-hidden='true'>&times;</button><strong>Error</strong> could not locate given ID</div>");
        
        }
        seq_p = $("#prot_seq").val();
        seq_r = $("#rna_seq").val();

        updateDrawing();

    } else {
        $(".alert_box").append("<div class='alert alert-danger alert-dismissable'><button type='button' class='close' data-dismiss='alert' aria-hidden='true'>&times;</button><strong>Error</strong> probably not a UniprotID.</div>");
    }

};