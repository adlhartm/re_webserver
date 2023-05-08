function updateCorResults(corList) {
    for(corID in corList) {
        cor = calculate_pearson(corID, corList[corID]);
        cor = "R = " + cor;
        $("#"+corID+".corr_output").text(cor);
    }    
}

function findBestShift(id1, id2, inv=false) {
    realLength1 = data.seq[id1].profile.length + data.seq[id1].window - 1;
    realLength2 = data.seq[id2].profile.length + data.seq[id2].window - 1;
    data.seq[id1].shift = 0;
    data.seq[id2].shift = 0;
    drawingProfiles(data);

    if (realLength1 > realLength2) {
        correlations = []
        maxShift = realLength1 - realLength2;
        for (var i=0; i<maxShift; i++) {
            data.seq[id2].shift = i;
            drawingProfiles(data);
            correlations.push(Number(calculate_pearson(id1, id2)));
        }
        if (inv) { data.seq[id2].shift = correlations.indexOf(d3.min(correlations)); }
        else { data.seq[id2].shift = correlations.indexOf(d3.max(correlations)); }
        drawingProfiles(data);
        updateShifts(data);
    }
    else if (realLength1 < realLength2) {
        correlations = []
        maxShift = realLength2 - realLength1;
        for (var i=0; i<maxShift; i++) {
            data.seq[id1].shift = i;
            drawingProfiles(data);
            correlations.push(Number(calculate_pearson(id1, id2)));
        }
        if (inv) { data.seq[id1].shift = correlations.indexOf(d3.min(correlations)); }
        else { data.seq[id1].shift = correlations.indexOf(d3.max(correlations)); }        
        drawingProfiles(data);
        updateShifts(data);
    }
    else {
        console.log("No shift possible due to same length")
    }
}

function updateLegend(id) {
    $("#legend"+id+" .legend-text").text(getPropertyString(id));
    
}

function extractGaps(seq, wndw) {
//not used anymore
    var gaps = [];
    var gaplessSeq = "";

    for (var i = 0, len = seq.length; i < len; i++) {
        if (seq[i] == "-") {
            if (i >= (wndw-1)/2) {
                gaps.push(i-(wndw-1)/2);
            }
        } else {
            gaplessSeq += seq[i];
        }
    }

    return {gap: gaps, seq: gaplessSeq};
}
function extractGapsCodon(seq, wndw) {
//not used anymore
    var gaps = [];
    var gaplessSeq = "";

    for (var i = 0, len = seq.length; i < len-3; i+=3) {
        if (seq.slice(i,i+3) == "---") {
            if (i/3 >= (wndw-1)/2) {
                gaps.push(i/3-(wndw-1)/2);
            }
        } else {
            gaplessSeq += seq.slice(i,i+3);
        }
    }

    return {gap: gaps, seq: gaplessSeq};
}
function profile_calculation(seq, wndw, scaleID, type, shift, method, min_existing) {
    var seqArray = [];
    var weights = [];
    var smoothArray = [];
    var smoothWeights = [];

    if (scaleID == "" || scaleID == "no scale selected") {
        for (var i = 0, len = seq.length; i < len; i++) {
            seqArray.push(0);
            weights.push(1)
        };
    }
    else {
        for (var i = 0, len = seq.length; i < len; i++) {
            var val = scale[type][scaleID][seq[i]];
            var weight = Number(!isNaN(val));
            seqArray.push(val);
            weights.push(weight);
        };
    }
    switch(method) {
        case "boxcar":
            switch (type) {
                case "protein":
                    var tmp = windowAverage(seqArray, weights, wndw,1);
                    smoothArray = tmp.arr
                    smoothWeights = tmp.weight

                    var xVals = [Math.round((wndw-1)/2)+shift], x=Math.round((wndw-1)/2);
                    while(x<(smoothArray.length+Math.round((wndw-1)/2) ) ){x+=1;xVals.push(x+shift)};
                    break;

                case "dna":
                case "rna":
                    var tmp = windowAverage(seqArray, weights, wndw,3);
                    smoothArray = tmp.arr;
                    smoothWeights = tmp.weight;


                    var xVals = [((wndw-1)/2)+shift], x=(wndw-1)/2;
                    while(x<(smoothArray.length+Math.round((wndw-1)/2) ) ){x+=1;xVals.push(x+shift)};
                    break;
            };
        break;
        case "savitzky-golay":
            switch (type) {
                case "protein":
                    smoothArray = savitzky_golay(seqArray,wndw,1);
                    var xVals = [Math.round((wndw-1)/2)+shift], x=Math.round((wndw-1)/2);
                    while(x<(smoothArray.length+Math.round((wndw-1)/2) ) ){x+=1;xVals.push(x+shift)};
                    break;

                case "dna":
                case "rna":
                    smoothArray = savitzky_golay(seqArray,wndw,3)
                    var xVals = [Math.round((wndw-1)/2)+shift], x=Math.round((wndw-1)/2);
                    while(x<(smoothArray.length+Math.round((wndw-1)/2) ) ){x+=1;xVals.push(x+shift)};
                    break;
            };

    };

    var data = []
    for (var i = 0; i < smoothArray.length; i++) {
        if (smoothWeights[i] >= min_existing){
            data.push({"x":xVals[i], "y":smoothArray[i], "w":smoothWeights[i]});
        }
        else {
            data.push({"x":xVals[i], "y":NaN, "w":smoothWeights[i]});
        }

    };

    return data;
};

function drawingProfiles(data, id = ".canvas") {

    data.seq.forEach(function(d) {
        if (d.active && d.visible) {
            d.profile = profile_calculation(d.sequence, d.window, d.scale, d.type, d.shift, d.smoothing_method, d.min_existing);
        }
    });


    input = {"seq": []};
    input.title = data.title;
    input.seq = data.seq.filter(function (n) {
        return n.visible && n.active;
    });
    d3.selectAll("svg > *").remove();
    d3.selectAll("#SVGholder > div").remove();
    var tooltip = d3.select('#SVGholder').append('div')
        .attr('class', 'SVGhidden SVGtooltip');
    var svg = d3.select(id),
        margin = {top: 30, right: 70, bottom: 40, left: 70},
        width = +svg.attr("width") - margin.left - margin.right,
        height = + svg.attr("height") - margin.top - margin.bottom,
        g = svg.append("g").attr("transform", "translate(" + margin.left +"," + margin.top +")")
               .attr("class", "canvas_g");

    var mask = svg.append("defs").append("clipPath")
                .attr("id", "clip")
                .append("rect")
                .attr("width", width)
                .attr("height", height)


    var x1_min = d3.min(input.seq.map(function(d) {return d["shift"]}));
    var x1_max = d3.max(input.seq.map(function(d) {return d["profile"].length + d["shift"] + d["window"]}));

    var x1 = d3.scaleLinear()
               .domain([x1_min, x1_max])
               .rangeRound([0,width]);

    var xAxis = d3.axisBottom(x1);
    var gX1 = g.append('g')
               .classed("axisBlack", true)
               .classed("x1", true)
               .classed("axis", true)
               .attr("transform", "translate(0," + height + ")")
               .call(xAxis);



    // setting xlabel
    if (input.seq.map(function(d) {return d["type"]}).includes("rna") || input.seq.map(function(d) {return d["type"]}).includes("dna")) {
        if (input.seq.map(function(d) {return d["type"]}).includes("protein")) {
            g.append("text")
                .attr("y", height + margin.bottom)
                .attr("x", width* (2/5))
                .text("Residue / Codon")
        }
        else {
            g.append("text")
                .attr("y", height + margin.bottom)
                .attr("x", width* (2/5))
                .text("Codon")
        }
    }
    else {
        if (input.seq.map(function(d) {return d["type"]}).includes("protein")) {
            g.append("text")
                .attr("y", height + margin.bottom)
                .attr("x", width* (2/5))
                .text("Residue")
        }
        else {
            g.append("text")
                .attr("y", height + margin.bottom)
                .attr("x", width* (2/5))
                .text("")
        }
    }

    if (input.seq.map(function(d) {return d["type"]}).includes("protein")) {

        if (data.relative) {
            rescale(data.seq, 'protein');

            var y1_min = findMinGlobal(data.seq);
            var y1_max = findMaxGlobal(data.seq);
            if (Math.abs(y1_min) >= Math.abs(y1_max)) {
                y1_max = y1_min*(-1)
            } else {
                y1_min = y1_max*(-1)
            };


        } else {
            var y1_min = findMin(data.seq, "protein");
            var y1_max = findMax(data.seq, "protein");
        };

        if (data.y1_reverse) {
            padding1 = (y1_min - y1_max)*0.1;
            var gap_line_position_protein = y1_min + padding1*0.8
            var y1 = d3.scaleLinear()
                        .domain([y1_max - padding1, y1_min + padding1])
                        .rangeRound([height, 0]);
        } else {
            var padding1 = (y1_max - y1_min)*0.1;
            var gap_line_position_protein = y1_max + padding1*0.8
            var y1 = d3.scaleLinear()
                    .domain([y1_min - padding1, y1_max + padding1])
                    .rangeRound([height, 0]);
        }

        var yAxis1 = d3.axisLeft(y1);

        var protLine = d3.line()
            .x(function(d) {return x1(d.x); })
            .y(function(d) {return y1(d.y); })
            .defined(function(d) {return !isNaN(d.y);});

        var protGapLine  = d3.line()
            .x(function(d) {return x1(d.x); })
            .y(function(d) {return y1(gap_line_position_protein); })
            .defined(function(d) {return (d.w != 1 && !isNaN(d.y));});



        var gY1 = g.append("g")
                   .classed("axisBlack", true)
                   .classed("yAxis", true)
                   .classed("axis", true)
                   .call(yAxis1);
        if (data.relative) {
            g.append("text")
                .attr("transform", "rotate(-90)")
                .attr("y", 0 - margin.left*(2/3))
                .attr("x", 0 - height *(3/5))
                .text("Protein [sd]")
        } else {
            g.append("text")
                .attr("transform", "rotate(-90)")
                .attr("y", 0 - margin.left*(2/3))
                .attr("x", 0 - height *(3/5))
                .text("Protein")
        };
    }
    if (input.seq.map(function(d) {return d["type"]}).includes("rna") || input.seq.map(function(d) {return d["type"]}).includes("dna")) {

        if (data.relative) {
            rescale(data.seq, 'rna');
            rescale(data.seq, 'dna');

            var y2_min = findMinGlobal(input.seq);
            var y2_max = findMaxGlobal(input.seq);
            if (Math.abs(y2_min) >= Math.abs(y2_max)) {
                y2_max = y2_min*(-1)
            } else {
                y2_min = y2_max*(-1)
            };
        } else {
            var y2_min = d3.min([findMin(input.seq, "rna"), findMin(input.seq, "dna")]);
            var y2_max = d3.max([findMax(input.seq, "rna"), findMax(input.seq, "dna")]);
        }

        if (data.y2_reverse) {
            padding2 = (y2_min - y2_max)*0.1;
            var gap_line_position_rna = y2_min + padding2*0.8
            var y2 = d3.scaleLinear()
                        .domain([y2_max - padding2, y2_min + padding2])
                        .rangeRound([height, 0]);
        } else {
            var padding2 = (y2_max - y2_min)*0.1;
            var gap_line_position_rna = y2_max + padding2*0.8

            var y2 = d3.scaleLinear()
                        .domain([y2_min - padding2, y2_max + padding2])
                        .rangeRound([height, 0]);
        }

        var yAxis2 = d3.axisRight(y2);

        var rnaLine  = d3.line()
            .x(function(d) {return x1(d.x); })
            .y(function(d) {return y2(d.y); })
            .defined(function(d) {return (!isNaN(d.y));})
            ;
        var rnaGapLine  = d3.line()
            .x(function(d) {return x1(d.x); })
            .y(function(d) {return y2(gap_line_position_rna); })
            .defined(function(d) {return (d.w != 1 && !isNaN(d.y));})
            ;


        var gY2 = g.append("g")
                   .classed("axisBlack", true)
                   .classed("yAxis", true)
                   .classed("axis", true)
                   .attr("transform", "translate(" + width + ",0)")
                   .call(yAxis2);
        if (data.relative) {
            g.append("text")
                .attr("y", -width - margin.right * (2/3))
                .attr("x", height*(2/5))
                .attr("transform", "rotate(90)")
                .text("RNA / DNA [sd]")
        } else {
            g.append("text")
                .attr("y", -width - margin.right * (2/3))
                .attr("x", height*(2/5))
                .attr("transform", "rotate(90)")
                .text("RNA / DNA")
        }
    }

    var zoom = d3.zoom()
                 .scaleExtent([0.25, 8])
                 .translateExtent([[-width, -height], [2 * width, 2 *height]])
                 .on("zoom", zoomed);

    var zoomRect = svg.call(zoom);

    function zoomed() {
        var xz = d3.event.transform.rescaleX(x1);

        gX1.call(xAxis.scale(xz));

        if (input.seq.map(function(d) {return d["type"]}).includes("rna") || (input.seq.map(function(d) {return d["type"]}).includes("dna") )) {
            rnaLine
                .x(function(d) {return xz(d.x); })
            rnaGapLine
                .x(function(d) {return xz(d.x); })
        }
        if (input.seq.map(function(d) {return d["type"]}).includes("protein")) {
            protLine
                .x(function(d) {return xz(d.x); })
            protGapLine
                .x(function(d) {return xz(d.x); })
        }
        g.selectAll(".dot")
                .attr('cx', function(d) {return xz(d.x)})

        d3.selectAll(".protein").attr("d", protLine)
        d3.selectAll(".rna").attr("d", rnaLine)
        d3.selectAll(".rnagaps").attr("d", rnaGapLine)
        d3.selectAll(".proteingaps").attr("d", protGapLine)

    }

    function drawLines(item, index) {
        g = d3.select(".canvas_g")
        if (item.type == 'rna' || item.type == "dna") {

            var currentLine = g.append("path")
             .classed('rna', true)
             .datum(item.profile)
             .attr("clip-path", "url(#clip)")
             .attr("fill", "none")
             .attr("stroke", item.color)
             .attr("stroke-linejoin", "round")
             .attr("stroke-linecap", "round")
             .attr("stroke-width", item.thickness)
             .attr("d", rnaLine)
             .on('mousemove', function(d) {
                d3.select(".SVGtooltip").classed('SVGhidden', false)
                   .attr('style', 'left:' + (d3.event.clientX + 20) + 'px; top:' + (d3.event.clientY - 20) + 'px; text-align: left;')
                   .html("<b>"+removeScaleFromName(input.seq[index].name)+"</b><br>" + trimString(scaleDescriptors[item.type][input.seq[index].scale]["name"], 80));
            })
            .on('mouseout', function() {
               d3.select(".SVGtooltip").classed('SVGhidden', true);
            })
            if (item.linestyle == 'dashed1') {
                currentLine.attr("stroke-dasharray", (item.thickness + " " + (Number(item.thickness)*2)));
            } else if (item.linestyle == 'dashed2') {
                currentLine.attr("stroke-dasharray", (item.thickness + " " + (Number(item.thickness)*3)));
            } else if (item.linestyle == 'dotted') {
                currentLine.attr("stroke-dasharray", ("1 " + (Number(item.thickness)*3)));
            }

            g.append("path")
             .classed('rnagaps', true)
             .datum(item.profile)
             .attr("clip-path", "url(#clip)")
             .attr("fill", "none")
             .attr("stroke", item.color)
             .attr("stroke-linejoin", "round")
             .attr("stroke-linecap", "round")
             .attr("stroke-width", item.thickness)
             .attr("d", rnaGapLine)
             .attr("stroke-dasharray", (item.thickness + " " + (Number(item.thickness)*2)))
             .on('mousemove', function(d) {
                d3.select(".SVGtooltip").classed('SVGhidden', false)
                   .attr('style', 'left:' + (d3.event.clientX + 20) + 'px; top:' + (d3.event.clientY - 20) + 'px; text-align: left;')
                   .html("<b>"+removeScaleFromName(input.seq[index].name)+"</b><br>some missing values in averaging window");
            })
            .on('mouseout', function() {
               d3.select(".SVGtooltip").classed('SVGhidden', true);
            })




        } else if (item.type == 'protein') {
            var currentLine = g.append("path")
             .classed('protein', true)
             .datum(item.profile)
             .attr("clip-path", "url(#clip)")
             .attr("fill", "none")
             .attr("stroke", item.color)
             .attr("stroke-linejoin", "round")
             .attr("stroke-linecap", "round")
             .attr("stroke-width", item.thickness)
             .attr("d", protLine)
             .on('mousemove', function(d) {
                 d3.select(".SVGtooltip").classed('SVGhidden', false)
                    .attr('style', 'left:' + (d3.event.clientX + 20) + 'px; top:' + (d3.event.clientY - 20) + 'px; text-align: left;')
                    .html("<b>"+removeScaleFromName(input.seq[index].name)+"</b><br>" + trimString(scaleDescriptors["protein"][input.seq[index].scale]["name"], 80));
             })
             .on('mouseout', function() {
                d3.select(".SVGtooltip").classed('SVGhidden', true);
             })
            if (item.linestyle == 'dashed1') {
                currentLine.attr("stroke-dasharray", (item.thickness + " " + (Number(item.thickness)*2)));
            } else if (item.linestyle == 'dashed2') {
                currentLine.attr("stroke-dasharray", (item.thickness + " " + (Number(item.thickness)*3)));
            } else if (item.linestyle == 'dotted') {
                currentLine.attr("stroke-dasharray", ("1 " + (Number(item.thickness)*3)));
            }

            g.append("path")
             .classed('proteingaps', true)
             .datum(item.profile)
             .attr("clip-path", "url(#clip)")
             .attr("fill", "none")
             .attr("stroke", item.color)
             .attr("stroke-linejoin", "round")
             .attr("stroke-linecap", "round")
             .attr("stroke-width", item.thickness)
             .attr("d", protGapLine)
             .attr("stroke-dasharray", (item.thickness + " " + (Number(item.thickness)*2)))
             .on('mousemove', function(d) {
                d3.select(".SVGtooltip").classed('SVGhidden', false)
                   .attr('style', 'left:' + (d3.event.clientX + 20) + 'px; top:' + (d3.event.clientY - 20) + 'px; text-align: left;')
                   .html("<b>"+removeScaleFromName(input.seq[index].name)+"</b><br>some missing values in averaging window");
            })
            .on('mouseout', function() {
               d3.select(".SVGtooltip").classed('SVGhidden', true);
            })
        }

    };

    input.seq.forEach(drawLines);

    updateCorResults(corList);
    if (document.getElementsByClassName('.correlation-data')) {
        addTable('.correlation-data', data);
    }
    return "complete";
};



var closeModal = function() {
    var alphabet = {
        "rna" : /^[AGTUC-]*$/gi,
        "dna": /^[AGTUC-]*$/gi,
        "protein": /^[ACDEFGHIKLMNPQRSTVWY-]*$/gi
    };

    var that = $(this)
    if ((data.seq[this.id].scale != "no scale selected") && (data.seq[this.id].scale) && (data.seq[this.id].sequence != "")) {

        if (alphabet[data.seq[this.id].type].test(data.seq[this.id].sequence)) {
             $('#setup' + this.id).modal('hide');
        }
        else {
            $("#setup" + this.id + " .alert_box2").empty();
            $("#setup" + this.id + " .alert_box2").append("<div class='alert alert-danger alert-dismissable'><button type='button' class='close' data-dismiss='alert' aria-hidden='true'>&times;</button><strong>Warning</strong> Invalid characters in " + typeNames[data.seq[this.id].type]+" sequence - if retained, these will be treated as gaps</div>");
            $("#setup" + this.id + " .alert_box").empty();
            $("#setup" + this.id + " .alert_box").append("<div class='alert alert-danger alert-dismissable'><button type='button' class='close' data-dismiss='alert' aria-hidden='true'>&times;</button><strong>Warning</strong> Invalid characters in " + typeNames[data.seq[this.id].type]+" sequence - if retained, these will be treated as gaps</div>");

        }
    }
    else {
        that.tooltip('show');
        setTimeout(function(){
            that.tooltip('hide')
        }, 2000);
    }

}




var clickScaleTab = function() {
    var that = $(this)
    var i = this.id.split("_")[1]

    if ($("#" + i + ".btn-type").html().trim() == "Type") {
        that.tooltip('show');
        setTimeout(function(){
            that.tooltip('hide')
        }, 2000);
    }

}

const Item = ({number}) =>`
<div class="sequence">
    <div class="title_bar">
    <button class="color-picker-btn" id="colorpicker${number}"></button>
        <div class="name-value">${number+1}:</div>
        <input class="name" type="text" id="${number}" value="Seq${number+1}" autocomplete="off" autocorrect="off" autocapitalize="off" spellcheck="false"></input>
        <div class="close-button" data-toggle="tooltip" data-placement="right" title="Remove data" id="${number}"><i class="material-icons md-title">close</i></div>
        <div class="button" id="${number}"><i class="material-icons md-24">expand_less</i></div>
        <div class="visible" id="${number}"><i class="material-icons md-title">visibility</i></div>
        <button type="button" class="settings" data-toggle="modal" data-target="#setup${number}"><i class="material-icons md-title-small">build</i></button>

        <div class="modal" id="setup${number}">
            <div class="modal-dialog modal-lg">
                <div class="modal-content">
                    <div class="modal-header">
                        <h4 class="modal-title">Options</h4>
                        <button type="button" class="close" id="${number}" data-dismiss="modal">&times;</button>
                    </div>
                    <div class="modal-body">
                        <ul class="nav nav-tabs" role="tablist">
                            <li class="nav-item">
                                <a class="nav-link active" data-toggle="tab" href="#sequence${number}">Sequence</a>
                            </li>
                            <li class="nav-item">
                                <a class="nav-link disabled" data-toggle="tab" data-toggle="tooltip" data-placement="right" id="scaleTab_${number}" title="Please select molecule type before continuing"  data-trigger="manual"  onclick="clickScaleTab.call(this)" href="#scale${number}">Scale</a>
                            </li>
                            <li class="nav-item">
                                <a class="nav-link" data-toggle="tab" href="#visuals${number}">Visuals</a>
                            </li>
                        </ul>

                        <!-- Tab panes -->
                        <div class="tab-content">

                            <div id="sequence${number}" class="container tab-pane active"><br>
                                <textarea class="form-inline" id="${number}" type="text" size="300" style="display: block; width: 100%; height: 8em; resize: none; font-family: 'Source Code Pro', monospace;"></textarea>

                                <div class="modal-footer">
                                    <div style="display:flex;align-items:center;justify-content: space-between; width: 100%;">
                                        <div style="display:flex;align-items:center;justify-content: space-between;">
                                            <div class="btn-group type-menu" style="float: left;">
                                                <button class="btn btn-grey btn-sm btn-type" id="${number}" type="button">
                                                    Type
                                                </button>
                                                <button type="button" class="btn btn-sm btn-grey dropdown-toggle dropdown-toggle-split" data-toggle="dropdown" aria-haspopup="true" aria-expanded="false">
                                                    <span class="sr-only">Toggle Dropdown</span>
                                                </button>
                                                <div class="dropdown-menu type-selector" id="${number}">
                                                    <a class="dropdown-item" id="protein_dropdown_item${number}" >Protein</a>
                                                    <a class="dropdown-item" id="rna_dropdown_item${number}">RNA</a>
                                                    <a class="dropdown-item" id="dna_dropdown_item${number}">DNA</a>

                                                </div>
                                            </div>
                                            <div class="btn-group database-menu" id="${number}" style="float: left;">
                                                <button class="btn btn-grey btn-sm btn-label" id="${number}" type="button">
                                                    Database
                                                </button>
                                                <button type="button" class="btn btn-sm btn-grey dropdown-toggle dropdown-toggle-split btn-flat" data-toggle="dropdown" aria-haspopup="true" aria-expanded="false">
                                                    <span class="sr-only">Toggle Dropdown</span>
                                                </button>
                                                <div class="dropdown-menu database-selector" id="${number}">
                                                    <a class="dropdown-item" id="uniprot_dropdown_item${number}">Uniprot</a>
                                                    <a class="dropdown-item" id="ena_dropdown_item${number}">ENA</a>
                                                </div>
                                            </div>
                                            <input class="sequence-id-input" id="${number}" style="margin-right: 7px;" type="text" placeholder="Enter sequence ID"/>
                                            <button type="button" class="btn btn-sm btn-blue btn-submit btn-flat" id="${number}" style="float: left;">Go!</button>
                                        </div>

                                    </div>
                                    <button type="button" class="btn btn-sm btn-red btn-clr btn-flat" id="${number}">Clear</button>
                                    <button type="button" class="btn btn-sm btn-blue btn-flat" id="${number}" data-toggle="tooltip" data-placement="right" title="Please select a scale and a sequence before continuing" data-trigger="manual" onclick="closeModal.call(this)">Done</button>

                                </div>
                                <div class="modal-footer">
                                    <div class="alert_box" style="display: inline-block; width: 100%; min-width: 350px; font-size: 15px"></div>
                                </div>
                            </div>

                            <div id="scale${number}" class="container tab-pane fade"><br>
                                <div style="width: 35%; height: 90%; float: left;">
                                    <select class="selectpicker" id="${number}" multiple>
                                        <option value="p_alpha">Alpha-propensity</option>
                                        <option value="p_beta">Beta-propensity</option>
                                        <option value="p_charge">Charge</option>
                                        <option value="p_comp">Composition</option>
                                        <option value="p_hydro">Hydrophobicity</option>
                                        <option value="p_pchem">general Physico-Chemical</option>
                                        <option value="p_rnaAf">RNA-Affinity</option>
                                        <option value="p_other">Other</option>
                                    </select>
                                    <select class="scaleList" id="${number}" size="10" style="width: 220px; margin-top: 20px;">

                                    </select>

                                </div>
                                <div style="width: 65%; height: 90%; min-height: 270px; float: left;">
                                    <div class="scaleDescription" id="${number}" style="font-size: 12px;">

                                    </div>
                                </div><br>
                                <div class="modal-footer">
                                    <button type="button" class="btn btn-sm btn-blue btn-flat" id="${number}" data-toggle="tooltip" data-placement="right" title="Please select a scale and a sequence before continuing" data-trigger="manual" onclick="closeModal.call(this)">Done</button>
                                </div>
                                <div class="alert_box2" style="display: inline-block; width: 100%; min-width: 350px; font-size: 15px"></div>
                            </div>

                            <div id="visuals${number}" class="container tab-pane fade"><br>
                                <div class="box-card">
                                    <div id="${number}" class="thickness-text" style="width: 14em;">Thickness: 2.0</div>
                                    <div class="thickness-slider-container">
                                        <div class="thickness-slider-inner">
                                            <button id="${number}" class="btn btn-sq-xs btn-less-thick">
                                                <i class="material-icons" md-24>remove</i>
                                            </button>
                                            <input class="thickness-slider" data-slider-id='thickness-slider' data-slider-tooltip="hide" type="text" data-slider-min="0.1" data-slider-max="10" data-slider-step="0.1" data-slider-value="2" style="width: 20em;" id='${number}'/>
                                            <button id="${number}" class="btn btn-sq-xs btn-more-thick">
                                                <i class="material-icons" md-24>add</i>
                                            </button>
                                        </div>
                                    </div>
                                </div>
                                <div class="modal-footer">
                                    <div style="display:flex;align-items:center; width: 100%;">
                                        <span class="sliderSwitch-label" style="margin-left: 0px; margin-top: 0px; width: 14em;">Linetype:</span>
                                        <select class="form-control linestyleSelector" id="${number}" style="width:40%;">
                                            <option value="solid" selected>Solid</option>
                                            <option value="dashed1">Dashed 1</option>
                                            <option value="dashed2">Dashed 2</option>
                                            <option value="dotted">Dotted</option>
                                        </select>
                                    </div>
                                </div>
                                <div class="modal-footer">
                                    <div style="display:flex;align-items:center; width: 100%;">
                                        <span class="min_existing-label" style="margin-left: 0px; width: 14em;">Minimum percentage of existing values in window:</span>
                                        <div class="input-group" style="width:18%">
                                            <input class="form-control min_existing-input" style="border-top-right-radius: 0.25rem;border-bottom-right-radius: 0.25rem" id="${number}" type="number" min="0" max="100" step="1" value="0"/>
                                            <div class="input-group-append">
                                                <span class="input-group-text" style="background-color:#fff; border:none; padding-left:0.45rem" >%</span>
                                            </div>
                                        </div>
                                    </div>
                                </div>
                                <div class="modal-footer">
                                    <button type="button" class="btn btn-sm btn-blue btn-flat" id="${number}" data-toggle="tooltip" data-placement="right" title="Please select a scale and a sequence before continuing" data-trigger="manual" onclick="closeModal.call(this)">Done</button>
                                </div>
                                <div class="alert_box2" style="display: inline-block; width: 100%; min-width: 350px; font-size: 15px"></div>
                            </div>
                        </div>
                    </div>
                </div>
            </div>
        </div>



    </div>
    <div class="box">
        <div id="${number}" class="firstRow-text">Window: 21</div>
        <div class="firstRow-slider">
            <div class="firstRow-slider-inner">
                <button id="${number}" class="btn btn-sq-xs btn-less-wndw">
                    <i class="material-icons" md-24>remove</i>
                </button>
                <input class="prot-window" data-slider-id='prot-window' data-slider-tooltip="hide" type="text" data-slider-min="1" data-slider-max="63" data-slider-step="2" data-slider-value="21" style="width: 20em;" id='${number}'/>
                <button id="${number}" class="btn btn-sq-xs btn-more-wndw">
                    <i class="material-icons" md-24>add</i>
                </button>
            </div>
        </div>

        <div id="${number}" class="thirdRow-text">Shift: 0</div>
        <div class="thirdRow-slider">
            <div class="thirdRow-slider-inner">
                <button id="${number}" class="btn btn-sq-xs btn-less">
                    <i class="material-icons md-24">remove</i>
                </button>
                <input class="shift" id="${number}" data-provide="slider" data-slider-id='shift-dist' data-slider-tooltip="hide" type="text" data-slider-min="0" data-slider-max="1" data-slider-step="1" data-slider-value="0" style="width: 20em;" />
                <button id="${number}" class="btn btn-sq-xs btn-more">
                    <i class="material-icons" md-24>add</i>
                </button>
            </div>
        </div>
    </div>
</div>
`;

const Legend = ({number}) =>`
    <div class="color-display">
    </div>
    <div class="legend-text">
    </div>
`;

data_prototype = {
                  "type" : "protein",
                  "scale" : "no scale selected",
                  "sequence" : "",
                  "window" : 21,
                  "shift" : 0,
                  "smoothing_method" : "boxcar",
                  "name" : "",
                  "organism" : "",
                  "thickness" : 2,
                  "min_existing": 0,
                  "linestyle" : "solid",
                  "color" : "",
                  "active" : true,
                  "visible" : true
                }

function addInterface() {
    $("#addInterfaceButton").tooltip('hide');
    $("#addInterfaceButton").text("Add additional sequence")
    var div = document.createElement('div');

    num = Object.keys(data.seq).length

    div.setAttribute('class', 'interfacebox');
    div.setAttribute('id', num);

    div.innerHTML = Item({number: num});

    data.seq[num] = $.extend(true, {}, data_prototype);
    data.seq[num].sequence = "";

    data.seq[num].color = getBiasedColorHSL(num, .60, .60, .2);

    data.seq[num].name = "Seq" + (num+1);
    data.seq[num].name = addScaleToName(num,data.seq[num].name)



    document.getElementById('interface_container').appendChild(div);

    new jscolor(document.getElementById('colorpicker' + num),
                {valueElement:'none', onFineChange:'updateColor(this, ' + num + ')', value:data.seq[num].color.slice(1)})

    $('#'+num+'.name').val(data.seq[num].name);

    window_p = $('.prot-window').each(function(element, value) {
        $(this).slider()
        .on("slide", function() {
            if (data.seq[this.id].type == "protein") {
                data.seq[this.id].window = Number(this.value);
                $("#"+this.id+".firstRow-text").text("Window: " + this.value);
            } else {
                data.seq[this.id].window = Number(this.value);
                $("#"+this.id+".firstRow-text").text("Window: " + Number(this.value));
            }
            drawingProfiles(data);
        })
        .on("slideStop", function() {
            if (data.seq[this.id].type == "protein") {
                data.seq[this.id].window = Number(this.value);
                $("#"+this.id+".firstRow-text").text("Window: " + this.value);
            } else {
                data.seq[this.id].window = Number(this.value);
                $("#"+this.id+".firstRow-text").text("Window: " + Number(this.value));
            }
            drawingProfiles(data);
        })

        .data('slider');
    });

    shift_p = $('.shift').each(function(element, value) {
        $(this).slider()
        .on("slide", function() {
            data.seq[this.id].shift = Number(this.value);
            $("#"+this.id+".thirdRow-text").text("Shift: " + this.value);
            drawingProfiles(data);
        })
        .on("slideStop", function() {
            data.seq[this.id].shift = Number(this.value);
            $("#"+this.id+".thirdRow-text").text("Shift: " + this.value);
            drawingProfiles(data);
        })
        .data('slider');
    });

    thickness_p = $('.thickness-slider').each(function(element, value) {
        $(this).slider()
        .on("slide", function() {
            data.seq[this.id].thickness = Number(this.value);
            $("#"+this.id+".thickness-text").text("Thickness: " + data.seq[this.id].thickness);
            drawingProfiles(data);
        })
        .on("slideStop", function() {
            data.seq[this.id].thickness = Number(this.value);
            $("#"+this.id+".thickness-text").text("Thickness: " + data.seq[this.id].thickness);
            drawingProfiles(data);
        })
        .data('slider');
    });

    $('.form-inline').each(function(element, value) {
        $(this).val(data.seq[this.id].sequence);
        });

    $('.form-inline').each(function(element, value) {
        $(this).bind('input propertychange', function() {
            data.seq[this.id].sequence = parseFASTA(this.value);
            drawingProfiles(data);
            updateShifts(data);
        });
    });
    $('.name').each(function(element, value) {
        $(this).on("change", function() {
            data.seq[this.id].name = this.value;
            updateCorrList(data);
            updateLegend(this.id);
        })
    });
    $('#'+num+' .linestyleSelector').on('change', function() {
        data.seq[this.id].linestyle = this.value;
        drawingProfiles(data);
    });
    $('#'+num+' .min_existing-input').on('change', function() {
        data.seq[this.id].min_existing = Number(this.value)/100;
        drawingProfiles(data);
    });
    $("div#"+num+".button").click(function(){
        if($(this).html() == '<i class="material-icons md-24">expand_less</i>'){
            $(this).html('<i class="material-icons md-24">expand_more</i>');
        }
        else{
            $(this).html('<i class="material-icons md-24">expand_less</i>');
        }
        $(this).parent().parent().find('.box').slideToggle(50, 'linear');
    });
    $("div#"+num+".visible").click(function(){
        if($(this).html() == '<i class="material-icons md-title">visibility</i>'){
            $(this).html('<i class="material-icons md-title">visibility_off</i>');
            data.seq[this.id].visible = false;
            drawingProfiles(data);
        }
        else{
            $(this).html('<i class="material-icons md-title">visibility</i>');
            data.seq[this.id].visible = true;
            drawingProfiles(data);
        }

    });

    $(".close-button").tooltip();

    $("#setup"+num).modal('show');

    $("#setup"+num).on('hidden.bs.modal', function() {
        $("#setup"+num + " .alert_box").empty();
    });

    $('#'+num+'.selectpicker').selectpicker('refresh');

    $('#'+num+'.selectpicker').on('change', function(){
        createScaleList(this.id, $(this).val(), data.seq[this.id].type);
    });

    $('.scaleList').on('change', function(){
        renderScaleDetails(this.id, $(this).val(), data.seq[this.id].type);
        data.seq[this.id].scale = $(this).val();
        data.seq[this.id].name = addScaleToName(this.id, data.seq[this.id].name);
        $('#'+this.id+'.name').val(data.seq[this.id].name);
        drawingProfiles(data);
    });

    addLegend(num);
}

function addInterfaceBasedOnData(data) {
    $("#addInterfaceButton").tooltip('hide');
    $("#addInterfaceButton").text("Add additional sequence")

    $("#interface_container").empty()

    for (num in data.seq) {
        num = Number(num);
        var div = document.createElement('div');
        div.setAttribute('class', 'interfacebox');
        div.setAttribute('id', num);

        div.innerHTML = Item({number: num});

        document.getElementById('interface_container').appendChild(div);

        new jscolor(document.getElementById('colorpicker' + num),
                    {valueElement:'none', onFineChange:'updateColor(this, ' + num + ')', value:data.seq[num].color.slice(1)});

        $('#'+num+'.name').val(data.seq[num].name);
        $('#'+num+'.name').val(addScaleToName(num,data.seq[num].name));
        $('#'+num+'.scale_label').html(data.seq[num].scale);

        $("#scaleTab_"+num).removeClass("disabled");
        $("#" + num + ".btn-type").html(trimString(typeNames[data.seq[num].type], 9));
        $("#" + num + ".btn-type").val(trimString(typeNames[data.seq[num].type], 9));

        if (data.seq[num].type == "protein") {
            var db = "Uniprot";
        }
        else {
            var db = "ENA";
        }

        $("#"+num+".database-menu").find('.btn-label').html(db);
        $("#"+num+".database-menu").find('.btn-label').val(db);


        addLegend(num);

        $("div#"+num+".button").click(function(){
            if($(this).html() == '<i class="material-icons md-24">expand_less</i>'){
                $(this).html('<i class="material-icons md-24">expand_more</i>');
            }
            else{
                $(this).html('<i class="material-icons md-24">expand_less</i>');
            }
            $(this).parent().parent().find('.box').slideToggle(50, 'linear');
        });
        $("div#"+num+".visible").click(function(){
            if($(this).html() == '<i class="material-icons md-title">visibility</i>'){
                $(this).html('<i class="material-icons md-title">visibility_off</i>');
                data.seq[this.id].visible = false;
                drawingProfiles(data);
            }
            else{
                $(this).html('<i class="material-icons md-title">visibility</i>');
                data.seq[this.id].visible = true;
                drawingProfiles(data);
            }

        });
        if (data.seq[num].type == "protein") {
            $("#"+num+".selectpicker").html(`
                                    <option value="p_alpha">Alpha-propensity</option>
                                    <option value="p_beta">Beta-propensity</option>
                                    <option value="p_charge">Charge</option>
                                    <option value="p_comp">Composition</option>
                                    <option value="p_hydro">Hydrophobicity</option>
                                    <option value="p_pchem">general Physico-Chemical</option>
                                    <option value="p_rnaAf">RNA-Affinity</option>
                                    <option value="p_other">Other</option>
                                    `)
            $("#"+num+".selectpicker").selectpicker('refresh')
        }
        else if (data.seq[num].type == "rna") {
            $("#"+num+".selectpicker").html(`
                                    <option value="r_comp">Composition</option>
                                    `)
            $("#"+num+".selectpicker").selectpicker('refresh')
        }
        else if (data.seq[num].type == "dna") {
            $("#"+num+".selectpicker").html(`
                                    <option value="d_comp">Composition</option>
                                    `)
            $("#"+num+".selectpicker").selectpicker('refresh')
        }
        $('#'+num+'.selectpicker').on('change', function(){
            createScaleList(this.id, $(this).val(), data.seq[this.id].type);
        });
        $('#'+num+'.selectpicker').selectpicker('refresh');
        $('#'+num+' .linestyleSelector').on('change', function() {
            data.seq[this.id].linestyle = this.value;
            drawingProfiles(data);
        });
        $('#'+num+' .min_existing-input').on('change', function() {
            data.seq[this.id].min_existing = Number(this.value)/100;
            drawingProfiles(data);
        });
        renderScaleDetails(num, data.seq[num].scale, data.seq[num].type);
    }

    window_p = $('.prot-window').each(function(element, value) {
        $(this).slider()
        .on("slide", function() {
            if (data.seq[this.id].type == "protein") {
                data.seq[this.id].window = Number(this.value);
                $("#"+this.id+".firstRow-text").text("Window: " + this.value);
            } else {
                data.seq[this.id].window = Number(this.value);
                $("#"+this.id+".firstRow-text").text("Window: " + Number(this.value));
            }
            drawingProfiles(data);
        })
        .on("slideStop", function() {
            if (data.seq[this.id].type == "protein") {
                data.seq[this.id].window = Number(this.value);
                $("#"+this.id+".firstRow-text").text("Window: " + this.value);
            } else {
                data.seq[this.id].window = Number(this.value);
                $("#"+this.id+".firstRow-text").text("Window: " + Number(this.value));
            }
            drawingProfiles(data);
        })
        .data('slider');
    });

    shift_p = $('.shift').each(function(element, value) {
        $(this).slider()
        .on("slide", function() {
            data.seq[this.id].shift = Number(this.value);
            $("#"+this.id+".thirdRow-text").text("Shift: " + this.value);
            drawingProfiles(data);
        })
        .on("slideStop", function() {
            data.seq[this.id].shift = Number(this.value);
            $("#"+this.id+".thirdRow-text").text("Shift: " + this.value);
            drawingProfiles(data);
        })
        .data('slider');
    });

    thickness_p = $('.thickness-slider').each(function(element, value) {
        $(this).slider()
        .on("slide", function() {
            data.seq[this.id].thickness = Number(this.value);
            $("#"+this.id+".thickness-text").text("Thickness: " + data.seq[this.id].thickness);
            drawingProfiles(data);
        })
        .on("slideStop", function() {
            data.seq[this.id].thickness = Number(this.value);
            $("#"+this.id+".thickness-text").text("Thickness: " + data.seq[this.id].thickness);
            drawingProfiles(data);
        })
        .data('slider');
    });

    $('.form-inline').each(function(element, value) {
        $(this).val(data.seq[this.id].sequence);
        });

    $('.form-inline').each(function(element, value) {
        $(this).bind('input propertychange', function() {
            data.seq[this.id].sequence = parseFASTA(this.value);
            drawingProfiles(data);
            updateShifts(data);
        });
    });
    $('.name').each(function(element, value) {
        $(this).on("change", function() {
            data.seq[this.id].name = this.value;
            updateCorrList(data);
            updateLegend(this.id);
        })
    });

    $('.scaleList').on('change', function(){
        renderScaleDetails(this.id, $(this).val(), data.seq[this.id].type);
        data.seq[this.id].scale = $(this).val();
        data.seq[this.id].name = addScaleToName(this.id, data.seq[this.id].name);
        $('#'+this.id+'.name').val(data.seq[this.id].name);
        drawingProfiles(data);
        updateShifts(data)
    });


    $(".close-button").tooltip();

    drawingProfiles(data);
    updateShifts(data);
    updateWindows(data);

}

function updateCorrList(data) {
    $(".corr-selector").empty();
    activeData = getActiveData(data);
    for (var i = data.seq.length - 1; i >= 0; i--) {
        $("#"+i+".corr-selector").append(getCorOptions(activeData, i));
    }
}

function updateShifts(data) {
    $(".shift").slider('setAttribute', 'max', findMaxLength(data.seq));
    $(".shift").slider('refresh');
    for (var i = data.seq.length - 1; i >= 0; i--) {
        if (data.seq[i].active) {
            $("#"+i+".shift").data('slider').setValue(data.seq[i].shift);
            $("#"+i+".thirdRow-text").text("Shift: " + data.seq[i].shift);
        }
    }
}

function updateWindows(data) {
    for (var i = data.seq.length - 1; i >= 0; i--) {
        if (data.seq[i].active) {
            $("#"+i+".prot-window").data('slider').setValue(data.seq[i].window);
            $("#"+i+".firstRow-text").text("Window: " + data.seq[i].window);

        }
    }
}

function resetScale(id) {
    $("#"+id+".scaleDescription").html("");
    createScaleList(id, [], data.seq[id].type)
    $("#"+id+".selectpicker").val("default")
    $("#"+id+".selectpicker").selectpicker("refresh")
    data.seq[id].scale = "no scale selected";
    data.seq[id].name = "Seq" + (Number(id)+1);
    data.seq[id].name = addScaleToName(id,data.seq[id].name)
    $('#'+id+'.name').val(data.seq[id].name);

}

function download_data_Uniprot(id) {

    var xhttp = new XMLHttpRequest();
    var parser = new DOMParser();
    var db_ID = $('#'+id+'.sequence-id-input').val();
    var seqName = "";

    if (db_ID.length == 6 || db_ID.length == 8) {
        xhttp.open("GET", "https://www.ebi.ac.uk/proteins/api/proteins/" + db_ID, true);
        xhttp.setRequestHeader("Content-type", "application/json");
        xhttp.onload = function (e) {
            if (xhttp.readyState === 4) {
                if (xhttp.status === 200) {
                    response = JSON.parse(xhttp.responseText);
                    data.seq[id].sequence = response.sequence.sequence;
                    $('#'+id+'.form-inline').val(response.sequence.sequence)

                    try {
                        seqName = response.protein.recommendedName.fullName.value;
                    } catch (err) {
                        seqName = "*Name not found*";
                    }

                    data.seq[id].name = addScaleToName(id, seqName);
                    data.seq[id].window = 21;
                    data.seq[id].type = 'protein';
                    $('#'+id+'.name').val(data.seq[id].name);

                    $("#setup" + id + " .alert_box").empty();
                    $("#setup" + id + " .alert_box").append("<div class='alert alert-success alert-dismissable'><button type='button' class='close' data-dismiss='alert' aria-hidden='true'>&times;</button><strong>Success</strong> loading data for " + seqName + "</div>");

                    $("#"+id+".selectpicker").html(`
                                        <option value="p_alpha">Alpha-propensity</option>
                                        <option value="p_beta">Beta-propensity</option>
                                        <option value="p_charge">Charge</option>
                                        <option value="p_comp">Composition</option>
                                        <option value="p_hydro">Hydrophobicity</option>
                                        <option value="p_pchem">general Physico-Chemical</option>
                                        <option value="p_rnaAf">RNA-Affinity</option>
                                        <option value="p_other">Other</option>
                                        `)
                    $("#"+id+".selectpicker").selectpicker('refresh')

                    drawingProfiles(data);
                    updateShifts(data);
                }
                else {
                    $("#setup" + id + " .alert_box").empty();
                    $("#setup" + id + " .alert_box").append("<div class='alert alert-danger alert-dismissable'><button type='button' class='close' data-dismiss='alert' aria-hidden='true'>&times;</button><strong>Error</strong> could not locate given ID</div>");
                }
            }
        }
        xhttp.onerror = function(e) {
            $("#setup" + id + " .alert_box").empty();
            $("#setup" + id + " .alert_box").append("<div class='alert alert-danger alert-dismissable'><button type='button' class='close' data-dismiss='alert' aria-hidden='true'>&times;</button><strong>Error</strong> could not connect to Uniprot database</div>");
        }
        xhttp.send();
    }
    else {
        $("#setup" + id + " .alert_box").empty();
        $("#setup" + id + " .alert_box").append("<div class='alert alert-danger alert-dismissable'><button type='button' class='close' data-dismiss='alert' aria-hidden='true'>&times;</button><strong>Value Error:</strong> ID has not the correct length (6 or 8) of an UniprotID " + db_ID + "</div>");
    }
}

function download_data_ENA(id) {
    var xhttp = new XMLHttpRequest();
    var db_ID = $('#'+id+'.sequence-id-input').val();

    var rna_request = "https://www.ebi.ac.uk/ena/browser/api/fasta/" + db_ID
    xhttp.open("GET", rna_request, true);
    xhttp.setRequestHeader("Content-type", "text/plain");
    xhttp.onload = function (e) {
        if (xhttp.readyState === 4) {
            if (xhttp.status === 200) {
                try {
                    RNA_response = xhttp.responseText;
                    RNA = parseFASTA(RNA_response)


                    if (data.seq[id].type == "rna") {
                        RNA = RNA.replace(/T/g,'U');
                        $("#"+id+".selectpicker").html(`
                                        <option value="r_comp">Composition</option>
                                        `)
                    }
                    else {
                        $("#"+id+".selectpicker").html(`
                                        <option value="d_comp">Composition</option>
                                        `)
                    }
                    $("#"+id+".selectpicker").selectpicker('refresh')

                    data.seq[id].sequence = RNA;
                    $('#'+id+'.form-inline').val(RNA);

                    data.seq[id].window = 21;
                    data.seq[id].name = addScaleToName(id,db_ID);
                    
                    $('#'+id+'.name').val(data.seq[id].name);

                 
                    drawingProfiles(data);
                    updateShifts(data);
                  
                }
                catch(e) {
                    $("#setup" + id + " .alert_box").empty();
                    $("#setup" + id + " .alert_box").append("<div class='alert alert-danger alert-dismissable'><button type='button' class='close' data-dismiss='alert' aria-hidden='true'>&times;</button><strong>Value Error</strong> could not locate given ID in the European Nucleotide Archive (" + db_ID + ")</div>");
                  
                }
            }
            else {
                $("#setup" + id + " .alert_box").empty();
                $("#setup" + id + " .alert_box").append("<div class='alert alert-danger alert-dismissable'><button type='button' class='close' data-dismiss='alert' aria-hidden='true'>&times;</button><strong>Value Error</strong> could not locate given ID in the European Nucleotide Archive (" + db_ID + ")</div>");
               
            }
        }
    };
    xhttp.onerror = function (e) {
        $("#setup" + id + " .alert_box").empty();
        $("#setup" + id + " .alert_box").append("<div class='alert alert-danger alert-dismissable'><button type='button' class='close' data-dismiss='alert' aria-hidden='true'>&times;</button><strong>Value Error</strong> could not locate given ID in the European Nucleotide Archive (" + db_ID + ")</div>");
         
    };
    xhttp.send();
}

function removeInterfaceBox(id) {
    var myNode = $("#"+id+".interfacebox");
    $("#"+id+".interfacebox").remove();
    $("#legend"+id).remove();
    data.seq[id].active = false;
    data.seq[id].visible = false;
    data.seq[id].profile = null;
    drawingProfiles(data);
    updateShifts(data);
    updateCorrList(data);
}

function addLegend(id) {
    /* only add legend if a container has been created */
    if (document.getElementById('legend-container')) {
        var div = document.createElement('div');
        
        div.setAttribute('class', 'legendbox');
        div.setAttribute('id', "legend" + id);
        
        div.innerHTML = Legend({number: id});

        document.getElementById('legend-container').appendChild(div);

        $("#legend"+id+" .color-display").css('background-color', data.seq[id].color);
        updateLegend(id);
    }

}


function partey() {
    data.seq.forEach(function(e) {
        e.color = getRandomColor();
        });
    drawingProfiles(data);
}

(function(console){

console.save = function(data, filename){

    if(!data) {
        console.error('Console.save: No data')
        return;
    }

    if(!filename) filename = 'console.json'

    if(typeof data === "object"){
        data = JSON.stringify(data, undefined, 4)
    }

    var blob = new Blob([data], {type: 'text/json'}),
        e    = document.createEvent('MouseEvents'),
        a    = document.createElement('a')

    a.download = filename
    a.href = window.URL.createObjectURL(blob)
    a.dataset.downloadurl =  ['text/json', a.download, a.href].join(':')
    e.initMouseEvent('click', true, false, window, 0, 0, 0, 0, 0, false, false, false, false, 0, null)
    a.dispatchEvent(e)
 }
})(console)

function saveData(name) {
    if (name) {
        console.save(data, name);
    }
    else {
        var Dt = new Date(),
        d = Dt.getDate(),
        m = Dt.getMonth() +1,
        y = Dt.getFullYear(),
        h = Dt.getHours(),
        n = Dt.getMinutes(),
        s = Dt.getSeconds();
        var name = 'bugreport' + '_' + d + '-' + m + '-' + y + '-' + h + '-' + n + '-' + s + '.json'
        console.save(data, name);
    }  
    
}

function collapseAllInterface() {
    activeData = getActiveData(data)
    for (entry in activeData) {
        if($("div#"+activeData[entry][1]+".button").html() == '<i class="material-icons md-24">expand_less</i>'){
            $("div#"+activeData[entry][1]+".button").html('<i class="material-icons md-24">expand_more</i>');
            $("div#"+activeData[entry][1]+".button").parent().parent().find('.box').slideToggle(50, 'linear');
        }
    }
}

function addTable(target, data) {
    var targetDiv = $(target + ' .cor-content');
    var targetHeader = $(target + ' .title_bar span')
    var table = document.createElement('table');
    table.classList.add('correlation-table')

    var localData = getActiveData(data);
    titleString = "Correlation - not enough sequences";
    if (localData.length > 1) {
        titleString = "Correlation"
        bottomString = "";
        switch (data.cor_method.upper) {
            case "pearson":
                upper = calculate_pearson;
                bottomString += "upper: Pearson; "
                break;
            case "spearman":
                upper = calculate_spearman;
                bottomString += "upper: Spearman; "
                break;
            case "determination":
                upper = calculate_determination;
                bottomString += "upper: R<sup>2</sup>; "
                break;
            case "rmsd":
                upper = calculate_RMSD;
                bottomString += "upper: RMSD; "
                break;
            default:
                upper = calculate_pearson;
                bottomString += "upper: Pearson; "
        };
        switch (data.cor_method.lower) {
            case "pearson":
                lower = calculate_pearson;
                bottomString += "lower: Pearson"
                break;
            case "spearman":
                lower = calculate_spearman;
                bottomString += "lower: Spearman"
                break;
            case "determination":
                lower = calculate_determination;
                bottomString += "lower: R<sup>2</sup>; "
                break;
            case "rmsd":
                lower = calculate_RMSD;
                bottomString += "lower: RMSD"
                break;
            default:
                lower = calculate_spearman;
                bottomString += "lower: Spearman"
        };
        for (var i = -1; i < localData.length; i++) {
            var tr = table.insertRow();
            if (i>=0) {
                for (var j = -1; j < localData.length; j++) {
                    
                    if (j>=0) {
                        if (j<i) {
                            var td = tr.insertCell();
                            td.appendChild(document.createTextNode(lower(localData[i][1], localData[j][1])));
                            if (data.cor_method.lower == "rmsd"){
                                if (data.seq[i].scale != data.seq[j].scale || data.seq[i].type != data.seq[j].type) {
                                    td.classList.add('nonsense');
                                }
                            }
                        }
                        if (j==i) {
                            var td = tr.insertCell();
                            td.appendChild(document.createTextNode("1"));
                            td.classList.add('identity')
                        }
                        if (j>i) {
                            var td = tr.insertCell();
                            td.appendChild(document.createTextNode(upper(localData[i][1], localData[j][1])));
                            if (data.cor_method.upper == "rmsd"){
                                if (data.seq[i].scale != data.seq[j].scale || data.seq[i].type != data.seq[j].type) {
                                    td.classList.add('nonsense');
                                }
                            }
                        }
                    }
                    else {
                        var th = tr.insertCell();
                        th.outerHTML = "<th>"+String(localData[i][1]+1)+"</th>";
                    }
                }
            }
            else {
                for (var j = -1; j < localData.length; j++) {
                    if (j>=0) {
                        var th = tr.insertCell();
                        th.outerHTML = "<th>"+String(localData[j][1]+1)+"</th>";
                    }
                    else {
                        var th = tr.insertCell();
                        th.outerHTML = "<th></th>";
                    }
                }
            }
        }
        
        targetDiv.html(table);
        var bottomDiv = document.createElement("div");
        bottomDiv.classList.add('cor_bottom_text');
        if (data.cor_method.upper == "rmsd" || data.cor_method.lower == "rmsd") {
            bottomString += "&#10;RMSD is only meaningful for identical scales"
        }
        bottomDiv.innerHTML = bottomString;
        targetDiv.append(bottomDiv);
        
    }
    else {
        targetDiv.html("");
    }
    targetHeader.text(titleString);
}

function saveCSV(dat, id) {
    console.save(createCSVfromArray(getArraysFromProfile(dat, id)), id+".csv");
}

function saveMultipleCSV(dat, idArray) {
    for (var i=0; i< idArray.length; i++) {
        idArray[i] = Number(idArray[i]);
    }
    var localData = {}
    localData[0] = $.extend(true, {}, dat.seq[idArray[0]].profile);
    var xArray = getArraysFromProfile(dat, idArray[0]).x;
    for (var i=1; i<idArray.length; i++) {
        xArray = xArray.concat(getArraysFromProfile(dat, idArray[i]).x);
        localData[i] = $.extend(true, {}, dat.seq[idArray[1]].profile);
    }
    xArray = getUniqueFromArray(xArray).sort(function(a,b) {return a-b});
    
    csvData = "x,";
    tmpData = [];
    for (var i=0; i<idArray.length; i++) {
        tmpData.push("y"+idArray[i]);
    }
    csvData += tmpData.join() + "\r\n";
    
    for (var i=0; i<xArray.length; i++) {
        tmpData = [];
        for (var j=-1; j<idArray.length; j++) {
            if (j==-1) {
                tmpData.push(xArray[i]);
            }
            else {
                for (var k=0; k<Object.keys(localData[j]).length; k++) {
                    if (localData[j][k].x == xArray[i]) {
                        tmpData.push(localData[j][k].y);
                        break;
                    }
                }
                if (!tmpData[j+1]) {
                    tmpData.push(NaN);
                }
            }
        }
        csvData += tmpData.join().replace(/NaN/gi,'') + "\r\n";
    }

    console.save(csvData, "visualizationData.csv");
}