function updateCorResults(corList) {
    for(corID in corList) {
        cor = calculate_correlation(corID, corList[corID]);
        cor = "R = " + cor;
        $("#"+corID+".corr_output").text(cor);
    }    
}

function findBestShift(id1, id2, inv=false) {
    f1 = 1;
    f2 = 1;
    if (data.seq[id1].type == "rna") {f1 = 3};
    if (data.seq[id2].type == "rna") {f2 = 3};
    realLength1 = data.seq[id1].profile.length + data.seq[id1].window/f1 - 1;
    realLength2 = data.seq[id2].profile.length + data.seq[id2].window/f2 - 1;
    data.seq[id1].shift = 0;
    data.seq[id2].shift = 0;
    drawingProfiles(data);

    if (realLength1 > realLength2) {
        correlations = []
        maxShift = realLength1 - realLength2;
        for (var i=0; i<maxShift; i++) {
            data.seq[id2].shift = i;
            drawingProfiles(data);
            correlations.push(Number(calculate_correlation(id1, id2)));
        }
        if (inv) { data.seq[id2].shift = correlations.indexOf(d3.min(correlations)); }
        else { data.seq[id2].shift = correlations.indexOf(d3.max(correlations)); }
        drawingProfiles(data);
        updateShifts(data);
        // 1 is longer than 2 => slide 2
    }
    else if (realLength1 < realLength2) {
        correlations = []
        maxShift = realLength2 - realLength1;
        for (var i=0; i<maxShift; i++) {
            data.seq[id1].shift = i;
            drawingProfiles(data);
            correlations.push(Number(calculate_correlation(id1, id2)));
        }
        if (inv) { data.seq[id1].shift = correlations.indexOf(d3.min(correlations)); }
        else { data.seq[id1].shift = correlations.indexOf(d3.max(correlations)); }        
        drawingProfiles(data);
        updateShifts(data);
        // 2 is longer than 1 = > slide 1
    }
    else {
        console.log("No shift possible due to same length")
    }
}

function updateLegend(id) {
    $("#legend"+id+" .legend-text").text(getPropertyString(id));
    
}
function profile_calculation(seq, wndw, scaleID, type, shift, method) {
    var seqArray = [];
    var smoothArray = [];
    if (scaleID == "") {
        for (var i = 0, len = seq.length; i < len; i++) {
            seqArray.push(0);
        };  
    }
    else {
        for (var i = 0, len = seq.length; i < len; i++) {
            seqArray.push(scale[type][scaleID][seq[i]]);
        };
    }
    switch(method) {
        case "boxcar":
            switch (type) {
                case "protein":
                    smoothArray = windowAverage(seqArray,wndw,1);
                    var xVals = [((wndw-1)/2)+shift], x=(wndw-1)/2;
                    while(x<(smoothArray.length+((wndw-1)/2) ) ){x+=1;xVals.push(x+shift)};
                    break;
                case "rna":
                    smoothArray = windowAverage(seqArray,wndw,3)
                    var xVals = [(((wndw/3)-1)/2)+shift], x=((wndw/3)-1)/2;
                    while(x<(smoothArray.length+(((wndw/3) -1)/2) ) ){x+=1; xVals.push(x+shift)};
                    break;
            };
        break;
        case "savitzky-golay":
            switch (type) {
                case "protein":
                    smoothArray = savitzky_golay(seqArray,wndw,1);
                    var xVals = [((wndw-1)/2)+shift], x=(wndw-1)/2;
                    while(x<(smoothArray.length+((wndw-1)/2) ) ){x+=1;xVals.push(x+shift)};
                    break;
                case "rna":
                    smoothArray = savitzky_golay(seqArray,wndw,3)
                    var xVals = [(((wndw/3)-1)/2)+shift], x=((wndw/3)-1)/2;
                    while(x<(smoothArray.length+(((wndw/3) -1)/2) ) ){x+=1; xVals.push(x+shift)};
                    break;
            };

    };

    var data = []
    for (var i = 0; i < smoothArray.length; i++) {
        data.push({"x":xVals[i], "y":smoothArray[i]});
    };
    return data;
};

function drawingProfiles(data, id = ".canvas") {
    
    data.seq.forEach(function(d) {
        if (d.active || d.visible) {
            d.profile = profile_calculation(d.sequence, d.window, d.scale, d.type, d.shift, d.smoothing_method);
        }
    });
    
    input = {"seq": []};
    input.title = data.title;
    input.seq = data.seq.filter(function (n) {
        return n.visible && n.active;
    });
    d3.selectAll("svg > *").remove();
    // d3.select(".canvas_g").remove();
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
    
    g.append("text")
        .attr("y", height + margin.bottom)
        .attr("x", width* (2/5))
        .text("Residue / Codon")

    if (input.seq.map(function(d) {return d["type"]}).includes("protein")) {

        var y1_min = findMin(data.seq, "protein");
        var y1_max = findMax(data.seq, "protein");
        var padding1 = (y1_max - y1_min)*0.1;

        var y1 = d3.scaleLinear()
                   .domain([y1_min - padding1, y1_max + padding1])
                   .rangeRound([height, 0]);

        var yAxis1 = d3.axisLeft(y1);
        
        var protLine = d3.line()
            .x(function(d) {return x1(d.x); })
            .y(function(d) {return y1(d.y); });

        var gY1 = g.append("g")
                   .classed("axisBlack", true)
                   .classed("yAxis", true)
                   .classed("axis", true)
                   .call(yAxis1);
        g.append("text")
            .attr("transform", "rotate(-90)")
            .attr("y", 0 - margin.left*(2/3))
            .attr("x", 0 - height *(3/5))
            .text("Protein")
    }
    if (input.seq.map(function(d) {return d["type"]}).includes("rna")) {

        var y2_min = findMin(input.seq, "rna");
        var y2_max = findMax(input.seq, "rna");
        padding2 = (y2_min - y2_max)*0.1;

        var y2 = d3.scaleLinear()
                   .domain([y2_max - padding2, y2_min + padding2])
                   .rangeRound([height, 0]);

        var yAxis2 = d3.axisRight(y2);

        var rnaLine  = d3.line()
            .x(function(d) {return x1(d.x); })
            .y(function(d) {return y2(d.y); });

        var gY2 = g.append("g")
                   .classed("axisBlack", true)
                   .classed("yAxis", true)
                   .classed("axis", true)
                   .attr("transform", "translate(" + width + ",0)")
                   .call(yAxis2);
        g.append("text")
            .attr("y", -width - margin.right * (2/3))
            .attr("x", height*(2/5))
            .attr("transform", "rotate(90)")
            .text("RNA")
    }

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
        var xz = d3.event.transform.rescaleX(x1);
        // var xy1 = d3.event.transform.rescaleY(y1);
        // var xy2 = d3.event.transform.rescaleY(y2);

        gX1.call(xAxis.scale(xz));
        // gY1.call(yAxis1.scale(xy1));
        // gY2.call(yAxis2.scale(xy2));
        if (input.seq.map(function(d) {return d["type"]}).includes("rna")) {
            rnaLine
                .x(function(d) {return xz(d.x); })
                // .y(function(d) {return xy2(d.y); })
        }
        if (input.seq.map(function(d) {return d["type"]}).includes("protein")) {
            protLine
                .x(function(d) {return xz(d.x); })
                // .y(function(d) {return xy1(d.y); })
        }
        d3.selectAll(".protein").attr("d", protLine)
        // pL.attr("d", protLine)
        d3.selectAll(".rna").attr("d", rnaLine)

    }

    function drawLines(item, index) {
        g = d3.select(".canvas_g")
        if (item.type == 'rna') {
            g.append("path")
             .classed('rna', true)
             .datum(item.profile)
             .attr("clip-path", "url(#clip)")
             .attr("fill", "none")
             .attr("stroke", item.color)
             .attr("stroke-linejoin", "round")
             .attr("stroke-linecap", "round")
             .attr("stroke-width", item.thickness)
             .attr("d", rnaLine);
        } else if (item.type == 'protein') {
            g.append("path")
             .classed('protein', true)
             .datum(item.profile)
             .attr("clip-path", "url(#clip)")
             .attr("fill", "none")
             .attr("stroke", item.color)
             .attr("stroke-linejoin", "round")
             .attr("stroke-linecap", "round")
             .attr("stroke-width", item.thickness)
             .attr("d", protLine);
        }
    };

    input.seq.forEach(drawLines);
    updateCorResults(corList);
    // addInterface(data);

    return "complete";

    
};

const Item = ({number}) =>`
<div class="sequence">
    <div class="title_bar">
        <div class="name-value">${number+1}:</div>
        <input class="name" type="text" id="${number}" value="Seq${number+1}" autocomplete="off" autocorrect="off" autocapitalize="off" spellcheck="false"></input>
        <div class="close-button" data-toggle="tooltip" data-placement="right" title="Remove data" id="${number}"><i class="material-icons md-36">close</i></div>
        <div class="button" id="${number}"><i class="material-icons md-36">expand_less</i></div>
        <div class="visible" id="${number}"><i class="material-icons md-24" style="margin-top: 6px">visibility</i></div>
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
        <div class="secondRow">
            <div class="btn-group" style="float: left;">
                <button class="btn btn-grey btn-sm btn-type" id="${number}" type="button">
                    Type
                </button>
                <button type="button" class="btn btn-sm btn-grey dropdown-toggle dropdown-toggle-split" data-toggle="dropdown" aria-haspopup="true" aria-expanded="false">
                    <span class="sr-only">Toggle Dropdown</span>
                </button>
                <div class="dropdown-menu type-selector" id="${number}">
                    <a class="dropdown-item">Protein</a>
                    <a class="dropdown-item">RNA</a>
                </div>
            </div>

            <button type="button" class="btn btn-grey btn-sm btn-sequence-input" data-toggle="modal" data-target="#modal${number}">Sequence</button>
            <div class="modal" id="modal${number}">
                <div class="modal-dialog modal-lg">
                    <div class="modal-content">
                        <div class="modal-header">
                            <h4 class="modal-title">Sequence</h4>
                            <button type="button" class="close" data-dismiss="modal">&times;</button>
                        </div>
                        <div class="modal-body">
                            <textarea class="form-inline" id="${number}" type="text" size="300" style="display: block; width: 100%; height: 8em; resize: none; font-family: 'Source Code Pro', monospace;"></textarea>
                        </div>
                        <div class="modal-footer">
                            
                            <div class="mr-auto">
                                <div class="btn-group database-menu" style="float: left;">
                                    <button class="btn btn-grey btn-sm btn-label" id="${number}" type="button">
                                        Database
                                    </button>
                                    <button type="button" class="btn btn-sm btn-grey dropdown-toggle dropdown-toggle-split btn-flat" data-toggle="dropdown" aria-haspopup="true" aria-expanded="false">
                                        <span class="sr-only">Toggle Dropdown</span>
                                    </button>
                                    <div class="dropdown-menu database-selector" id="${number}">
                                        <a class="dropdown-item">Uniprot</a>
                                        <a class="dropdown-item">ENA</a>
                                    </div>
                                </div>
                                <input class="sequence-id-input" id="${number}" style="margin-right: 7px;" type="text" placeholder="Enter sequence ID"/>
                                <button type="button" class="btn btn-sm btn-blue btn-submit btn-flat" id="${number}" style="float: right;">Submit</button>
                            </div>
                            <button type="button" class="btn btn-sm btn-red btn-clr btn-flat" id="${number}">Clear</button>
                            
                        </div>
                        <div class="modal-footer">
                            <div class="alert_box" style="display: inline-block; width: 100%; min-width: 350px; font-size: 15px"></div>
                        </div>
                    </div>
                </div>
            </div>

            <div class="btn-group" style="float: left;">
                <button class="btn btn-grey btn-sm btn-scale" id="${number}" type="button">
                    Scale
                </button>
                <button type="button" class="btn btn-sm btn-grey dropdown-toggle dropdown-toggle-split" data-toggle="dropdown" aria-haspopup="true" aria-expanded="false">
                    <span class="sr-only">Toggle Dropdown</span>
                </button>
                <div class="dropdown-menu scale-selector" id="${number}">
                    <a class="dropdown-item" id="ADE">ADE</a>
                   <a class="dropdown-item" id="CYT">CYT</a>
                   <a class="dropdown-item" id="GUA">GUA</a>
                   <a class="dropdown-item" id="URA">URA</a>
                </div>
            </div>
            

            <button class="color-picker-btn" id="colorpicker${number}"></button>
            
            <div class="btn-group" style="float: left;">
                <button class="btn btn-grey btn-sm" type="button">
                    Compare
                </button>
                <button type="button" class="btn btn-sm btn-grey dropdown-toggle dropdown-toggle-split" data-toggle="dropdown" aria-haspopup="true" aria-expanded="false">
                    <span class="sr-only">Toggle Dropdown</span>
                </button>
                <div class="dropdown-menu corr-selector" id="${number}">
                </div>
            </div>

            <label id="${number}" class="corr_output"></label>
            
            
                        
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
                  "scale" : "GUA",
                  "sequence" : "",
                  "window" : 21,
                  "shift" : 0,
                  "smoothing_method" : "boxcar",
                  "name" : "",
                  "organism" : "",
                  "thickness" : 2,
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
    // div.className = 'interfacebox';
    div.innerHTML = Item({number: num});

    data.seq[num] = $.extend(true, {}, data_prototype);
    data.seq[num].sequence = "";
    // data.seq[num].sequence = getRandomSequence();
    data.seq[num].color = getBiasedColorHSL(num, .60, .60, .2);
    // data.seq[num].color = getRandomColorHSL(null, .40, .50, .3);
    // data.seq[num].color = getRandomColor();
    data.seq[num].name = "Seq" + (num+1);
    
    document.getElementById('interface_container').appendChild(div);

    new jscolor(document.getElementById('colorpicker' + num),
                {valueElement:'none', onFineChange:'updateColor(this, ' + num + ')', value:data.seq[num].color.slice(1)})

    window_p = $('.prot-window').each(function(element, value) {
        $(this).slider()
        .on("slide", function() {
            if (data.seq[this.id].type == "protein") {
                data.seq[this.id].window = Number(this.value);
                $("#"+this.id+".firstRow-text").text("Window: " + this.value);
            } else {
                data.seq[this.id].window = 3*Number(this.value);
                $("#"+this.id+".firstRow-text").text("Window: " + 3*Number(this.value));
            }
            drawingProfiles(data);
        })
        .on("slideStop", function() {
            if (data.seq[this.id].type == "protein") {
                data.seq[this.id].window = Number(this.value);
                $("#"+this.id+".firstRow-text").text("Window: " + this.value);
            } else {
                data.seq[this.id].window = 3*Number(this.value);
                $("#"+this.id+".firstRow-text").text("Window: " + 3*Number(this.value));
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

    $('.form-inline').each(function(element, value) {
        $(this).val(data.seq[this.id].sequence);
        });

    $('.form-inline').each(function(element, value) {
        $(this).bind('input propertychange', function() {
            data.seq[this.id].sequence = parseFASTA(this.value);
            drawingProfiles(data);
        });
    });
    $('.name').each(function(element, value) {
        $(this).on("change", function() {
            data.seq[this.id].name = this.value;
            updateCorrList(data);
            updateLegend(this.id);
        })
    });

    $("div#"+num+".button").click(function(){
        if($(this).html() == '<i class="material-icons md-36">expand_less</i>'){
            $(this).html('<i class="material-icons md-36">expand_more</i>');
        }
        else{
            $(this).html('<i class="material-icons md-36">expand_less</i>');
        }
        $(this).parent().parent().find('.box').slideToggle();
    });
    $("div#"+num+".visible").click(function(){
        if($(this).html() == '<i class="material-icons md-24" style="margin-top: 6px">visibility</i>'){
            $(this).html('<i class="material-icons md-24" style="margin-top: 6px">visibility_off</i>');
            data.seq[this.id].visible = false;
            drawingProfiles(data);
        }
        else{
            $(this).html('<i class="material-icons md-24" style="margin-top: 6px">visibility</i>');
            data.seq[this.id].visible = true;
            drawingProfiles(data);
        }
        
    });

    $(".close-button").tooltip();

    $("#modal"+num).modal('show');

    $("#modal"+num).on('hidden.bs.modal', function() {
        $("#modal"+num + " .alert_box").empty();
    });


    
    drawingProfiles(data);
    updateShifts(data);
    updateInterface(data);
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
        addLegend(num);

        $("div#"+num+".button").click(function(){
            if($(this).html() == '<i class="material-icons md-36">expand_less</i>'){
                $(this).html('<i class="material-icons md-36">expand_more</i>');
            }
            else{
                $(this).html('<i class="material-icons md-36">expand_less</i>');
            }
            $(this).parent().parent().find('.box').slideToggle();
        });
        $("div#"+num+".visible").click(function(){
            if($(this).html() == '<i class="material-icons md-24" style="margin-top: 6px">visibility</i>'){
                $(this).html('<i class="material-icons md-24" style="margin-top: 6px">visibility_off</i>');
                data.seq[this.id].visible = false;
                drawingProfiles(data);
            }
            else{
                $(this).html('<i class="material-icons md-24" style="margin-top: 6px">visibility</i>');
                data.seq[this.id].visible = true;
                drawingProfiles(data);
            }
            
        });
    }

    window_p = $('.prot-window').each(function(element, value) {
        $(this).slider()
        .on("slide", function() {
            if (data.seq[this.id].type == "protein") {
                data.seq[this.id].window = Number(this.value);
                $("#"+this.id+".firstRow-text").text("Window: " + this.value);
            } else {
                data.seq[this.id].window = 3*Number(this.value);
                $("#"+this.id+".firstRow-text").text("Window: " + 3*Number(this.value));
            }
            drawingProfiles(data);
        })
        .on("slideStop", function() {
            if (data.seq[this.id].type == "protein") {
                data.seq[this.id].window = Number(this.value);
                $("#"+this.id+".firstRow-text").text("Window: " + this.value);
            } else {
                data.seq[this.id].window = 3*Number(this.value);
                $("#"+this.id+".firstRow-text").text("Window: " + 3*Number(this.value));
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

    $('.form-inline').each(function(element, value) {
        $(this).val(data.seq[this.id].sequence);
        });

    $('.form-inline').each(function(element, value) {
        $(this).bind('input propertychange', function() {
            data.seq[this.id].sequence = parseFASTA(this.value);
            drawingProfiles(data);
        });
    });
    $('.name').each(function(element, value) {
        $(this).on("change", function() {
            data.seq[this.id].name = this.value;
            updateCorrList(data);
            updateLegend(this.id);
        })
    });

    $(".close-button").tooltip();

    drawingProfiles(data);
    updateShifts(data);
    updateWindows(data);
    updateInterface(data);

}

function updateInterface(data) {
    $(".scale-selector").empty();
    for (var i = data.seq.length - 1; i >= 0; i--) {
        if (data.seq[i].active) {
            if (data.seq[i].type.toLowerCase() == 'rna') {
                $("#"+ i +".firstRow-text").text("Window: " + 3 * Number($("#" + i + ".interfacebox .firstRow-slider .slider").text()));
                data.seq[i].window = 3 * Number($("#" + i + ".interfacebox .firstRow-slider .slider").text());
                $("#"+i+".scale-selector").append(`<a class="dropdown-item" id="ADE">ADE-content</a>
                                                <a class="dropdown-item" id="CYT">CYT-content</a>
                                                <a class="dropdown-item" id="GUA">GUA-content</a>
                                                <a class="dropdown-item" id="URA">URA-content</a>
                                                <a class="dropdown-item" id="PUR">PUR-content</a>
                                                <a class="dropdown-item" id="PYR">PYR-content</a>`);
                                                // <a class="dropdown-item" id="voro_interface">Protein Interface (Voro2018)</a>`)
                if  (!['ADE', 'CYT', 'GUA', 'URA', 'PUR', 'PYR'].includes(data.seq[i].scale)) {
                    data.seq[i].scale = "";
                }
            } else if (data.seq[i].type.toLowerCase() == 'protein') {
                $("#"+ i +".firstRow-text").text("Window: " + Number($("#" + i + ".interfacebox .firstRow-slider .slider").text()));
                data.seq[i].window = Number($("#" + i + ".interfacebox .firstRow-slider .slider").text());
                $("#"+i+".scale-selector").append(`<a class="dropdown-item" id="ADE">ADE-affinity (NAR2013)</a>
                                                <a class="dropdown-item" id="CYT">CYT-affinity (NAR2013)</a>
                                                <a class="dropdown-item" id="GUA">GUA-affinity (NAR2013)</a>
                                                <a class="dropdown-item" id="URA">URA-affinity (NAR2013)</a>
                                                <a class="dropdown-item" id="PUR">PUR-affinity (NAR2013)</a>
                                                <a class="dropdown-item" id="FAC1">Hydrophobicity (Factor 1)</a>
                                                <a class="dropdown-item" id="Charge">Charge</a>
                                                <a class="dropdown-item" id="PR">Polar Requirement (Woese 1973)</a>
                                                <a class="dropdown-item" id="PR_2008">Polar Requirement (Mathew 2008)</a>`)
                                                // <a class="dropdown-item" id="voro003_ADE">ADE-affinity (Voro2018)</a>
                                                // <a class="dropdown-item" id="voro003_CYT">CYT-affinity (Voro2018)</a>
                                                // <a class="dropdown-item" id="voro003_GUA">GUA-affinity (Voro2018)</a>
                                                // <a class="dropdown-item" id="voro003_URA">URA-affinity (Voro2018)</a>
                                                // <a class="dropdown-item" id="voro003_RNAgeneral">RNA Interface (Voro2018)</a>
                                                
                if (!['ADE', 'CYT', 'GUA', 'URA', 'PUR', 'FAC1', 'Charge', 'PR'].includes(data.seq[i].scale)) {
                    data.seq[i].scale = "";
                }
            } else {
                console.error("TYPE ERROR: data.seq["+i+"].type = " + data.seq[i].type)
            }
            updateCorrList(data);

            if (data.seq[i].scale) {
                $("#" + i + ".btn-scale").html(trimString(scaleNames[data.seq[i].type][data.seq[i].scale], 9)); 
            }
            else {
                $("#" + i + ".btn-scale").html(trimString("No Scale", 9));
            }
            $("#" + i + ".btn-type").html(trimString(typeNames[data.seq[i].type], 9));

            updateLegend(i);

        }
    }
}

function updateCorrList(data) {
    $(".corr-selector").empty();
    activeData = getActiveData(data);
    for (var i = data.seq.length - 1; i >= 0; i--) {
        $("#"+i+".corr-selector").append(getCorOptions(activeData, i));
    }
}

function updateShifts(data) {
    $(".shift").slider('setAttribute', 'max', findMaxLength(data.seq)-3);
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

function download_data_Uniprot(id) {

    var xhttp = new XMLHttpRequest();
    var parser = new DOMParser();
    var db_ID = $('#'+id+'.sequence-id-input').val();
    var seqName = "";

    if (db_ID.length == 6) {
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

                    data.seq[id].name = seqName;
                    data.seq[id].window = 21;
                    data.seq[id].type = 'protein';
                    $('#'+id+'.name').val(seqName);

                    $("#modal" + id + " .alert_box").append("<div class='alert alert-success alert-dismissable'><button type='button' class='close' data-dismiss='alert' aria-hidden='true'>&times;</button><strong>Success</strong> loading data for " + seqName + "</div>");

                    drawingProfiles(data);
                    updateInterface(data);
                    updateShifts(data);
                }
                else {
                    $("#modal" + id + " .alert_box").append("<div class='alert alert-danger alert-dismissable'><button type='button' class='close' data-dismiss='alert' aria-hidden='true'>&times;</button><strong>Error</strong> could not locate given ID</div>");
                }
            }
        }
        xhttp.onerror = function(e) {
            $("#modal" + id + " .alert_box").append("<div class='alert alert-danger alert-dismissable'><button type='button' class='close' data-dismiss='alert' aria-hidden='true'>&times;</button><strong>Error</strong> could not connect to Uniprot database</div>");
        }
        xhttp.send();
    }
    else {
        $("#modal" + id + " .alert_box").append("<div class='alert alert-danger alert-dismissable'><button type='button' class='close' data-dismiss='alert' aria-hidden='true'>&times;</button><strong>Value Error:</strong> ID has not the correct length (6) if an UniprotID " + db_ID + "</div>");
    }
}

function download_data_ENA(id) {

    var xhttp = new XMLHttpRequest();
    var parser = new DOMParser();
    var db_ID = $('#'+id+'.sequence-id-input').val();
        
    var rna_request = "https://www.ebi.ac.uk/ena/data/view/" + db_ID.split('.')[0] + "&display=xml"
    xhttp.open("GET", rna_request, true);
    xhttp.setRequestHeader("Content-type", "text/xml");
    xhttp.onload = function (e) {
        if (xhttp.readyState === 4) {
            if (xhttp.status === 200) {
                try {
                    RNA_response = parser.parseFromString(xhttp.responseText, "text/xml");
                    RNA = RNA_response.getElementsByTagName("sequence")[0].childNodes[0].nodeValue;
                    RNA = RNA.replace(/\s/g,'');
                    RNA = RNA.toUpperCase();
                    RNA = RNA.replace(/T/g,'U');

                    data.seq[id].sequence = RNA;
                    $('#'+id+'.form-inline').val(RNA);
                    
                    data.seq[id].name = db_ID;
                    data.seq[id].window = 63;
                    data.seq[id].type = 'rna';
                    $('#'+id+'.name').val(db_ID);

                    drawingProfiles(data);
                    updateInterface(data);
                    updateShifts(data);
                }
                catch(e) {
                    $("#modal" + id + " .alert_box").append("<div class='alert alert-danger alert-dismissable'><button type='button' class='close' data-dismiss='alert' aria-hidden='true'>&times;</button><strong>Value Error</strong> could not locate given ID in the European Nucleotide Archive (" + db_ID + ")</div>");
                }
            }
            else {
                $("#modal" + id + " .alert_box").append("<div class='alert alert-danger alert-dismissable'><button type='button' class='close' data-dismiss='alert' aria-hidden='true'>&times;</button><strong>Value Error</strong> could not locate given ID in the European Nucleotide Archive (" + db_ID + ")</div>");
            }
        }
    };
    xhttp.onerror = function (e) {
        $("#modal" + id + " .alert_box").append("<div class='alert alert-danger alert-dismissable'><button type='button' class='close' data-dismiss='alert' aria-hidden='true'>&times;</button><strong>Value Error</strong> could not locate given ID in the European Nucleotide Archive (" + db_ID + ")</div>");
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
    /* only fire if a container has been created*/
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
    activeDatat = getActiveData(data)
    for (entry in activeData) {
        if($("div#"+activeData[entry][1]+".button").html() == '<i class="material-icons md-36">expand_less</i>'){
            $("div#"+activeData[entry][1]+".button").html('<i class="material-icons md-36">expand_more</i>');
            $("div#"+activeData[entry][1]+".button").parent().parent().find('.box').slideToggle();
        }
    }
}
// var zoom = d3.zoom()
//     .scaleExtent([0.25, 8])
//     .translateExtent([[-width, -height], [2 * width, 2 *height]])
//     .on("zoom", zoomed);

// var zoomRect = svg.append("rect")
//     .attr("width", width)
//     .attr("height", height)
//     .attr("fill", "none")
//     .attr("pointer-events", "all")
//     .attr("transform", "translate(" + margin.left +"," + margin.top +")")
//     .call(zoom);
