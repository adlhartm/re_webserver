function xVal(arr) {
    n = arr.length;
    xVals = new Array(n);
    for (i=0;i<n;i++) {
        xVals[i] = arr[i].x;
    }

    return xVals;
}

function yVal(arr) {
    n = arr.length;
    yVals = new Array(n);
    for (i=0;i<n;i++) {
        yVals[i] = arr[i].y;
    }

    return yVals;
}

function findOverlap(arr1, arr2) {
// not used anymore - replaced by calculateOverlap(prof1, prof2, shift1, shift2) for compatibility with gaps
    if (arr1.length >= 1 && arr2.length >= 1) {
        if (arr1[0].x < arr2[0].x) {
            r_min = 0;
            p_min = arr2[0].x - arr1[0].x;
        } else if (arr1[0].x > arr2[0].x) {
            r_min = arr1[0].x - arr2[0].x;
            p_min = 0;
        } else {
            r_min = 0;
            p_min = 0;
        }

        if (arr1[arr1.length-1].x < arr2[arr2.length-1].x) {
            r_max = arr2.length - (arr2[arr2.length-1].x - arr1[arr1.length-1].x);
            p_max = arr1.length;
        } else if (arr1[arr1.length-1].x > arr2[arr2.length-1].x) {
            r_max = arr2.length;
            p_max = arr1.length - (arr1[arr1.length-1].x - arr2[arr2.length-1].x);
        } else {
            r_max = arr2.length;
            p_max = arr1.length;
        }

        return [p_min, p_max, r_min, r_max]
    }
    else {
        return [0,0,0,0];
    }

}

function calculateOverlap(prof1, prof2, shift1, shift2) {
    p1 = prof1.slice();
    p2 = prof2.slice();
    var newProf1 = [];
    var newProf2 = [];
    j=0;
    loopLength = p1.length;
    for (var i=0; i<loopLength; i++) {
        while((p1[i].x > p2[j].x) && (j < p2.length-1)) {
            j+=1;
        }
        while((p1[i].x < p2[j].x) && (i<p1.length-1)) {
            i++;
        }
        // do not ignore y-values of 0
        if ((!isNaN(p1[i].y) && !isNaN(p2[j].y)) && (p1[i].x == p2[j].x)) {
            newProf1.push(p1[i]);
            newProf2.push(p2[j]);
        }
    }

    return {p1: newProf1, p2: newProf2}
}

function windowAverage(arr, weight, win, step) {
    var output = {"arr":[],
                  "weight":[]
    };
    win = win*step;
    if (arr.length >= win) {
        tmp = arr.splice(0,win);
        tmp_weight = weight.splice(0, win);
        output.arr.push(d3.sum(tmp)/d3.sum(tmp_weight))
        output.weight.push(d3.sum(tmp_weight)/win)
    } else {
        return output;
    }
    while (arr.length >= step) {
        tmp.splice(0,step);
        tmp = tmp.concat(arr.splice(0,step));
        tmp_weight.splice(0, step);
        tmp_weight = tmp_weight.concat(weight.splice(0,step));
        output.arr.push(d3.sum(tmp)/d3.sum(tmp_weight));
        output.weight.push(d3.sum(tmp_weight)/win)
    }

    return output;
}

function savitzky_golay(arr, win, step) {
    // just a placeholder currently...
    return windowAverage(arr, win, step);
}

function translate(rnaSeq, shift = 0) {
    let i=shift;
    let prot = '';
    while (i<rnaSeq.length-2){
        prot += code[rnaSeq.slice(i,i+3)];
        i += 3;
    }
    return prot;
}

function RMSD(seq1, seq2) {
    let RMSD = 0;
    for (let i=0; i<seq1.length; i++) {
        RMSD += (seq1[i] - seq2[i])**2;
    }
    RMSD = RMSD / seq1.length;
    RMSD = Math.sqrt(RMSD);
    return RMSD;
}

function getUniqueFromArray(array) {
    var a = array.concat();
    for (var i=0; i<a.length; ++i) {
        for (var j=i+1; j<a.length; ++j) {
            if (a[i] === a[j]) {
                a.splice(j--, 1);
            }
        }
    }
    return a;
}

function hasGaps(profile) {
    for (var i=0; i<profile.length; i++) {
        if (profile[i].w != 1) {
            return true
        };
    };
    return false
}

function findMin(dat, type) {
    tmplst = []
    for (var i=0; i<dat.length; i++) {
        if (dat[i].type == type & dat[i].active & dat[i].visible) {
            tmplst.push(d3.min(dat[i].profile, function(d) {return d.y}));
        };
    };
    return d3.min(tmplst);
};

function findMax(dat, type) {
    tmplst = []
    for (var i=0; i<dat.length; i++) {
        if (dat[i].type == type & dat[i].active & dat[i].visible) {
            tmplst.push(d3.max(dat[i].profile, function(d) {return d.y}));
        };
    };
    return d3.max(tmplst);
};

function findMinGlobal(dat) {
    tmplst = []
    for (var i=0; i<dat.length; i++) {
        if (dat[i].active & dat[i].visible) {
            tmplst.push(d3.min(dat[i].profile, function(d) {return d.y}));
        };
    };
    return d3.min(tmplst);
};

function findMaxGlobal(dat) {
    tmplst = []
    for (var i=0; i<dat.length; i++) {
        if (dat[i].active & dat[i].visible) {
            tmplst.push(d3.max(dat[i].profile, function(d) {return d.y}));
        };
    };
    return d3.max(tmplst);
};

function findMean(dat) {
    if (dat.active & dat.visible) {
        return d3.mean(dat.profile, function(d) {return d.y});
    } else {
        return null;
    };
};

function findDeviation(dat) {
    if (dat.active & dat.visible) {
        return d3.deviation(dat.profile, function(d) {return d.y});
    } else {
        return null;
    };
};

function rescale(dat, type){
    for (var i=0; i<dat.length; i++) {
        if (dat[i].type == type & dat[i].active & dat[i].visible) {
            mean = findMean(dat[i]);
            sd = findDeviation(dat[i]);
            for(var j=0; j<dat[i].profile.length; j++) {
                dat[i].profile[j].y = (dat[i].profile[j].y - mean)/sd;
            }
        } 
    }
};

function getArraysFromProfile(dat, id) {
    p = dat.seq[id].profile;
    xA = [];
    yA = [];
    for (d in p) {
        xA.push(p[d].x);
        yA.push(p[d].y);
    }
    return {'x':xA, 'y':yA};
}

function createCSVfromArray(xyDat) {
    csv = "x,y\r\n";
    for (var i=0; i<xyDat.x.length; i++) {
        csv += xyDat.x[i] + "," + xyDat.y[i] + "\r\n";
    }
    return csv;
}

function findMaxLength(dat) {
    tmplst = []
    for (var i=0; i<dat.length; i++) {
        if (dat[i].active) {
            tmplst.push(dat[i].profile.length);
        };
    };
    if (tmplst.length <= 1) {
        return 0;
    }
    else {
        return d3.max(tmplst);
    }
};

function getActiveData(dat) {
    tmplst = [];
    for (var i=0; i < dat.seq.length; i++) {
        if (dat.seq[i].active && data.seq[i].sequence != "") {
            tmplst.push([dat.seq[i].name, i]);
        }
    }
    return tmplst;
}

function getCorOptions(activeData, num) {
    datString = "";
    for (var i=0; i < activeData.length; i++) {
        if (activeData[i][1] != num) {
            datString += '<a class="dropdown-item" id="';
            datString += activeData[i][1];
            datString += '">'
            datString += activeData[i][1]+1 + ": " + activeData[i][0];
            datString += '</a>'
        }
    }
    return datString;
}

function trimString(str, len) {
    if (str.length <= len) {
        return str;
    }
    else {
        return str.substring(0, len-3) + "...";
    }
}

function parseFASTA(str) {
    re = /(>.*\n|\n)/g;
    return str.replace(re, '').toUpperCase()
}

function getPropertyString(id) {
    content = Number(id)+1 + ": " + data.seq[id].name + " (" + typeNames[data.seq[id].type] + ", " + scaleNames[data.seq[id].type][data.seq[id].scale] + ")";
    return trimString(content, 40)
}

function getRandomSequence() {
    var seq = "M";
    var letters = "ACDEFGHIKLMNPQRSTVWY"

    var max_len = Math.floor(Math.random() * 450) + 50;

    for (var i = 0; i < max_len; i++) {
        seq += letters.charAt(Math.floor(Math.random() * letters.length));
    }

    return seq;
}

function getRandomColor() {
    var color = "#";
    var letters = "0123456789ABCDEF";
    for (var i = 0; i < 6; i++) {
        color += letters.charAt(Math.floor(Math.random() * 16));
    }

    return color;
}

function decimalToHex(d) {
    hex = Number(d).toString(16).toUpperCase();

    while ( hex.length < 2) {
        hex = "0" + hex;
    }

    return hex;
}

function hslToRgb(h, s, l){
    var r, g, b;

    if(s == 0){
        r = g = b = l; // achromatic
    }else{
        var hue2rgb = function hue2rgb(p, q, t){
            if(t < 0) t += 1;
            if(t > 1) t -= 1;
            if(t < 1/6) return p + (q - p) * 6 * t;
            if(t < 1/2) return q;
            if(t < 2/3) return p + (q - p) * (2/3 - t) * 6;
            return p;
        }

        var q = l < 0.5 ? l * (1 + s) : l + s - l * s;
        var p = 2 * l - q;
        r = hue2rgb(p, q, h + 1/3);
        g = hue2rgb(p, q, h);
        b = hue2rgb(p, q, h - 1/3);
    }
    return "#" + decimalToHex(Math.round(r * 255)) + decimalToHex(Math.round(g * 255)) + decimalToHex(Math.round(b * 255));
}

function getBiasedColorHSL(i, S = null, L = null, noise = null) {
    H = ((i*3) % 10) / 10;
    if (noise) {
        if (! S) { S = Math.random(); }
        else { S = S + (Math.random() * noise) - (noise/2);
            if (S < 0) { S = 0;}
            if (S > 1) { S = 1;}
        }
        if (! L) { L = Math.random(); }
        else { L = L + (Math.random() * noise) - (noise/2);
            if (L < 0) { L = 0;}
            if (L > 1) { L = 1;}
        }
    }
    else {
        if (! S) {
            S = Math.random();
        }
        if (! L) {
            L = Math.random();
        }
    }

    return hslToRgb(H,S,L);
}

function getRandomColorHSL(H = null, S = null, L = null, noise = null) {
    if (noise) {

        if (! H) { H = Math.random(); }
        else {
            H = H + (Math.random() * noise) - (noise/2);
            if (H < 0) { H = 0;}
            if (H > 1) { H = 1;}
        }
        if (! S) { S = Math.random(); }
        else { S = S + (Math.random() * noise) - (noise/2);
            if (S < 0) { S = 0;}
            if (S > 1) { S = 1;}
        }
        if (! L) { L = Math.random(); }
        else { L = L + (Math.random() * noise) - (noise/2);
            if (L < 0) { L = 0;}
            if (L > 1) { L = 1;}
        }
    }
    else {
        if (! H) {
            H = Math.random();
        }
        if (! S) {
            S = Math.random();
        }
        if (! L) {
            L = Math.random();
        }
    }

    return hslToRgb(H,S,L);
}

function calculate_pearson(a, b) {
    overlap = calculateOverlap(data.seq[a].profile, data.seq[b].profile, data.seq[a].shift, data.seq[b].shift);
    if (overlap.p1.length == 0) {
        return NaN;
    }
    else {
        return spearson.round(spearson.correlation.pearson(yVal(overlap.p1), yVal(overlap.p2), true),2).toString();
    }
}


function calculate_determination(a, b) {
    overlap = calculateOverlap(data.seq[a].profile, data.seq[b].profile, data.seq[a].shift, data.seq[b].shift);
    if (overlap.p1.length == 0) {
        return NaN;
    }
    else {
        return spearson.round(Math.pow(spearson.correlation.pearson(yVal(overlap.p1), yVal(overlap.p2), true),2),2).toString();
    }
}

function calculate_RMSD(a, b) {
    overlap = calculateOverlap(data.seq[a].profile, data.seq[b].profile, data.seq[a].shift, data.seq[b].shift);
    return spearson.round(RMSD(yVal(overlap.p1), yVal(overlap.p2)),2).toString();
}

function calculate_spearman(a, b) {
    overlap = calculateOverlap(data.seq[a].profile, data.seq[b].profile, data.seq[a].shift, data.seq[b].shift);
    return spearson.round(spearson.correlation.spearman(yVal(overlap.p1), yVal(overlap.p2), true),2).toString();
}

function downloadLinkSVG(id) {
    // try not necessary since a modern browser can be assumed...
    try {
        var isFileSaverSupported = !!new Blob();
    } catch (e) {
        console.log("blob not supported", e)
    }

    var svg_obj = document.getElementById(id);
    var serializer = new XMLSerializer();
    var svg_string = serializer.serializeToString(svg_obj);
    var svg_blob = new Blob([svg_string], {type:"image/svg+xml;charset=utf-8"});
    saveAs(svg_blob, "myProfile.svg");
}

function downloadLinkPNG() {
    saveSvgAsPng(document.getElementsByTagName("svg")[0], "plot.png");
}

function setOptionsList() {
    active = getActiveData(data);
    html_content = "<select id='specific_selector' class='custom-select' multiple>";
    for (entry in active) {
        html_content += "<option value='" + active[entry][1] + "'>" + String(Number(active[entry][1])+1) + ": " + active[entry][0] + "</option>"
    }
    html_content += "</select>"
    $("#specific").html(html_content);
}

function readSingleFile(e) {
    var file = e.target.files[0];
    if (!file) {
        return;
    }
    var reader = new FileReader();
    reader.onload = function(e) {
        var contents = e.target.result;
        importer(contents);
    };
    reader.readAsText(file);
}

function importer(newData) {
    newData = JSON.parse(newData);
    if (assertCorrect(newData)) {
        data = newData;
        loadVizState();
    }
    else {
        alert("It seems like the file provided is not correctly formated.");
    }
}

function assertCorrect(d) {
    listOfKeys = ["active", "color", "name", "organism", "profile", "scale", "sequence", "shift", "smoothing_method", "thickness","min_existing", "type", "visible", "window"];
    if (d.hasOwnProperty('seq')) {
        for (i in d.seq) {
            for (key in listOfKeys) {
                if (!d.seq[i].hasOwnProperty(listOfKeys[key])) {
                    return false;
                }
            }
        }
        return true;
    }
    else {
        return false;
    }
}

function loadVizState() {
    // drawingProfiles(data);
    addInterfaceBasedOnData(data);
    // updateShifts(data);
    // updateInterface(data);
    // addLegend(num);
}

function createStringFromScaleValues(scaleValue, scaleType) {
    if (scaleType=="protein") {
        retStr = "";
        names = ["A", "R", "N", "D", "C", "Q", "E", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V"];
    }
    else if (scaleType=="rna") {
        retStr = "";
        names = ["A", "C", "G", "U"];
    }
    else if (scaleType=="dna") {
        retStr = "";
        names = ["A", "C", "G", "T"];
    }
    else {
        console.log("Type ", scaleType, " not supported");
        return;
    }

    var unavailable = false
    for (i in names) {
        var name = names[i];
        var value = scale[scaleType][scaleValue][name];

        if (isNaN(value)) {
            retStr += name + ":* ";
            unavailable = true;
        }
        else {
            retStr += name + ":" + scale[scaleType][scaleValue][name] + " "
        }
    }

    if (unavailable == true) {
        retStr += "<br>*not available  "
    }

    return retStr.slice(0,-1);

}

function createScaleList(targetID, selectedGroups, type) {
    if (type=="protein") {
        longOptions = { "p_alpha": "Alpha-propensity",
                        "p_beta": "Beta-propensity",
                        "p_charge": "Charge",
                        "p_comp": "Composition",
                        "p_hydro": "Hydrophobicity",
                        "p_pchem": "general Physico-Chemical",
                        "p_rnaAf": "RNA-Affinity",
                        "p_other": "Other"};

        newContent = "";
        for (var j=0; j < selectedGroups.length; j++) {
            var group = selectedGroups[j];
            newContent += "<option disabled>&#8195;"+longOptions[group]+"</option>"
            for (var i=0; i < classification[group].length; i++) {
                if (classification[group][i] in scale.protein) {
                    newContent += "<option>" + classification[group][i] + "</option>";
                }
            }
        }
        $("#"+targetID+".scaleList").html(newContent);
    }
    else if (type=="rna") {
        longOptions = { "r_comp" : "Composition",
                        "r_energy" : "Energy related",
                        "r_pchem" : "general Physico-Chemical"};

        newContent = "";
        for (var j=0; j < selectedGroups.length; j++) {
            var group = selectedGroups[j];
            newContent += "<option disabled>&#8195;"+longOptions[group]+"</option>"
            for (var i=0; i < classification[group].length; i++) {
                if (classification[group][i] in scale.rna) {
                    newContent += "<option>" + classification[group][i] + "</option>";
                }
            }
        }
        $("#"+targetID+".scaleList").html(newContent);
    }
    else if (type=="dna") {
        longOptions = { "d_comp" : "Composition"};

        newContent = "";
        for (var j=0; j < selectedGroups.length; j++) {
            var group = selectedGroups[j];
            newContent += "<option disabled>&#8195;"+longOptions[group]+"</option>"
            for (var i=0; i < classification[group].length; i++) {
                if (classification[group][i] in scale.dna) {
                    newContent += "<option>" + classification[group][i] + "</option>";
                }
            }
        }
        $("#"+targetID+".scaleList").html(newContent);
    }
}

function renderScaleDetails(targetID, selectedScale, type) {
    newContent = "";
    newContent += "<h6 style='font-size: 14px;'>"+scaleDescriptors[type][selectedScale].name+"</h6>"
    if (scaleDescriptors[type][selectedScale].aaIndexID) {
        newContent += "<a href='https://www.genome.jp/dbget-bin/www_bget?aaindex:"+
                    scaleDescriptors[type][selectedScale].aaIndexID+
                    "' target='_blank'>"+scaleDescriptors[type][selectedScale].aaIndexID+
                    "</a>";
    }
    newContent += "<br><br><span class='scaleView'>"+createStringFromScaleValues(selectedScale, type)+"</span><br><br>";
    newContent += scaleDescriptors[type][selectedScale].summary + "<br><br>";
    if (scaleDescriptors[type][selectedScale].author) {
        newContent += scaleDescriptors[type][selectedScale].author + "; " + scaleDescriptors[type][selectedScale].journal + "<br>";
    }
    if (scaleDescriptors[type][selectedScale].PMID) {
        newContent += "PMID: <a href='https://www.ncbi.nlm.nih.gov/pubmed/" + scaleDescriptors[type][selectedScale].PMID + "' target='_blank'>" + scaleDescriptors[type][selectedScale].PMID + "</a><br>";
    }
    $("#"+targetID+".scaleDescription").html(newContent);
}

function relativeUnitsCheck() {
    if (document.getElementById('relativeUnitsCheckbox').checked) {
        data.relative = true;
    } else {
        data.relative = false;
    }
    drawingProfiles(data);
}
function invertY1Check() {
    if (document.getElementById('invertY1Checkbox').checked) {
        data.y1_reverse = true;
    } else {
        data.y1_reverse = false;
    }
    drawingProfiles(data);
}
function invertY2Check() {
    if (document.getElementById('invertY2Checkbox').checked) {
        data.y2_reverse = true;
    } else {
        data.y2_reverse = false;
    }
    drawingProfiles(data);
}
function addScaleToName(id, name) {
    name = name.split(/@(.+)/)
    if (name.length == 1) {
        name = name[0].trim();
    } else {
        name = name[1].trim();  
    }
    s = data.seq[id].scale;
    return s + " @ " + name;
}
function removeScaleFromName(name) {
    name = name.split(/@(.+)/)
    if (name.length == 1) {
        name = name[0].trim();
    } else {
        name = name[1].trim();  
    }
    return name;
}