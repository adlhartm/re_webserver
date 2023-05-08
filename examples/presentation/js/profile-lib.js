function findMin(dat, type) {
    tmplst = []
    for (var i=0; i<dat.length; i++) {
        if (dat[i].type == type & dat[i].active) {
            tmplst.push(d3.min(dat[i].profile, function(d) {return d.y}));
        };
    };
    return d3.min(tmplst);
};

function findMax(dat, type) {
    tmplst = []
    for (var i=0; i<dat.length; i++) {
        if (dat[i].type == type & dat[i].active) {
            tmplst.push(d3.max(dat[i].profile, function(d) {return d.y}));
        };
    };
    return d3.max(tmplst);
};

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
        if (dat.seq[i].active) {
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
    return str.replace(re, '')
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

    // return hsl(' + H + ', ' + S + '%, ' + L + '%)';
}

function calculate_correlation(a, b) {
    profile_1 = data.seq[a].profile;
    profile_2 = data.seq[b].profile;
    limits = findOverlap(profile_1, profile_2);
    pSlice = data.seq[a].profile.slice(limits[0], limits[1]);
    rSlice = data.seq[b].profile.slice(limits[2], limits[3]);
    return spearson.round(spearson.correlation.pearson(yVal(pSlice), yVal(rSlice), true),2).toString();
}

function downloadLinkSVG(canvas) {
    try {
        // Blob is not supported in the following browsers:
        // 
        // Firefox <20
        // Opera <15
        // Safari <6
        // 
        // This could be fixed by using https://github.com/eligrey/Blob.js/blob/master/Blob.js
        // 
        var isFileSaverSupported = !!new Blob();
    } catch (e) {
        console.log("blob not supported", e)
    }

    var html = d3.select(canvas)
        .attr("title", "pro.viz.export")
        .attr("version", 1.1)
        .attr("xmlns", "http://www.w3.org/2000/svg")
        .node().parentNode.innerHTML;
    
    var blob = new Blob([html], {type: "image/svg+xml"});
    saveAs(blob, "myProfile.svg");
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
    listOfKeys = ["active", "color", "name", "organism", "profile", "scale", "sequence", "shift", "smoothing_method", "thickness", "type", "visible", "window"];
    if (d.hasOwnProperty('seq')) {
        for (i in d.seq) {
            for (key in listOfKeys) {
                if (!d.seq[i].hasOwnProperty(listOfKeys[key])) {
                    // console.log(i, d.seq[i], key, listOfKeys[key]);
                    return false;
                }
            }
        }
        return true;
    }
    else {
        // console.log(d);
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

    for (name in names) {
        retStr += names[name] + ":" + scale[scaleType][scaleValue][names[name]] + " ";
    }
    return retStr.slice(0,-1);

}

function createScaleList(targetID, selectedGroups) {
    newContent = "";
    for (var j=0; j < selectedGroups.length; j++) {
        var group = selectedGroups[j];
        for (var i=0; i < classification[group].length; i++) {
            if (classification[group][i] in scale.protein) {
                newContent += "<option>" + classification[group][i] + "</option>";
            }
        }
    }
    $("#"+targetID+".scaleList").html(newContent);
}

function renderScaleDetails(targetID, selectedScale) {
    // <h6 style="font-size: 14px;">Signal sequence helical potential</h6></br>
    // <a href="https://www.genome.jp/dbget-bin/www_bget?aaindex:ARGP820102">ARGP820102</a><br><br>

    // <span style="font-family: 'Source Code Pro', monospace;">
    // A:1.18 R:0.20 N:0.23 D:0.05 C:1.89 Q:0.72 E:0.11 G:0.49 H:0.31 I:1.45 L:3.23 K:0.06 M:2.67 F:1.96 P:0.76 S:0.97 T:0.84 W:0.77 Y:0.39 V:1.08</span><br><br>

    // Structural prediction of membrane-bound proteins<br><br>
    // Argos, P., Rao, J.K.M. and Hargrave, P.A.; Eur. J. Biochem. 128, 565-575 (1982) <br>
    // PMID: <a href="https://www.ncbi.nlm.nih.gov/pubmed/7151796">7151796</a><br>
    newContent = "";
    newContent += "<h6 style='font-size: 14px;'>"+scaleDescriptors['protein'][selectedScale].name+"</h6>"
    if (scaleDescriptors['protein'][selectedScale].aaIndexID) {
        newContent += "<a href='https://www.genome.jp/dbget-bin/www_bget?aaindex:"+
                    scaleDescriptors['protein'][selectedScale].aaIndexID+
                    "'>"+scaleDescriptors['protein'][selectedScale].aaIndexID+
                    "</a>";
    }
    newContent += "<br><br><span class='scaleView'>"+createStringFromScaleValues(selectedScale, "protein")+"</span><br><br>";
    newContent += scaleDescriptors['protein'][selectedScale].summary + "<br><br>";
    newContent += scaleDescriptors['protein'][selectedScale].author + "; " + scaleDescriptors['protein'][selectedScale].journal + "<br>";
    if (scaleDescriptors['protein'][selectedScale].PMID) {
        newContent += "PMID: <a href='https://www.ncbi.nlm.nih.gov/pubmed/" + scaleDescriptors['protein'][selectedScale].PMID + "'>" + scaleDescriptors['protein'][selectedScale].PMID + "</a><br>";
    }
    $("#"+targetID+".scaleDescription").html(newContent);
}