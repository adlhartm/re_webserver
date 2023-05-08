var window_p = $('.prot-window').each(function(element, value) {
    $(this).slider()
        .on("slide", function() {
            if (data.seq[this.id].type == "protein") {
                data.seq[this.id].window = Number(this.value);
                $("#"+this.id+".sliderWindow-value").text(this.value);
            } else {
                data.seq[this.id].window = Number(this.value);
                $("#"+this.id+".sliderWindow-value").text(Number(this.value));
            }
            drawingProfiles(data);
        })
        .on("slideStop", function() {
            if (data.seq[this.id].type == "protein") {
                data.seq[this.id].window = Number(this.value);
                $("#"+this.id+".sliderWindow-value").text(this.value);
            } else {
                data.seq[this.id].window = Number(this.value);
                $("#"+this.id+".sliderWindow-value").text(Number(this.value));
            }
            drawingProfiles(data);
        })
        .data('slider');
});
var shift_p = $('.shift').each(function(element, value) {
    $(this).slider()
        .on("slide", function() {
            data.seq[this.id].shift = Number(this.value);
            drawingProfiles(data);
        })
        .on("slideStop", function() {
            data.seq[this.id].shift = Number(this.value);
            drawingProfiles(data);
        })
        .data('slider');
});

updateShifts(data);

// insert the sequence into the text field
$('.form-inline').each(function(element, value) {
    $(this).val(data.seq[this.id].sequence);
});
// redraw the graph upon changes in the sequence text field
$('.form-inline').each(function(element, value) {
    $(this).bind('input propertychange', function() {
        data.seq[this.id].sequence = this.value;
        drawingProfiles(data);
        updateShifts(data);
    });    
});


// editing the name of sequences
$('.name').each(function(element, value) {
    $(this).on("change", function() {
        data.seq[this.id].name = this.value;
        updateLegend(this.id);
    })
});


// seleceting scales in the main interface
$(document).on('click', '.scale-selector a', function() {
    data.seq[$(this).parent().attr('id')].scale = this.id;
    cur_id =  $(this).parent().attr('id');
    $("#" + cur_id + ".btn-scale").html(trimString(scaleNames[data.seq[cur_id].type][data.seq[cur_id].scale], 9));
    drawingProfiles(data);
    updateLegend(cur_id);
})


// selecting different types of polymers
$(document).on('click', '.type-selector a', function() {
    i = $(this).parent().attr('id');
    data.seq[i].type = this.text.toLowerCase();
    $("#" + i + ".btn-type").html(trimString(typeNames[data.seq[i].type], 9));
    if (this.text.toLowerCase()=="protein") {
        $("#"+i+".selectpicker").html(`
                                <option value="p_alpha">Alpha-propensity</option>
                                <option value="p_beta">Beta-propensity</option>
                                <option value="p_charge">Charge</option>
                                <option value="p_comp">Composition</option>
                                <option value="p_hydro">Hydrophobicity</option>
                                <option value="p_pchem">Physcial-Chemical</option>
                                <option value="p_rnaAf">RNA-Affinity</option>
                                <option value="p_other">Other</option>
                                `)
        $("#"+i+".selectpicker").selectpicker('refresh')
    }
    else if (this.text.toLowerCase()=="rna") {
        $("#"+i+".selectpicker").html(`
                                <option value="r_comp">Composition</option>
                                `)
        $("#"+i+".selectpicker").selectpicker('refresh')
    }
    drawingProfiles(data);
})

// selecting the correlation partner
$(document).on('click', '.corr-selector a', function() {
    corList[$(this).parent().attr('id')] = this.id;
    updateCorResults(corList);
})

// click on shift + button
$(document).on('click', '.btn-more-wndw', function() {
    curVal = $("#"+this.id+".prot-window").data('slider').getValue();
    if (curVal < 63) { // magic number - CHANGE!!!
        curVal += 1;
        $("#"+this.id+".prot-window").data('slider').setValue(curVal);
        // if (data.seq[this.id].type == "rna") {
        //     data.seq[this.id].window = 3 * (curVal);
        //     $("#"+this.id+".firstRow-text").text("Window: " + 3*Number(curVal));
        // }
        // else {
        data.seq[this.id].window = curVal;
        $("#"+this.id+".firstRow-text").text("Window: " + curVal);
        // }
        drawingProfiles(data);
    }
})
$(document).on('click', '.btn-less-wndw', function() {
    curVal = $("#"+this.id+".prot-window").data('slider').getValue();
    if (curVal > 1) { // magic number - CHANGE!!!
        curVal -= 1;
        $("#"+this.id+".prot-window").data('slider').setValue(curVal);
        // if (data.seq[this.id].type == "rna") {
        //     data.seq[this.id].window = 3 * (curVal);
        //     $("#"+this.id+".firstRow-text").text("Window: " + 3*Number(curVal));
        // }
        // else {
            data.seq[this.id].window = curVal;
            $("#"+this.id+".firstRow-text").text("Window: " + curVal);
        // }
        drawingProfiles(data);
    }
})
$(document).on('click', '.btn-more', function() {
    if (data.seq[this.id].shift < findMaxLength(data.seq)-3) {
        data.seq[this.id].shift += 1;
        $("#"+this.id+".shift").data('slider').setValue(data.seq[this.id].shift);
        $("#"+this.id+".thirdRow-text").text("Shift: " + data.seq[this.id].shift);
        drawingProfiles(data);
    }
})
// click on shift - button
$(document).on('click', '.btn-less', function() {
    if (data.seq[this.id].shift > 0) {
        data.seq[this.id].shift -= 1;
        $("#"+this.id+".shift").data('slider').setValue(data.seq[this.id].shift);
        $("#"+this.id+".thirdRow-text").text("Shift: " + data.seq[this.id].shift);
        drawingProfiles(data);
    }
})
// click on thickness - button
$(document).on('click', '.btn-more-thick', function() {
    if (data.seq[this.id].thickness < 10) {
        data.seq[this.id].thickness *= 10;
        data.seq[this.id].thickness += 1;
        data.seq[this.id].thickness /= 10;
        $("#"+this.id+".thickness-slider").data('slider').setValue(data.seq[this.id].thickness);
        $("#"+this.id+".thickness-text").text("Thickness: " + data.seq[this.id].thickness);
        drawingProfiles(data);
    }
})
// click on thickness - button
$(document).on('click', '.btn-less-thick', function() {
    if (data.seq[this.id].thickness > 0.1) {
        data.seq[this.id].thickness *= 10;
        data.seq[this.id].thickness -= 1;
        data.seq[this.id].thickness /= 10;
        $("#"+this.id+".thickness-slider").data('slider').setValue(data.seq[this.id].thickness);
        $("#"+this.id+".thickness-text").text("Thickness: " + data.seq[this.id].thickness);
        drawingProfiles(data);
    }
})


// deal with clicks on submit button in database interface
$(document).on('click', '.btn-submit', function() {
    if ($('#'+this.id+'.btn-label').val() == "Uniprot") {
        download_data_Uniprot(this.id);
    }
    else if ($('#'+this.id+'.btn-label').val() == "ENA") {
        download_data_ENA(this.id);   
    }
    else {
        $("#setup" + this.id + " .alert_box").empty();
        $("#setup" + this.id + " .alert_box").append("<div class='alert alert-warning alert-dismissable'><button type='button' class='close' data-dismiss='alert' aria-hidden='true'>&times;</button><strong>Warning</strong> Please select a database first</div>");
    }
})

$(document).on('keyup', '.sequence-id-input', function (e) {
    if (e.keyCode == 13) {
        if ($('#'+this.id+'.btn-label').val() == "Uniprot") {
            download_data_Uniprot(this.id);
        }
        else if ($('#'+this.id+'.btn-label').val() == "ENA") {
            download_data_ENA(this.id);   
        }
        else {
            $("#setup" + this.id + " .alert_box").empty();
            $("#setup" + this.id + " .alert_box").append("<div class='alert alert-warning alert-dismissable'><button type='button' class='close' data-dismiss='alert' aria-hidden='true'>&times;</button><strong>Warning</strong> Please select a database first</div>");
        }
    }
});

$(document).on('click', '.btn-clr', function() {
    $('#'+this.id+'.form-inline').val("");
    data.seq[this.id].sequence = "";
    drawingProfiles(data);
})


// switch label and value of dropdown database selector
$(document).on('click', '.database-selector a', function() {
    $(this).parents('.database-menu').find('.btn-label').html($(this).text());
    $(this).parents('.database-menu').find('.btn-label').val($(this).text());
})


// deal with color changes in the colorpicker
function updateColor(jscolor, id) {
    data.seq[id].color = "#" + jscolor
    $("#legend"+id+" .color-display").css('background-color', data.seq[id].color)
    drawingProfiles(data);
}

$(document).on('click', '.close-button', function() {
    $('.tooltip').tooltip('hide');
    delete corList[this.id];
    // delete all references to this.id in corList
    for(var d in corList) {
        if(corList.hasOwnProperty(d) && corList[d] == this.id) {
            delete corList[d];
            $("#" + d + ".corr_output").text("");
        }
    }
    
    removeInterfaceBox(this.id);
    $("#" + this.id + ".corr_output").text("");

    if (getActiveData(data).length == 0) {
        $("#addInterfaceButton").text("Run Analysis")
    }
    
})

$("#addInterfaceButton").tooltip();

$('#exportModal').on('shown.bs.modal', function(){
    setOptionsList();
});

$('.selectpicker').on('change', function(){
    createScaleList(this.id, $(this).val(), data.seq[this.id].type);
});

$('.scaleList').on('change', function(){
    renderScaleDetails(this.id, $(this).val(), data.seq[this.id].type);
});
/*
document.getElementById('importButton').onclick = function() {
    document.getElementById('dataImport').click();
}
document.getElementById('dataImport')
    .addEventListener('change', readSingleFile, false);
*/