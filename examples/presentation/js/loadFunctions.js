var window_p = $('.prot-window').each(function(element, value) {
    $(this).slider()
        .on("slide", function() {
            if (data.seq[this.id].type == "protein") {
                data.seq[this.id].window = Number(this.value);
                $("#"+this.id+".sliderWindow-value").text(this.value);
            } else {
                data.seq[this.id].window = 3*Number(this.value);
                $("#"+this.id+".sliderWindow-value").text(3*Number(this.value));
            }
            drawingProfiles(data);
        })
        .on("slideStop", function() {
            if (data.seq[this.id].type == "protein") {
                data.seq[this.id].window = Number(this.value);
                $("#"+this.id+".sliderWindow-value").text(this.value);
            } else {
                data.seq[this.id].window = 3*Number(this.value);
                $("#"+this.id+".sliderWindow-value").text(3*Number(this.value));
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
    data.seq[$(this).parent().attr('id')].type = this.text.toLowerCase();
    i = $(this).parent().attr('id');
    $("#" + i + ".btn-type").html(trimString(typeNames[data.seq[i].type], 9));
    updateInterface(data);
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
        if (data.seq[this.id].type == "rna") {
            data.seq[this.id].window = 3 * (curVal);
            $("#"+this.id+".firstRow-text").text("Window: " + 3*Number(curVal));
        }
        else {
            data.seq[this.id].window = curVal;
            $("#"+this.id+".firstRow-text").text("Window: " + curVal);
        }
        drawingProfiles(data);
    }
})
$(document).on('click', '.btn-less-wndw', function() {
    curVal = $("#"+this.id+".prot-window").data('slider').getValue();
    if (curVal > 1) { // magic number - CHANGE!!!
        curVal -= 1;
        $("#"+this.id+".prot-window").data('slider').setValue(curVal);
        if (data.seq[this.id].type == "rna") {
            data.seq[this.id].window = 3 * (curVal);
            $("#"+this.id+".firstRow-text").text("Window: " + 3*Number(curVal));
        }
        else {
            data.seq[this.id].window = curVal;
            $("#"+this.id+".firstRow-text").text("Window: " + curVal);
        }
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


// deal with clicks on submit button in database interface
$(document).on('click', '.btn-submit', function() {
    if ($('#'+this.id+'.btn-label').val() == "Uniprot") {
        download_data_Uniprot(this.id);
    }
    else if ($('#'+this.id+'.btn-label').val() == "ENA") {
        download_data_ENA(this.id);   
    }
    else {
        console.error("VALUE ERROR: database could not be determined "+$('#'+this.id+'.btn-label').val())
    }
})

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

// minimize/maximize interface elements upon licking the +/- toggles
$(".button").click(function(){
    if($(this).html() == "-"){
        $(this).html("+");
    }
    else{
        $(this).html("-");
    }
    $(this).parent().parent().find('.box').slideToggle();
});
$("#addInterfaceButton").tooltip();

$('#exportModal').on('shown.bs.modal', function(){
    setOptionsList();
});

$('.selectpicker').on('change', function(){
    createScaleList(this.id, $(this).val());
});

$('.scaleList').on('change', function(){
    renderScaleDetails(this.id, $(this).val());
});
/*
document.getElementById('importButton').onclick = function() {
    document.getElementById('dataImport').click();
}
document.getElementById('dataImport')
    .addEventListener('change', readSingleFile, false);
*/