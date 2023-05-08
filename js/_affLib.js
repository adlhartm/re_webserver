function yVal(arr) {
    // function to extract y-vals because i'm too stupid
    // to do this with a lambda, apparently
    n = arr.length;
    yVals = new Array(n);
    for (i=0;i<n;i++) {
        yVals[i] = arr[i].y;
    }

    return yVals;
}

function findOverlap(arr1, arr2) {
    // arr1 prot
    // arr2 rna
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

function windowAverage(arr, win, step) { // returns window avaraged list
    var output = []
    win = win*step;
    if (arr.length >= win) {
        tmp = arr.splice(0,win);
        output.push(d3.sum(tmp)/win)
    } else {
        return output;
    }
    while (arr.length >= step) {
        tmp.splice(0,step);
        tmp = tmp.concat(arr.splice(0,step));
        output.push(d3.sum(tmp)/win);
    }

    return output;
}

function savitzky_golay(arr, win, step) {
    // just a placeholder currently...
    return windowAverage(arr, win, step);
}
// function windowAverage(arr, win, step) { // returns window avaraged list

//     // somehow implement the codon steps correctly here...
//     console.log(arr.length);
//     if (step > 1) {
//         tmp = [];
//         for (i=0; i<(arr.length); i += step) {
//             tmp.push(d3.sum(arr.slice(i,i+step)));
//         }
//         arr = tmp;
//     }
//     console.log(arr.length);

//     win = 17;
//     norm = 323;
//     SG = [-21,-6,7,18,27,34,39,42,43,42,39,34,27,18,7,-6,-21];
    
//     // win = 7;
//     // norm = 21;
//     // norm = 28;
//     // norm = 42;
//     // SG = [-2,3,6,7,6,3,-2];
//     // SG = [-3,-2,-1,0,1,2,3];
//     // SG = [5,0,-3,-4,-3,0,5];
    
    
//     // norm = 408;
//     // norm = 3876;
    
    
//     // SG = [-8,-7,-6,-5,-4,-3,-2,-1,0,1,2,3,4,5,6,7,8];
//     // SG = [40,25,12,1,-8,-15,-20,-23,-24,-23,-20,-15,-8,1,12,25,40];

//     var output = []
//     if (arr.length >= win) {
//         tmp = arr.slice(0,win);
//         for (i = 0; i < win; i++) {
//             tmp[i] = (tmp[i] * SG[i])
//         }
//         output.push(d3.sum(tmp)/norm)
//     } else {
//         return output;
//     }
//     step = 1;
//     pos = step;
//     while (arr.length >= pos + win + step) {
//         tmp = arr.slice(pos,pos + win);
//         pos += step;
//         for (i = 0; i < win; i++) {
//             tmp[i] = (tmp[i] * SG[i])
//         }
//         output.push(d3.sum(tmp)/norm);
//     }

//     return output;
// }

// function protein_calculation(seq, protID, wndw) { // returns a window smoothed prot sequence
//     var seqArray = [];
//     var smoothArray = [];

//     for (var i = 0, len = seq.length; i < len; i++) {
//         seqArray.push(scale['protein'][protID][seq[i]]);
//     };

//     smoothArray = windowAverage(seqArray,wndw,1)

//     var xVals = [(wndw-1)/2], x=(wndw-1)/2;
//     while(x<(smoothArray.length+((wndw-1)/2) ) ){x+=1;xVals.push(x)}

//     var protData = []
//     for (var i = 0; i < smoothArray.length; i++) {
//         protData.push({"x":xVals[i], "y":smoothArray[i]});
//     };
//     return protData;
// }

// function rna_calculation(seq, rnaID, wndw) { // returns a window smoothed rna sequence
//     var seqArray = [];
//     var smoothArray = [];

//     for (var i = 0, len = seq.length; i < len; i++) {
//         seqArray.push(scale['rna'][rnaID][seq[i]]);
//     };

//     smoothArray = windowAverage(seqArray,wndw,3)

//     var xVals = [((wndw/3)-1)/2], x=((wndw/3)-1)/2;
//     while(x<(smoothArray.length+(((wndw/3) -1)/2) ) ){x+=1; xVals.push(x)}

//     var rnaData = []
//     for (var i = 0; i < smoothArray.length; i++) {
//         rnaData.push({"x":xVals[i], "y":smoothArray[i]})
//     };
//     return rnaData;
// }

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