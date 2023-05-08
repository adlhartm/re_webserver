var $buoop = {
    required:{e:-4,f:52,o:-3,s:10,c:64},
    // reminder: 1, reminderClosed: 24,
    // text: "This webtool requires modern browser technology not support by your browser, {brow_name}. Please update your webbrowser to a supported version: &lt;a{up_but}&gt;update&lt;/a&gt; or &lt;a{ignore_but}&gt;ignore&lt;/a&gt;.",
    jsshowurl: "https://browser-update.org/update.show.min.js",
    newwindow: true, noclose:true, no_permanent_hide: false,
    // test: true,
    api:2018.12

};
// var $buoop = {required:{e:-4,f:59,o:-3,s:-1,c:69},insecure:true,api:2018.12 }; 
function $buo_f(){   
    var e = document.createElement("script"); 
    e.src = "https://browser-update.org/update.min.js"; 
    document.body.appendChild(e);
};
try {document.addEventListener("DOMContentLoaded", $buo_f,false)}
catch(e){window.attachEvent("onload", $buo_f)}
