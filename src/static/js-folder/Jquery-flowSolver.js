var imageIsSup = false

 //on change hide all divs linked to select and show only linked to selected option
 $("#flowEquation").change(function () {
   //Saves in a variable the wanted div
   if ($("#flowEquation option:selected").val() == "Brinkman") {
     $("#rockMatrixPerm").collapse("show");
     $("#extraMeshRefinement").collapse("hide");
   } 
   else if($("#flowEquation option:selected").val() == "Stokes") {
    $("#rockMatrixPerm").collapse("hide");
    $("#extraMeshRefinement").collapse("show");
   }
   else {
     $("#rockMatrixPerm").collapse("hide");
     $("#extraMeshRefinement").collapse("hide");
   }
 });

  //on change hide all divs linked to select and show only linked to selected option
  $("#flowEquation").change(function () {
    //Saves in a variable the wanted div
    if ($("#flowEquation option:selected").val() == "Brinkman") {
      $("#rockMatrixPerm").collapse("show");
    } else {
      $("#rockMatrixPerm").collapse("hide");
    }
  });
 

 $("#macro-form").change(function () {
  showSimulateButton();
});

function showSimulateButton() {
  var flowEquation = $("#flowEquation option:selected").val();
  var radio = $("input[type=radio]:checked").length;
  var inpressure =$("#inletPressure").val();
  var outpressure =$("#outletPressure").val();
  var rockmatrix = $('#rockMatrixPerm-input').val()

  var inputImage = document.getElementById("file-picker").files.length;
  var btn = document.getElementById("convert-btn-div");

  if (flowEquation == "Stokes" && radio == 2 && outpressure != "" && inpressure != "" && inputImage == 1 ) {
    btn.style = "display:block";
  } else if (flowEquation == "Brinkman" && radio == 2 && outpressure != "" && inpressure != ""  && rockmatrix != ""  && inputImage == 1) {
    btn.style = "display:block";
  } else {
    btn.style = "display:none";
  }
}

$("#file-picker").change(function () {
  $("#img-true-size").removeAttr('src')
  $("#alert-msg")[0].style = "display:none;";
  checkEXT(this);

  if (imageIsSup == true){
  if (this.files[0].size > 1000000){
    $("#alert-msg").text("The image file is too big. Please upload an binary image file under 1 MB.");
    $("#alert-msg")[0].style = "display:block";
    imageIsSup = false
  }
}
  if (imageIsSup == true){
    readURL(this);
 }
}
);

 function checkEXT(input) {
  var ext = input.files[0].name.substr(input.files[0].name.lastIndexOf(".") + 1).toLowerCase();
  var erroMsg = $("#alert-msg");

  if (ext == "jpg" || ext == "png" || ext == "jpeg") {
  erroMsg.text("");
    imageIsSup = true;
  } else {

    $("#img-true-size").removeAttr('src');
    erroMsg.text("This file extension is no supported. JPG, JPEG and PNG are supported.");
    erroMsg[0].style = "display:block";
    imageIsSup = false;
  }
 }

 function checkBinary(){

  var img = $("#img-true-size")
  debug1 = img

  const canvas = document.createElement('canvas');
  const ctx = canvas.getContext('2d');  
  var cvs = document.createElement('canvas');

  cvs.width = img[0].naturalWidth;
  cvs.height = img[0].naturalHeight; // Display image in the client

  ctx.drawImage(img[0],0,0,cvs.width,cvs.height);

  var idt = ctx.getImageData(0,0,cvs.width,cvs.height);

  var d = idt.data
  var counter = 0
  var i

  for (i=0; i<d.length; i++){
    if(d[i] != 0 && d[i] != 255) {
      counter++}}

  if (counter != 0){
    // var warning_text = "The image file is not binary. Upload a binary image or convert your image" + <a href="#" class="alert-link">here</a> + '.'
    var warning_text = "The image file is not binary. Please upload a binary image."
    $("#alert-msg").text(warning_text);
    $("#alert-msg")[0].style = "display:block";
    imageIsSup = false;
    img.removeAttr('src')
  //  $("#img-converted-size").removeAttr('src')
  }
  else {
    img[0].width = 400
  }
}

  // Display image in the client
 function readURL(input) {
  if (input.files && input.files[0]) {
    debug2 = input
    var reader = new FileReader();

    reader.onload = function (e) {
     //  $("#img-converted-size").attr("src", e.target.result);
      $("#img-true-size").attr("src", e.target.result);
    };

    reader.readAsDataURL(input.files[0]); // convert to base64 string
  }
}