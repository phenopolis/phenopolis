/*
    PP - Phenopolis's JavaScript module

    Define a global PP (acronym for Phenopolis) object containing all
    PP associated methods:
*/

// define global PP object
var PP;
if (!PP) {
  PP = {};
}

// PP module
(function() {

  //
  PP.ajaxError = function(data, msg) {
    console.log(msg);
    console.log(data);
    // Cry & scream...
  };

  PP.initOmimPlot = function () {

  };

  PP.initHomsTable = function () {
    $.ajax({
      type: 'GET',
      url: '',
      dataType: 'json',
      timeout: 120000,
      success: function(data) { PP.HomsSuccess(data); },
      error: function(data, msg) { PP.ajaxError(data, msg); }
    });
  };

  PP.HomsSuccess = function(data) {
    var homs = data; // TODO
    for (var i = 0; i < homs.length; i++) {
      // generate tr
      $('<tr>').append(
         $('<td>').html(),
         $('<td>').text()
      ).appendTo('#homs_table_body');
    }
    $('#homs #progress_row').remove();
    $('#homs_table').show();
    PP.initTableSorter('#homs_table');
  };

  PP.initCompHetsTable = function () {

  };

  PP.initVariantsTable = function () {

  };

  PP.initExomiserTable = function () {

  };

}());
// End of PP module