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

  PP.fetchPatients = function () {
      $.ajax({
          type: 'GET',
          url: '/my_patients_json',
          dataType: 'json',
          success: function (data) {
              var patients = data.result;
              for (var i = 0; i < patients.length; i++) {
                  $('<tr>').append(
                    $('<td>').html(PP.create_url('/individual', patients[i].individual)), 
                    $('<td>').text(patients[i].gender), 
                    $('<td>').html(PP.generatePhenotypesHtml(patients[i].phenotypes)),
                    $('<td>').text(+patients[i].phenotypeScore.toFixed(2)), 
                    $('<td>').text(+patients[i].hom_count.toFixed(2)), 
                    $('<td>').text(+patients[i].het_count.toFixed(2)),
                    $('<td>').html(PP.generateGeneHtml(patients[i].genes)) // candidate_genes
                  ).appendTo('#my_patients_table_body');
              }
              PP.showNumberOfPatients(patients.length);
              $('#my_patients #progress_row').remove();
              $('#my_patients_table').show();
              PP.initTableSorter('#my_patients_table');
              $('#my_patients_table').on('filterEnd', PP.showNumberOfPatients(patients.length));
          },
          error: function (error, msg) { PP.ajaxError(error, msg); }
      });
  };

  PP.generateIndexes = function(columns) {
    var idx             = {};
    idx.individual      = columns.indexOf('individual');
    idx.gender          = columns.indexOf('gender');
    idx.phenotypes      = columns.indexOf('phenotypes');
    idx.phenotypeScore  = columns.indexOf('phenotypeScore');
    idx.hom_count  = columns.indexOf('hom_count');
    idx.het_count  = columns.indexOf('het_count');
    idx.genes = columns.indexOf('genes');
    return (idx);
  };

  PP.generatePhenotypesHtml = function(features) {
    phenotypesHtml = '';
    for (var i = 0; i < features.length; i++) {
      phenotypesHtml = phenotypesHtml + ' ' + PP.create_url('/hpo', features[i].name, features[i].termId, 'chip');
    }
    return (phenotypesHtml);
  };

  PP.generateGeneHtml = function(genes) {
    var genesHtml = '';
    for (var i = 0; i < genes.length; i++) {
      genesHtml = genesHtml + ' ' + PP.create_url('/gene', genes[i], genes[i], 'chip');
    }
    return (genesHtml);
  };

  PP.showNumberOfPatients = function(patientsSize) {
    $('#total-patients').text(' (Total: ' + patientsSize + ')');
  };
}());
// End of PP module
