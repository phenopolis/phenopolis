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

  PP.fetchPatients = function() {
    $.ajax({
      type: 'GET',
      url: '/my_patients_json',
      dataType: 'json',
      success: function(data) {
        var patients = data.result.data;
        var columns  = data.result.columns;
        var idx      = PP.generateIndexes(columns);
        for (var i = 0; i < patients.length; i++) {
          $('<tr>').append(
            $('<td>').html(PP.create_url('/individual', patients[i][ idx.individual ] )), // individual
            $('<td>').text( patients[i][ idx.gender ] ), // gender
            $('<td>').html(PP.generatePhenotypesHtml( patients[i][ idx.phenotypes ] )), // phenotypes
            $('<td>').text( +patients[i][ idx.phenotypeScore ].toFixed(2) ), // phenotypeScore
            $('<td>').text( +patients[i][ idx.hom_count ].toFixed(2) ), // hom_count
            $('<td>').text( +patients[i][ idx.het_count ].toFixed(2) ), // het_count
            $('<td>').html(PP.generateGeneHtml(patients[i][idx.genes])) // candidate_genes
          ).appendTo('#my_patients_table_body');
        }
        PP.showNumberOfPatients(patients.length);
        $('#my_patients #progress_row').remove();
        $('#my_patients_table').show();
        PP.initTableSorter('#my_patients_table');
        $('#my_patients_table').on('filterEnd', PP.showNumberOfPatients(patients.length) );
      },
      error: function(error, msg) { PP.ajaxError(error, msg); }
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
      phenotypesHtml = phenotypesHtml + ' ' + PP.create_url('/hpo', features[i].data.name, features[i].data.termId, 'chip');
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
