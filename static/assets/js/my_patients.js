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
      timeout: 120000,
      success: function(data) {
        var patients = data.result;
        for (var i = 0; i < patients.length; i++) {
          $('<tr>').append(
            $('<td>').html(PP.create_url('/individual', patients[i].external_id)),
            $('<td>').text(patients[i].sex),
            $('<td>').html(PP.generateFeaturesHtml(patients[i].features)),
            $('<td>').text(patients[i].specificity.score.toFixed(2)),
            $('<td>').text(patients[i].rare_homozygous_variants_count),
            $('<td>').text(patients[i].rare_compound_hets_count),
            $('<td>').text(patients[i].rare_variants_count),
            $('<td>').html(PP.generateGeneHtml(patients[i].genes))
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

  PP.generateFeaturesHtml = function(features) {
    var features_html = '';
    for (var i = 0; i < features.length; i++) {
      features_html = features_html + ' ' + PP.create_url('/hpo', features[i].label, features[i].id, 'chip');
    }
    return (features_html);
  };

  PP.generateGeneHtml = function(genes) {
    var genes_html = '';
    for (var i = 0; i < genes.length; i++) {
      genes_html = genes_html + ' ' + PP.create_url('/gene', genes[i], genes[i], 'chip');
    }
    return (genes_html);
  };

  PP.showNumberOfPatients = function(patientsSize) {
    $('#total-patients').text(' (Total: ' + patientsSize + ')');
  };
}());
// End of PP module