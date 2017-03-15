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

  // init the table directly under the header
  PP.initPatientFeatureTable = function (patient_id) {
    $.ajax({
      type: 'GET',
      url: 'individual_json/' + patient_id, // TODO - Parse error
      dataType: 'json',
      timeout: 120000,
      success: function(data) { PP.PatientFeatureSuccess(data); },
      error: function(data, msg) { PP.ajaxError(data, msg); }
    });
  };

  //
  PP.PatientFeatureSuccess = function(data) {
    PP.addPatientFeaturesInfo(data.result.features); // Or observed_features
    PP.addPatientConsanguinityInfo(data.result.family_history);
    PP.addPatientInheritanceModeInfo(data.reult); //todo
    PP.addPatientGenesInfo(data.result.genes);
    $('#info_table').show();
  };

  //
  PP.addPatientFeaturesInfo = function(features) {
    var features_html = '';
    for (var i = 0; i < features.length; i++) {
      features_html = features_html + ' ' + PP.create_url('/hpo', features[i].label, features[i].id, 'chip');
    }
    $('#patient_features').html(features_html);
  };

  //
  PP.addPatientConsanguinityInfo = function(family_history) {
    if (family_history === undefined ) return;
    if (family_history.consanguinity === undefined || family_history.consanguinity === null ) return;
    $('#patient_consanguinity').text(family_history.consanguinity);
  };

  //
  PP.addPatientInheritanceModeInfo = function() {
    $('#patient_inheritance_mode').text('TODO');
  };

  //
  PP.addPatientGenesInfo = function(genes) {
    var gene_html = '';
    for (var i = 0; i < genes.length; i++) {
      gene_html = gene_html + ' ' + PP.create_url('/gene', genes[i].gene, genes[i].gene, 'chip');
    }
    $('#patient_genes').html(gene_html);
  };


  PP.initOmimPlot = function () {

  };

  PP.initHomsTable = function(patient_id) {
    $.ajax({
      type: 'GET',
      url: 'homozygous_variants_json/' + patient_id,
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

  PP.initCompHetsTable = function(patient_id) {
    $.ajax({
      type: 'GET',
      url: 'compound_het_variants_json/' + patient_id,
      dataType: 'json',
      timeout: 120000,
      success: function(data) { PP.CompHetSuccess(data); },
      error: function(data, msg) { PP.ajaxError(data, msg); }
    });
  };

  PP.CompHetSuccess = function(data) {

  };

  PP.initVariantsTable = function(patient_id) {
    $.ajax({
      type: 'GET',
      url: 'rare_variants_json/' + patient_id,
      dataType: 'json',
      timeout: 120000,
      success: function(data) { PP.VariantsSuccess(data); },
      error: function(data, msg) { PP.ajaxError(data, msg); }
    });
  };

  PP.VariantsSuccess = function(data) {

  };

  PP.initExomiserTab = function(patient_id) {
    // $.ajax({
    //   type: 'GET',
    //   url: 'homozygous_variants_json/',
    //   dataType: 'json',
    //   timeout: 120000,
    //   success: function(data) { PP.HomsSuccess(data); },
    //   error: function(data, msg) { PP.ajaxError(data, msg); }
    // });

  };

}());
// End of PP module