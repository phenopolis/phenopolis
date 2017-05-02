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
  PP.initPatientFeatureTable = function(patient_id) {
    $.ajax({
      type: 'GET',
      url: '/individual_json/' + patient_id,
      dataType: 'json',
      timeout: 120000,
      success: function(data) { PP.PatientFeatureSuccess(data); },
      error: function(data, msg) { PP.ajaxError(data, msg); }
    });
  };

  //
  PP.PatientFeatureSuccess = function(data) {
    PP.addPatientGenderInfo(data.result.gender); // Or observed_features
    PP.addPatientFeaturesInfo(data.result.observed_features); // Or observed_features
    PP.addPatientConsanguinityInfo(data.result.family_history);
    PP.addPatientGenesInfo(data.result.genes);
    PP.submitEditedIndividual();
    $('.info_table').show();
    $('#edit_icon').show();
    $('.modal').modal();
  };

  PP.addPatientGenderInfo = function(gender) {
    var gender_full;
    if (gender === 'M') {
      gender_full = 'Male';
    } else if (gender === 'F') {
      gender_full = 'Female';
    } else {
      gender_full = 'Unknown';
    }
    $('#'+gender_full+'_edit').prop('checked', true)
    $('#patient_gender').text(gender_full);
  };

  //
  PP.addPatientConsanguinityInfo = function(family_history) {
    if (family_history !== undefined) {
      $('#consanguinity_unknown_edit').prop('checked', true);
      return $('#patient_consanguinity').closest('tr').hide();
    }
    if (family_history.consanguinity === undefined || family_history.consanguinity === null) {
      $('#consanguinity_unknown_edit').prop('checked', true);
      return $('#patient_consanguinity').closest('tr').hide();
    }
    if (family_history.consanguinity == 'Y') {
      $('#consanguinity_yes_edit').prop('checked', true);
    } else if (family_history.consanguinity == 'N') {
      $('#consanguinity_no_edit').prop('checked', true);
    }
    $('#patient_consanguinity').text(family_history.consanguinity);
  };

  //
  PP.addPatientFeaturesInfo = function(features) {
    var features_html = '';
    var features_array = [];
    for (var i = 0; i < features.length; i++) {
      features_html = features_html + ' ' + PP.create_url('/hpo', features[i].label, features[i].id, 'chip');
      features_array.push( {tag: features[i].label} );
    }
    $('#patient_features').html(features_html);
    PP.setupChipsAutoComplete('#features_edit', features_array, '/phenotype_suggestions')
  };

  //
  PP.typeAheadSelectBind = function (wrapper, endpoint) {
    $(wrapper).bind('typeahead:select', '.chips-initial input', function(ev, suggestion) {
      currentData = $(wrapper+' .chips-autocomplete').material_chip('data');
      currentData.push( { tag: suggestion })
      $(wrapper+' .chips input').val('')
      // Delete current Material_chip binding and create new one with the updated data... 
      $(wrapper+' .chips-initial').material_chip('destroy')
      $(wrapper+' .chips-initial').material_chip({data: currentData});
      var inputAutoComplete = PP.initialiseBloodHound(endpoint);
      PP.initialiseTypeahead(wrapper+' .chips input', inputAutoComplete);
      PP.setupDropdownStyling(wrapper);
    });
  }

  // Makes the dropdown full width of the input
  PP.setupDropdownStyling = function(wrapper) {
    $(wrapper).bind('typeahead:active', '.chips-initial input', function() {
      var position = $(wrapper + ' .twitter-typeahead').position().left - 24;
      var width = $(wrapper + ' .chips-initial').width()
      $(wrapper + ' .autocomplete-content').css('margin-left', '-'+(position) +'px')
      $(wrapper + ' .autocomplete-content').css('width', (width) +'px');
    });
  }

  PP.focusOnInputAfterSelect = function(wrapper) {
    $(wrapper).bind('typeahead:select', '.chips-initial input', function() {
      $(wrapper + ' .chips-initial input').focus();
    });
  }

  //
  PP.setupChipsAutoComplete = function(wrapper, data, endpoint) {
    $(wrapper+' .chips-initial').material_chip({data: data})
    var inputAutoComplete = PP.initialiseBloodHound(endpoint);
    PP.initialiseTypeahead(wrapper+' .chips-initial input', inputAutoComplete);
    PP.typeAheadSelectBind(wrapper, endpoint);
    PP.setupDropdownStyling(wrapper);
    PP.focusOnInputAfterSelect(wrapper);
  }

  //
  PP.addPatientGenesInfo = function(genes) {
    var gene_html = '';
    var genes_array = [];
    for (var i = 0; i < genes.length; i++) {
      gene_html = gene_html + ' ' + PP.create_url('/gene', genes[i].gene, genes[i].gene, 'chip');
      genes_array.push( {tag: genes[i].gene} );
    }
    $('#patient_genes').html(gene_html);
    PP.setupChipsAutoComplete('#genes_edit', genes_array, '/gene_suggestions')
  };

  PP.submitEditedIndividual = function() {
    $('#save_modal').on('click', function(){
      $('#confirm_edit_modal').modal({ endingTop: '20%' });
      $('#confirm_edit_modal').modal('open');
    });
    $('#submit_edit_modal').on('click', function() {
      $.ajax({
        type: 'POST',
        url: '/update_patient_data',
        data: $('#edit_patient_data_form').serialize(),
        dataType: 'json',
        timeout: 120000,
        success: function(data) {
          // update page 
        },
        error: function(data, msg) {
          // 
        }
      });
    })
  }

  //
  PP.initOmimPlot = function(patient_id) {
    $.ajax({
      type: 'GET',
      url: '/venn_json/' + patient_id,
      dataType: 'json',
      timeout: 120000,
      success: function(data) { PP.vennSuccess(data); },
      error: function(data, msg) { PP.ajaxError(data, msg); }
    });
  };

  //
  PP.vennSuccess = function(data) {
    sets = PP.generateSets(data.result);
    PP.createVennDiagram(sets);
    $('#omim_tab #progress_row').remove();
  };

  //
  PP.generateSets = function(data) {
    sets = [];
    for (var i = 0; i < data.length; i++) {
      sets.push({ "sets": data[i].key, "size": data[i].value.length });
    }
    return sets;
  };

  PP.createVennDiagram = function(sets) {
    // draw venn diagram
    var div = d3.select("#venn");
    div.datum(sets).call(venn.VennDiagram().width(500).height(500));

    div.selectAll("path")
        .style("stroke-opacity", 0)
        .style("stroke", "#fff")
        .style("stroke-width", 3);

    // add a tooltip
    var tooltip = d3.select("body").append("div")
      .attr("class", "venntooltip");

    // add listeners to all the groups to display tooltip on mouseover
    div.selectAll("g")
      .on("mouseover", function(d, i) {
        // sort all the areas relative to the current item
        venn.sortAreas(div, d);

        // Display a tooltip with the current size
        tooltip.transition().duration(400).style("opacity", 0.9);
        tooltip.text(d.size + " Genes");

        // highlight the current path
        var selection = d3.select(this).transition("tooltip").duration(400);
        selection.select("path")
          .style("stroke-width", 3)
          .style("fill-opacity", d.sets.length == 1 ? 0.5 : 0.1)
          .style("stroke-opacity", 1);
      })
      .on("mousemove", function() {
        tooltip.style("left", (d3.event.pageX) + "px")
          .style("top", (d3.event.pageY - 28) + "px");
      })
      .on("mouseout", function(d, i) {
        tooltip.transition().duration(400).style("opacity", 0);
        var selection = d3.select(this).transition("tooltip").duration(400);
        selection.select("path")
          .style("stroke-width", 0)
          .style("fill-opacity", d.sets.length == 1 ? 0.25 : 0.0)
          .style("stroke-opacity", 0);
      });
  };

  PP.initHomsTable = function(patient_id) {
    $.ajax({
      type: 'GET',
      url: '/homozygous_variants_json/' + patient_id,
      dataType: 'json',
      timeout: 120000,
      success: function(data) { PP.variantTableSuccess(data); },
      error: function(data, msg) { PP.ajaxError(data, msg); }
    });
  };

  PP.initCompHetsTable = function(patient_id) {
    // $.ajax({
    //   type: 'GET',
    //   url: '/compound_het_variants_json/' + patient_id,
    //   dataType: 'json',
    //   timeout: 120000,
    //   success: function(data) { PP.variantTableSuccess(data); },
    //   error: function(data, msg) { PP.ajaxError(data, msg); }
    // });
  };

  PP.initVariantsTable = function(patient_id) {
    // $.ajax({
    //   type: 'GET',
    //   url: '/rare_variants_json/' + patient_id,
    //   dataType: 'json',
    //   timeout: 120000,
    //   success: function(data) { PP.variantTableSuccess(data); },
    //   error: function(data, msg) { PP.ajaxError(data, msg); }
    // });
  };

  PP.variantTableSuccess = function(data) {
    var d = data.result;
    for (var i = 0; i < d.length; i++) {
      // generate tr
      $('<tr>').append(
        $('<td>').html(PP.generateGeneNames(d[i].transcript_consequences)),
        $('<td>').text(' '),
        $('<td>').text(' '),
        $('<td>').html(PP.generateHpoLinks(d[i].HPO)),
        $('<td>').text(' '),
        $('<td>').html(PP.generateVariantHtml(d[i].synonym_variant_id)),
        $('<td>').text(' '),
        $('<td>').text(' '),
        $('<td>').text(' '),
        $('<td>').text(d[i].most_severe_consequence),
        $('<td>').html(PP.createHgvspTd(d[i].canonical_hgvsp, d[i].canonical_hgvsc)),
        $('<td>').text(' '),
        $('<td>').text(' '),
        $('<td>').html(PP.generateIndividualHtml(d[i].hom_samples)),
        $('<td>').text()
      ).appendTo('#homs_table_body');
    }
    $('#homs_tab #progress_row').remove();
    $('#homs_table').show();
    PP.initTableSorter('#homs_table');
  };

  PP.generateGeneNames = function(transcript_consequences) {
    geneNames = [];
    for (var i = 0; i < transcript_consequences.length; i++) {
      url = PP.create_url('/gene', transcript_consequences[i].gene_symbol);
      if (geneNames.indexOf(url) === -1) {
        geneNames.push(url);
      }
    }
    return (geneNames.join('<br>'));
  };

  PP.generateHpoLinks = function(hpoTerms) {
    hpoHtml = '';
    for (var i = 0; i < hpoTerms.length; i++) {
      hpoHtml = hpoHtml + PP.create_url('/hpo', hpoTerms[i].hpo_term, hpoTerms[i].hpo_id, 'chip');
    }
    return (hpoHtml);
  };

  PP.generateVariantHtml = function(variant) {
    return PP.create_url('/variant', variant) + '<br><a class="text-small" href="/sequence?variant_id=' + variant + '&symbol=&build=grch37">Primer Design</a>';
  };

  PP.generateIndividualHtml = function(individuals) {
    return individuals.join(' ');
  };

  PP.createHgvspTd = function(hgvsp, hgvsc) {
    if (hgvsp === undefined) return;
    if (hgvsc === undefined) return;
    return '<span data="' + hgvsc[0] + '">' + hgvsp[0].split(':')[1] + '</span>';
  };


  PP.initExomiserTab = function(patient_id) {
    // $.ajax({
    //   type: 'GET',
    //   url: '/homozygous_variants_json/',
    //   dataType: 'json',
    //   timeout: 120000,
    //   success: function(data) { PP.HomsSuccess(data); },
    //   error: function(data, msg) { PP.ajaxError(data, msg); }
    // });

  };

}());
// End of PP module
