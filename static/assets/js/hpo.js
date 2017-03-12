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

  ////
  //  Fetch Individual Table (HPO Page)
  ////

  // Adds the side navbar dropdown
  PP.initIndividualsTable = function(hpo_id) {
    $.ajax({
      type: 'GET',
      url: '/hpo_individuals_json/' + hpo_id,
      dataType: 'json',
      timeout: 120000,
      success: function(data) { PP.individualsSuccess(data); },
      error: function(data, msg) { PP.ajaxError(data, msg); }
    });
  };

  //
  PP.individualsSuccess = function(data) {
    var individuals = data.result.individuals;
    for (var i = 0; i < individuals.length; i++) {
      $('<tr>').append(
        $('<td>').html(PP.generateIndividualLinkHtml(individuals[i].external_id)), // individual_cell
        $('<td>').text(individuals[i].sex), //gender_cell
        $('<td>').html(PP.generateIndividualPhenotypeCell(individuals[i].features)), //phenotype_cell
        $('<td>').text(individuals[i].specificity.score), // phenotype_score_cell          
        $('<td>').text(individuals[i].rare_homozygous_variants_count), // rare_hom_cell
        $('<td>').text(individuals[i].rare_compound_hets_count), // rare_comp_cell
        $('<td>').text(individuals[i].rare_variants_count), // rare_cell
        $('<td>').html(PP.generateIndividualGeneCell(individuals[i].genes)), // candidate_genes_cell
        $('<td>').html(PP.generateIndividualSolvedVariantsCell(individuals[i].solved_variants, individuals[i].genes)) // candidate_variants_cell
      ).appendTo('#individuals_table_body');
    }
    $('#individual_tab #progress_row').remove();
    $('#individuals_table').show();
    PP.initTableSorter('#individuals_table');
    PP.addIndividualPopOver('#individuals_table', 'right');
  };

  //
  PP.generateIndividualLinkHtml = function(external_id) {
    if (external_id === 'hidden') {
      return '<a href="#!">' + external_id + '</a>';
    } else {
      return PP.create_url('/individual', external_id, external_id, 'individual_popover');
    }
  };

  //
  PP.generateIndividualPhenotypeCell = function(features) {
    var phenotype_html = '';
    for (var j = 0; j < features.length; j++) {
      phenotype_html = phenotype_html + ' ' + PP.create_url('/hpo', features[j].label, features[j].id, 'chip');
    }
    return (phenotype_html);
  };

  //
  PP.generateIndividualGeneCell = function(genes) {
    var genes_html = '';
    for (var j = 0; j < genes.length; j++) {
      genes_html = genes_html + ' ' + PP.create_url('/gene', genes[j]);
    }
    return (genes_html);
  };

  //
  PP.generateIndividualSolvedVariantsCell = function(solved_variants, genes) {
    var solved_variants_html = '';
    for (var j = 0; j < genes.length; j++) {
      if (solved_variants === undefined) return;
      if (solved_variants.genes === undefined) return;
      if (solved_variants.genes[j] === undefined) return;
      if (solved_variants.genes[j].hom === undefined) {
        solved_variants_html = solved_variants_html + ' [ , ' + PP.create_url('/variant', solved_variants.genes[j].hom) + '] ';
      } else if (solved_variants.genes[j].het === undefined) {
        solved_variants_html = solved_variants_html + ' [' + PP.create_url('/variant', solved_variants.genes[j].het) + ', ] ';
      } else {
        solved_variants_html = solved_variants_html + ' [' + PP.create_url('/variant', solved_variants.genes[j].het) + ', ' + PP.create_url('/variant', solved_variants.genes[j].hom) + '] ';
      }
    }
    return (solved_variants_html);
  };

  PP.addLoadingTooltips = function(id, direction) {
    var yellowPreloaderHtml = '<div class="preloader-wrapper active"><div class="spinner-layer spinner-yellow-only"><div class="circle-clipper left"><div class="circle"></div></div><div class="gap-patch"><div class="circle"></div></div><div class="circle-clipper right"><div class="circle"></div></div></div></div>';
    $(id).tooltip({ tooltip: yellowPreloaderHtml, position: direction, html: true });
  };

  //
  PP.addIndividualPopOver = function(table, direction) {
    PP.addLoadingTooltips('.individual_popover', direction);
    $(table).bind('pagerComplete', function() {
      PP.addLoadingTooltips('.individual_popover', direction);
    });

    $(table).hoverIntent({
      selector: '.individual_popover',
      interval: 500, // number = milliseconds delay before trying to call over          
      timeout: 500, // number = milliseconds delay before onMouseOut    
      over: function() {
        if (!$(this).is('[data-tooltip]')) {
          var individual_id = $(this).text();
          var e = this;
          $.get('/phenogenon_json/' + individual_id, function(data) {
            var content = PP.generateTooltipContent(data, individual_id);
            $(e).attr('data-tooltip', content);
            $(e).tooltip({ position: direction, html: true });
            $(e).trigger("mouseenter.tooltip");
          });
        }
      },
      out: function() {
        $(this).trigger("mouseleave.tooltip");
      }
    });
  };

  PP.generateTooltipContent = function(data, individual_id) {
    var content = '<span class=yellow-text>' + individual_id + '</span><br>';
    if (data.result === undefined) {
      return (content + '<span> No Data Found</span>');
    }
    var status, features;
    if (data.result.solved !== undefined) {
      status = data.result.solved.status;
    } else {
      status = 'Unknown Status';
    }
    if (data.result.featues !== undefined) {
      features = PP.getSortedIndividualFeatures(data.result.features);
    } else {
      features = 'Unknown Features';
    }    
    content = content + '<span>Status: ' + status + ' </span><br>' +
      '<span>Features: ' + features + ' </span>';
    return (content);
  };


  PP.getSortedIndividualFeatures = function(features) {
    var sortedFeatures = features.sort();
    var featuresHtml = '';
    for (var i = 0; i < features.length; i++) {
      featuresHtml = featuresHtml + PP.create_url('/hpo', features[j].label, features[j].label, 'chip');
    }
    return (featuresHtml);
  };

  ////
  //  Fetch Literature Genes & Phenogenon Table (HPO Page)
  ////
  //
  PP.initLitGenesAndPhenogenonTables = function(hpo_id) {
    $.ajax({
      type: 'GET',
      url: '/phenogenon_json/' + hpo_id,
      dataType: 'json',
      timeout: 120000,
      success: function(data) {
        PP.litGenesSuccess(data);
        PP.dominantPhenogenonSuccess(data);
        PP.recessivePhenogenonSuccess(data);
      },
      error: function(data, msg) { PP.ajaxError(data, msg); }
    });
  };

  //
  PP.litGenesSuccess = function(data) {
    var lit_genes = data.result.lit_genes;
    for (var i = 0; i < lit_genes.length; i++) {
      $('<tr>').append(
        $('<td>').html(PP.create_url('/gene', lit_genes[i].gene_name)), // Gene
        $('<td>').text(lit_genes[i].phenogenon.het.unrelated_dominant_all_p_val), // dominant pvalue
        $('<td>').text(lit_genes[i].phenogenon.hom_comp.unrelated_recessive_p_val) // recessive p value
      ).appendTo('#phenogenon_lit_genes_body');
    }
    $('#literature_genes #progress_row').remove();
    $('#phenogenon_lit_genes_table').show();
    PP.initTableSorter('#phenogenon_lit_genes_table');
  };

  //
  PP.dominantPhenogenonSuccess = function(data) {
    var dominant_genes = data.result.dominant_genes;
    for (var i = 0; i < dominant_genes.length; i++) {
      // generate tr
      $('<tr>').append(
        $('<td>').html(PP.generatePhenogenonGeneNameCell(dominant_genes[i])),
        $('<td>').text(dominant_genes[i].p_val)
      ).appendTo('#dominant_table_body');
    }
    $('#phenogenon_dominant #progress_row').remove();
    $('#gene_dominant_table').show();
    PP.initTableSorter('#gene_dominant_table');
  };

  //
  PP.recessivePhenogenonSuccess = function(data) {
    var recessive_genes = data.result.recessive_genes;
    for (var i = 0; i < recessive_genes.length; i++) {
      // generate tr
      $('<tr>').append(
        $('<td>').html(PP.generatePhenogenonGeneNameCell(recessive_genes[i])),
        $('<td>').text(recessive_genes[i].p_val)
      ).appendTo('#recessive_table_body');
    }
    $('#phenogenon_recessive #progress_row').remove();
    $('#gene_recessive_table').show();
    PP.initTableSorter('#gene_recessive_table');
  };

  PP.generatePhenogenonGeneNameCell = function(genes) {
    var gene_name;
    if (genes.known) {
      gene_name = PP.create_url('/gene', genes.gene_name, genes.gene_name, 'known_gene');
    } else {
      gene_name = PP.create_url('/gene', genes.gene_name);
    }
    return (gene_name);
  };

  ////
  //  Fetch Skats Table (HPO Page)
  ////
  //
  PP.initSkatsTable = function(hpo_id) {
    $.ajax({
      type: 'GET',
      url: '/hpo_skat_json/' + hpo_id,
      dataType: 'json',
      timeout: 120000,
      success: function(data) { PP.skatsSuccess(data); },
      error: function(data, msg) { PP.ajaxError(data, msg); }
    });
  };

  // 
  PP.skatsSuccess = function(data) {
    var skats = data.result.individuals;
    for (var i = 0; i < skats.length; i++) {
      // generate tr
      $('<tr>').append(
        $('<td>').html(PP.create_url('/gene', skats[i].Symbol)),
        $('<td>').text(skats[i].pli.toFixed(1)),
        $('<td>').text(skats[i].mode),
        $('<td>').text(skats[i].FisherPvalue),
        $('<td>').text(skats[i].SKATO),
        $('<td>').text(skats[i].OddsRatio),
        $('<td>').html(PP.generateSkatsVariantsCell(skats[i].variants))
      ).appendTo('#skats_table_body');
    }
    $('#skats #progress_row').remove();
    $('#skats_table').show();
    PP.initTableSorter('#skats_table');
    PP.addVariantPopOver('#skats_table', 'left');
  };

  //
  PP.generateSkatsVariantsCell = function(variants) {
    var variant_html = '';
    for (var j = 0; j < variants.length; j++) {
      variant_html = variant_html + ' ' + PP.create_url('/variant', variants[j], variants[j], 'chip variantId_popover');
    }
    return (variant_html);
  };

  //
  PP.addVariantPopOver = function(table, direction) {
    PP.addLoadingTooltips('.variantId_popover', direction);
    $(table).bind('pagerComplete', function() {
      PP.addLoadingTooltips('.variantId_popover', direction);
    });

    $(table).hoverIntent({
      selector: '.variantId_popover',
      interval: 500, // number = milliseconds delay before trying to call over          
      timeout: 500, // number = milliseconds delay before onMouseOut    
      over: function() {
        if (!$(this).is('[data-tooltip]')) {
          var variant_id = $(this).text();
          var e = this;
          $.get('/variant_json/' + variant_id, function(data) {
            var content = '<span class=yellow-text>' + variant_id + '</span><br>' +
              '<span>Hom Count: ' + data.result.HOM_COUNT + ' </span><br>' +
              '<span>Het Count: ' + data.result.HET_COUNT + ' </span><br>' +
              '<span>Miss Count: ' + data.result.MISS_COUNT + ' </span>';
            $(e).attr('data-tooltip', content);
            $(e).tooltip({ position: direction, html: true });
            $(e).trigger("mouseenter.tooltip");
          });
        }
      },
      out: function() {
        $(this).trigger("mouseleave.tooltip");
      }
    });
  };

  //
  PP.initPhenogenonTabs = function() {
    $($('ul.tabs .tab')[2]).on('click', function() {
      if (!$('#phenogenon_tabs').hasClass('tabs_initialised')) {
        setTimeout(
          function() {
            $('#phenogenon_tabs').tabs();
            $('#phenogenon_tabs').addClass('tabs_initialised');
          }, 1);
      }
    });
  };

}());
// End of PP module
