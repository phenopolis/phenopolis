$(document).ready(function() {
    if ($(window).width() < 768) {
        $('#gene_plot_container').css('width', $(window).width() + "px");
    } else {
        $('#gene_plot_container').css('width', $(window).width()*10/12 + "px");
    }
    precalc_coding_coordinates(window.transcript, window.coverage_stats, 'pos');
    precalc_coding_coordinates(window.transcript, window.variants_in_transcript, 'pos');

    // only show variants that have a coding coordinate
    window.variants_in_transcript = _.filter(window.variants_in_transcript, function(variant) {
        return variant.pos_coding != undefined;
    });

    // only show coding rects that have a coding coordinate
    window.coverage_stats = _.filter(window.coverage_stats, function(d) {
        return d.pos_coding != undefined;
    });
    $('#avg_coverage').html(coverage_sum('-log10pvalue'));
    $('#avg_coverage_x').html(coverage_sum('30')*100 + '%');

    var new_data = create_new_data(window.coverage_stats, 'pos_coding');
    var new_data_skip_utr = create_new_data(window.coverage_stats, 'pos_coding_noutr');
    if (window.coverage_stats != null) {
        gene_chart(window.coverage_stats, new_data_skip_utr, window.variants_in_transcript, window.transcript);
        if (window.variants_in_transcript.length) {
            update_variants();
        }
        $('#loading_coverage').hide();
    } else {
        $('#gene_plot').hide();
        $('#not_covered').show();
    }
    // Change coverage plot
    $('.coverage_metric_buttons').change(function () {
        var v = $(this).attr('id').replace('_covmet_button', '');
        $('.coverage_subcat_selectors').hide();
        if (v == 'covered') {
           $('#over_x_select_container').show();
            v = $('#over_x_select').val().replace('X', ''); 
        } else {
            $('#average_select_container').show();
            v = $("#average_select").val();
        }
        var detail = $('.display_coverage_metric_buttons.active').attr('id').replace('display_coverage_', '').replace('_button', '');
        var include_utrs = $('#include_utrs_checkbox').is(':checked');
        var plot_data = include_utrs ? new_data : new_data_skip_utr;
        change_coverage_chart(window.coverage_stats, plot_data, window.variants_in_transcript, window.transcript, detail, v, !include_utrs, '#gene_plot_container');
        refresh_links();
    });
    $('#over_x_select').change(function () {
        var detail = $('.display_coverage_metric_buttons.active').attr('id').replace('display_coverage_', '').replace('_button', '');
        var include_utrs = $('#include_utrs_checkbox').is(':checked');
        var plot_data = include_utrs ? new_data : new_data_skip_utr;
        $('#avg_coverage_type_x').html($(this).val());
        $('#avg_coverage_x').html(coverage_sum($(this).val().replace('X', ''))*100 + '%');
        change_coverage_chart(window.coverage_stats, plot_data, window.variants_in_transcript, window.transcript, detail, $(this).val().replace('X', ''), !include_utrs, '#gene_plot_container');
        refresh_links();
    });
    $('#average_select').change(function () {
        var detail = $('.display_coverage_metric_buttons.active').attr('id').replace('display_coverage_', '').replace('_button', '');
        var include_utrs = $('#include_utrs_checkbox').is(':checked');
        var plot_data = include_utrs ? new_data : new_data_skip_utr;
        $('#avg_coverage_type').html($(this).val());
        $('#avg_coverage').html(coverage_sum($(this).val()));
        change_coverage_chart(window.coverage_stats, plot_data, window.variants_in_transcript, window.transcript, detail, $(this).val(), !include_utrs, '#gene_plot_container');
        refresh_links();
    });
    $('#include_utrs_checkbox').change(function () {
        setTimeout(function() {
            var detail = $('.display_coverage_metric_buttons.active').attr('id').replace('display_coverage_', '').replace('_button', '');
            var v = $('.coverage_metric_buttons.active').attr('id').replace('_covmet_button', '');
            v = (v == 'covered') ? $('#over_x_select').val().replace('X', '') : $("#average_select").val(); 
            var include_utrs = $('#include_utrs_checkbox').is(':checked');
            var plot_data = include_utrs ? new_data : new_data_skip_utr;
            change_coverage_chart(window.coverage_stats, plot_data, window.variants_in_transcript, window.transcript, detail, v, !include_utrs, '#gene_plot_container');
            refresh_links();
        }, 10);
    });

// Change exon diagram
    $('#inverted_checkbox').change(function () {
        setTimeout(function () {
            var v = $('#inverted_checkbox').is(':checked');
            change_track_chart_variant_size(window.variants_in_transcript, v, '#gene_plot_container');
            refresh_links();
        }, 10);
    });

    $('.consequence_display_buttons, #filtered_checkbox').change(function () {
        setTimeout(function() {
            update_variants();
            refresh_links();
        }, 10);
    });

    $('.display_coverage_metric_buttons').change(function () {
        var detail = $(this).attr('id').replace('display_coverage_', '').replace('_button', '');
        var v = $('.coverage_metric_buttons.active').attr('id').replace('_covmet_button', '');
        if (v == 'covered') {
            $('#over_x_select_container').show();
            v = $('#over_x_select').val().replace('X', ''); 

        } else {
            $('#average_select_container').show();
            v = $("#average_select").val();
        }
        var include_utrs = $('#include_utrs_checkbox').is(':checked');
        var plot_data = include_utrs ? new_data : new_data_skip_utr;
        change_coverage_chart(window.coverage_stats, plot_data, window.variants_in_transcript, window.transcript, detail, v, !include_utrs, '#gene_plot_container');
        refresh_links();
    });
    refresh_links();
});


