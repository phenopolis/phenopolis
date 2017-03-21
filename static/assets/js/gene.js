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
    /*
     *  Gene Page  FUNCTIONS 
     */

    PP.generateGenePageTable = function(input_data, hpo_data, mode, remove_hom=false, remove_bad_vars=true) {
        if ($.isEmptyObject(input_data))  return true; // Don't do anything if empty input_data
        var data_copy = jQuery.extend(true, {}, input_data); // Create a copy of the input_data object

        returned = PP.transformData(data_copy, hpo_data, mode, remove_hom, remove_bad_vars);
        var data = returned.data;
        data = PP.dictToSortedArray(data, 'p_value');
        for (var i = 0; i < data.length; i++) {
          $('<tr>').append(
            $('<td>').html(PP.create_url('/hpo/', data[i].name, data[i].id)),
            $('<td>').text(data[i].phi),
            $('<td>').text(data[i].p_value),
            $('<td>').text(data[i].pat_h),
            $('<td>').text(data[i].pat_gh)
          ).appendTo('#'+mode+'-table_body');
          if (hpo_data.indexOf(data[i].id) != -1) {  //highlight
            $('#'+mode+'-table_body tr:last-child').addClass('highlight');
          }
        }
        $('#' + mode + '_small').text('(Pat_a:' + returned.pat_a + ', Pat_g:' + returned.pat_g + ')');
        PP.initTableSorter('gene-hpo-' + mode + '-table');
    };

    // remove_bad_vars determine if dubious variants getting removed.
    // since ppl with the variants marked as dubious have better changes to have the rare variants than average ppl, it's better to remove the patients out of the system. 
    // For dominant cases, if one patient has two qualified variants, but one of them is dubious, when calculating for dominant mode, shall we include this patient?
    // For now, remove the patient completely from the system.
    // collect qualified patients/vars from HP:0000001 (which includes info on all HP:*)
    PP.transformData = function(input_data, hpo_data, mode, exclude_hom, remove_bad_vars) {
        var qualified_patients = {},
            qualified_vars = {},
            dubious_patients = {};
        var unrelated_flag = $('#remove-related').is(':checked');
        var exac_cutoff = $('#'+mode+'-exac_slider').val();
        var cadd_cutoff = $('#'+mode+'-cadd_slider').val();
        $.each(input_data['HP:0000001'].data, function(k1, v1) {
            if (unrelated_flag && !v1.unrelated) return true; // check related?
            // count good vars
            var count = 0,
                subcount = 0,
                dubious = 0;
            $.each(v1.var, function(k2, v2) {
                if (v2.exac <= exac_cutoff && (!v2.cadd || v2.cadd >= cadd_cutoff)) {
                    if (!remove_bad_vars || !('filter' in v2) || v2.filter == 'PASS') {
                        count += 1;
                        qualified_vars[v2.variant_id] = 1;
                    } else {
                        subcount += 1;
                    }
                }
            });
            // qualified?
            if (subcount && ((mode == 'het' && !(exclude_hom) && count === 0) ||
                    (mode == 'het' && exclude_hom && count <= 1) || count < 2)) {
                // need to remove this patient out of the equation
                dubious = 1;
                dubious_patients[k1] = 1;
            }
            if (dubious) {
                return true;
            }
            if ((mode == 'hom' && count > 1) ||
                (mode == 'het' && !(exclude_hom) && count) ||
                (mode == 'het' && exclude_hom && count == 1)) {
                qualified_patients[k1] = 1;
            }
        });
        // re-calcluate input_data according to exac_cutoff
        var output = {};
        var num_of_vars = Object.keys(qualified_vars).length;
        var pat_g = Object.keys(qualified_patients).length;
        $.each(input_data, function(k1, v1) {
            // k1:hpo, k2:patient
            output[k1] = v1;
            output[k1].bad_gh = 0;
            var newData = {};
            $.each(v1.data, function(k2, v2) {
                if (k2 in dubious_patients) {
                    output[k1].bad_gh += 1;
                }
                if (!(k2 in qualified_patients)) {
                    return true;
                }
                var newVars = $.grep(v2.var, function(v3, i) {
                    return v3.variant_id in qualified_vars;
                });
                newData[k2] = { 'hpo': v2.hpo, 'var': newVars, 'unrelated': v2.unrelated };
            });
            output[k1].pat_gh = Object.keys(newData).length;
            if (output[k1].pat_gh === 0) {
                // empty node
                delete output[k1];
                return true;
            }
            output[k1].data = newData;
        });

        var num_of_hpos = Object.keys(output).length;
        // after filtering, recalculate the p-value and q-value
        var pat_a = 0;
        $.each(output, function(key, value) {
            var pat_h = value.unrelated_pat_h - value.bad_gh;
            if (!(unrelated_flag)) {
                pat_h = value.related_pat_h - value.bad_gh;
            }
            value.pat_h = pat_h;
            pat_a = (unrelated_flag ? value.unrelated_pat_a : value.related_pat_a) - Object.keys(dubious_patients).length;
            var r_ga = pat_g / pat_a;
            value.pat_a = pat_a;
            value.p_hg = 1 - PP.binomial_cdf(value.pat_gh - 1, pat_g, pat_h / pat_a);
            value.p_gh = 1 - PP.binomial_cdf(value.pat_gh - 1, pat_h, r_ga);
            value.phi = PP.phi([pat_a - pat_g - pat_h + value.pat_gh, pat_h - value.pat_gh, pat_g - value.pat_gh, value.pat_gh]);
            value.p_fisher = PP.exact22(pat_a - pat_g - pat_h + value.pat_gh, pat_h - value.pat_gh, pat_g - value.pat_gh, value.pat_gh);
            value.p_value = value.p_fisher.right;
            if (value.id == 'HP:0000001') {
                value.width = 0.4;
            } else {
                value.width = 0.1 * Math.log(1 / value.p_value);
            }
            //value.observed_freq = value.pat_gh + '/' + num_of_patients;
            value.observed_freq = value.pat_gh + '/' + pat_h;
            value.expected_freq = pat_g + '/' + pat_a;
            if (value.id == 'HP:0000001') {
                value.fillcolor = 'green';
            } else {
                value.fillcolor = '#44ccff';
                if (hpo_data.indexOf(value.id) != -1) {
                    value.fillcolor = 'red';
                }
            }
        });
        // compute q value
        output = PP.compute_qval(output);
        result = { 'pat_a': pat_a, 'pat_g': pat_g, 'num_of_vars': num_of_vars, 'data': output };
        return result;
    };

    PP.compute_qval = function(data) {
        //Benjamini–Hochberg
        var p_values_arr = [];
        $.each(data, function(key, val) {
            p_values_arr.push(val.p_value);
        });
        p_values_arr = p_values_arr.sort();  // sort p value

        // calculate q value
        var p_values_dict = {};
        var total_num_of_observations = p_values_arr.length;
        $.each(p_values_arr, function(i, val) {
            if (!(val in p_values_dict)) {
                p_values_dict[val] = val * total_num_of_observations / (i + 1);
            }
        });
        // assign q values
        $.each(data, function(key, val) {
            val.q_value = p_values_dict[val.p_value];
        });
        return data;
    };

    // transfrom dict to array in order to be drawn by drawTable. Sort by p value by default
    PP.dictToSortedArray = function(dict, sortedBy) {
        result = [];
        $.each(dict, function(k, v) {
            // format phi and p_value
            v.phi = PP.formatFloat(v.phi, 4);
            v.p_value = PP.formatFloat(v.p_value, 4);
            result.push(v);
        });
        //sort
        result.sort(function(a, b) {
            return ((a[sortedBy] < b[sortedBy]) ? -1 : ((a[sortedBy] > b[sortedBy]) ? 1 : 0));
        });
        return result;
    };

    // format float
    PP.formatFloat = function(v, n) {
        var str = v.toExponential();
        var [b, e] = str.split('e');
        b = (+b).toFixed(n - 1);
        var result = b + 'e' + e;
        return +result;
    };

}());
// End of PP module
