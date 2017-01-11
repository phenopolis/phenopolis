// this file includes functions to transform dot file data to dynamically display gene-hpo result.

// please check uclex_browser/draw_hpo_graph.py for details
function remove_duplicates(data){
    // remove duplicates in an array
    var seen = {};
    var result = [];
    $.each(data, function(i,v){
        if (!seen[v]){
            result.push(v);
            seen[v] = true;
        }
    });
    return result;
}
function compute_qval(data){
    //Benjamini–Hochberg
    var p_values_arr = [];
    var p_values_dict = {};
    $.each(data, function(key,val){
        p_values_arr.push(val['p_value']);
    });
    // sort p value
    p_values_arr = p_values_arr.sort();
    // calculate q value
    var total_num_of_observations = p_values_arr.length;
    $.each(p_values_arr,function(i,val){
        if (!(val in p_values_dict)){
            p_values_dict[val] = val * total_num_of_observations / (i+1);
        }
    });
    // assign q values
    $.each(data, function(key,val){
        val.q_value = p_values_dict[val.p_value];
    });
    return data;
}
function transform(input_data, hpo_freq, mode, exac_cutoff, cadd_cutoff, exclude_hom){
    // total number of patient hpo appearance in the whole cohort.
    var total = 0;
    $.each(hpo_freq, function(key, value){
        total = total + value.size;
    });
    
    // re-calcluate input_data according to exac_cutoff
    var output = {};
    // hpo_patients = total patient hpo appearance in this hpo graph
    var hpo_patients = [];
    var num_of_vars = 0;
    if (mode == 'hom_comp'){
        // recalculate observed cases for each hpo term
        $.each(input_data, function(key,value){
            // note that each patient has to have at least 2 variants that fit the cutoff
            var observed_cases = [];
            var val = {};
            $.each(value.data, function(k, list){
                var ok = [];
                $.each(list['var'], function(i, af){
                    if (af.cadd == null){
                        // frameshift don't have af.cadd. make it big
                        af.cadd = 100;
                    }
                    if (af.exac <= exac_cutoff && af.cadd >= cadd_cutoff){
                        ok= $.merge(ok, [af]);
                    }
                });
                if (ok.length > 1){
                    // this patient is good
                    num_of_vars = num_of_vars + ok.length;
                    observed_cases.push(k);
                    val[k] = {'hpo':list.hpo, 'var':ok};
                }
            });
            if (observed_cases.length){
                output[key] = value;
                output[key].data = val;
                // change patient_size
                output[key]['patient_size'] = observed_cases.length;
                hpo_patients = $.merge(hpo_patients, observed_cases);
            }
        });
    } else if (mode == 'het'){
        // recalculate observed cases for each hpo term
        // if exclude_dominant is true, remove patients who have more than one variants
        $.each(input_data, function(key,value){
            var observed_cases = {};
            var val= {};
            $.each(value.data, function(k,list){
                $.each(list['var'], function(i, af){
                    if (af.cadd == null){
                        // frameshift don't have af.cadd. make it big
                        af.cadd = 100;
                    }
                    if (af.exac <= exac_cutoff && af.cadd >= cadd_cutoff){
                        if (k in observed_cases){
                            observed_cases[k]+=1
                        }else{
                            observed_cases[k]=1;
                        }
                        num_of_vars = num_of_vars + 1;
                        val[k] = val[k] || {'hpo':list.hpo,'var':[]};
                        val[k]['var'].push(af);
                    }
                });
            });
            // remove more key with value > 1 if exclude_dominant true
            if (exclude_hom){
                $.each(observed_cases, function(kk,vv){
                    if (vv > 1){
                        delete observed_cases[kk];
                        delete val[kk];
                    }
                });
            }
            observed_cases = Object.keys(observed_cases);
            if (observed_cases.length){
                output[key] = value;
                output[key].data = val;
                // change patient_size
                output[key]['patient_size'] = observed_cases.length;
                hpo_patients = $.merge(hpo_patients,observed_cases);
            }
        });
    } else {
        console.log('ERROR: need to pass "hom_comp" or "het" as mode to transform!');
        return output;
    }
    // get unique paitents
    hpo_patients = remove_duplicates(hpo_patients);
    pat_g = hpo_patients.length;
    var num_of_hpos = Object.keys(output).length;
    // after filtering, recalculate the p-value and q-value

    $.each(output, function(key, value){
        pat_h = hpo_freq[key]['raw'].split('/')[0];
        pat_a = hpo_freq[key]['raw'].split('/')[1];
        r_ga = pat_g/pat_a;
        //value['p_value'] = 1-binomial_cdf(value['patient_size']-1,num_of_patients,hpo_freq[key]['freq']);
        value['p_value'] = 1-binomial_cdf(value['patient_size']-1,pat_h,r_ga);
        value['width'] = 0.1 * Math.log(1/value['p_value']);
        //value['observed_freq'] = value['patient_size'] + '/' + num_of_patients;
        value['observed_freq'] = value['patient_size'] + '/' + pat_h;
        value['expect_freq'] = pat_g + '/' + pat_a;
    });
    // compute q value
    output = compute_qval(output);
    result = [pat_g, num_of_vars, output];
    return result;
}

function convert_to_dot(data,style){
    // convert data to dot string for visulisation
    var dot = 'digraph {\nsize="6,4";\n';
    // add legend if style=wedged, ie from irdc_summary
    var centre_color = {'OXF':{'color':'green'},'LON':{'color':'red'},'LDS':{'color':'steelblue'},'MAN':{'color':'orange'}};
    if (style == 'wedged'){
        $.each(centre_color,function(k,e){
            dot = dot +'"' + k + '" [style="filled", fixedsized="true", fontsize="90", shape="square", width="3", fillcolor="'+e.color+'",label="'+k+'",color="transparent"];\n';
        });
    }

    $.each(data, function(key,d){
        d.label = d.name;
        if (d.width < 0.4){
            d.label = '';
        }
        if (d.p_value == 0){
            d.width = 0.001
        }
        if (d.q_value <= 0.01){
            d.label = d.name + '**';
        }else if (d.q_value <= 0.05){
            d.label = d.name + '*';
        }
        var fillcolor = [];
        if (style == 'filled'){
            // draw gene-hpo graph
            fillcolor = [d.fillcolor];
        }else if (style == 'wedged'){
            // draw irdc summary
            $.each(centre_color,function(k,e){
                var total = d.observed_freq.split('/')[0];
                var propotion = d[k].split('/')[0]/total;
                fillcolor.push(e.color+';'+propotion);
            });
            // calculate proportion

        }
        dot = dot + '"' + d.id +  '" [style="'+style+'", fixedsize="true", fontsize="6", shape="circle", width="' + d.width +'", fillcolor="' + fillcolor.join(':') +'", label="' + d.label +'", id="' + d.id + '", color="transparent"];\n';
        //dot = dot + '"' + d.id +  '" [style="wedged", fixedsize="true", fontsize="6", shape="circle", width="' + d.width +'", fillcolor="red;0.3:green;0.6:orange", label="' + d.label +'", id="' + d.id + '", color="transparent"];\n';
    });
    // add links
    $.each(data, function(key,value){
        if ('is_a' in value){
            $.each(value['is_a'], function(i, is_a){
                dot = dot + '"' + is_a + '" -> "' + key +'" [color=grey, lty="solid"];\n';
            });
        }
    });
    return dot + '}';
};
