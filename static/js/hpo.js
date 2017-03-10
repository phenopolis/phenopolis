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
var d3color = d3.scale.linear()
    .domain([-0.2,0,0.2])
    .range(['steelblue','lightgrey','red']);

function transform(input_data, options){
    // remove_bad_vars determine if dubious variants getting removed.
    // since ppl with the variants marked as dubious have better changes to have the rare variants than average ppl, it's better to remove the patients out of the system. 
    // For dominant cases, if one patient has two qualified variants, but one of them is dubious, when calculating for dominant mode, shall we include this patient?
    // For now, remove the patient completely from the system.

    // what kind of cutoff i'm looking at?
    var exac = options.mode == 'd' ? 'exac_af' : 'exac_hom';

    var qualified_patients = {},
        qualified_vars = {},
        dubious_vars = {},
        dubious_patients = {};
    //$.each(input_data.data.['HP:0000001'].data, function(k1,v1){
    // get qualified and dubious vars
    $.each(input_data.variants, function(k1,v1){
        v1.qualified = false;
        if (v1[exac] <= options.cutoffs.exac && (!v1.cadd || v1.cadd >= options.cutoffs.cadd)){
            if (!options.remove_bad_vars || v1.filter == 'PASS'){
                qualified_vars[k1] = 1;
                v1.qualified = true;
            } else {
                dubious_vars[k1] = 1;
            }
        }
    });
        
    $.each(input_data.patients, function(k1,v1){
        // related?
        if (options.unrelated_flag && !v1.unrelated){
            return true;
        }
        v1.qualified = false;
        // count good vars
        var count = 0;
        var subcount = 0;
        var dubious = 0; // dubious patient?
        $.each(v1.variants, function(k2,v2){
            if (qualified_vars.hasOwnProperty(v2)){
                count += 1;
            } else if (dubious_vars.hasOwnProperty(v2)){
                subcount += 1;
            }
        });
        // qualified?
        if (subcount && ( (options.mode == 'd' && !(options.exclude_hom) && count == 0 ) 
                        || (options.mode == 'd' && options.exclude_hom && count <= 1)
                        || count < 2)
            ){
            // need to remove this patient out of the equation
            dubious = 1;
            dubious_patients[k1] = 1
        }
        if (dubious){ return true; }
        if ( (options.mode == 'r' && count > 1) ||
             (options.mode == 'd' && !(options.exclude_hom) && count) ||
             (options.mode == 'd' && options.exclude_hom && count == 1) ){
                 qualified_patients[k1] = 1;
                 v1.qualified = true;
        }
    });
    // re-calcluate input_data according to cutoffs
    var output = { },
        num_of_vars = Object.keys(qualified_vars).length,
        pat_g = Object.keys(qualified_patients).length;
    $.each(input_data.data, function(k1,v1){
        // k1:hpo, v2:patient
        output[k1] = {
            'id':k1,
            'name':v1.name,
            'is_a':v1.is_a,
            'unrelated_pat_h':v1['unrelated_pat_h'],
            'related_pat_h':v1['related_pat_h'],
        };
        output[k1].bad_gh = 0;
        output[k1].pat_gh = 0;
        $.each(v1.p[options.mode], function(k2,v2){
            if (dubious_patients.hasOwnProperty(v2)){
                output[k1].bad_gh += 1;
            }
            if (!(qualified_patients.hasOwnProperty(v2))){
                return true;
            }
            output[k1].pat_gh += 1;
        });
        if (output[k1].pat_gh==0){
            // empty node
            delete output[k1];
            return true;
        }
                
    });
    var num_of_hpos = Object.keys(output).length;
    // after filtering, recalculate the p-value and q-value
    var pat_a = (options.unrelated_flag ? input_data.unrelated_pat_a : input_data.related_pat_a) - Object.keys(dubious_patients).length;
    $.each(output, function(key, value){
        var pat_h = value['unrelated_pat_h']-value.bad_gh;
        if (!(options.unrelated_flag)){
            pat_h = value['related_pat_h']-value.bad_gh;
        }
        value['pat_h'] = pat_h;
        var r_ga = pat_g/pat_a;
        value['pat_a'] = pat_a;
        //value['pat_h'] = pat_h;
        value['phi'] = phi([pat_a-pat_g-pat_h+value['pat_gh'],pat_h-value['pat_gh'],pat_g-value['pat_gh'],value['pat_gh']]);
        value['p_fisher'] = exact22(pat_a-pat_g-pat_h+value['pat_gh'],pat_h-value['pat_gh'],pat_g-value['pat_gh'],value['pat_gh']);
        value['p_value'] = value['p_fisher']['right'];
        if (value.id == 'HP:0000001'){
            value['width'] = 0.4;
        } else {
            value['width'] = 0.1 * Math.log(1/value['p_value']);
        }
        //value['observed_freq'] = value['pat_gh'] + '/' + num_of_patients;
        value['observed_freq'] = value['pat_gh'] + '/' + pat_h;
        value['expected_freq'] = pat_g + '/' + pat_a;
        if (value.id == 'HP:0000001'){
            value['fillcolor'] = 'green';
        }else{ 
            value['fillcolor'] = '#44ccff';
            if (options.hpo_terms.indexOf(value.id) != -1){
                value['fillcolor'] = 'red';
            }
        }
    });
    // compute q value
    output = compute_qval(output);
    result = {'pat_a':pat_a,'pat_g':pat_g, 'num_of_vars':num_of_vars, 'data':output};
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
        var fillcolor = [d.fillcolor];
        if (style == 'filled'){
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
function transform_overview(my_data,breakdown){
    // work on overview of uclex patients
    // breakdown on sex? solved?
    var cutoff = 1/20; // if a node has less than cutoff ratio of total, don't return
    var temp = {};
    var total = my_data['HP:0000001'].count;
    $.each(my_data,function(k1,v1){
        if (v1.count / total < cutoff){return true;}
        temp[k1] = {
            'count':v1.count,
            'freq':v1.freq,
            'is_a':v1.is_a,
            'name':v1.name
            };
        if (breakdown == 'cohort'){
            temp[k1].cohort = v1.cohort;
        } else if (breakdown == 'sex'){
            temp[k1].sex = {'M':{},'F':{},'U':{}};
            $.each(v1.cohort,function(k2,v2){
                $.each(v2,function(k3,v3){
                    temp[k1].sex[k3][k2] = v3;
                });
            });
        } else if (breakdown == 'solved'){
            temp[k1].solved = {'solved':{},'unsolved':{},'candidate':{}};
            $.each(v1.cohort,function(k2,v2){
                $.each(v2,function(k3,v3){
                    $.each(v3,function(k4,v4){
                        if (!(k2 in temp[k1].solved[k4])){
                            temp[k1].solved[k4][k2] = {'M':0,'F':0,'U':0};
                        }
                        temp[k1].solved[k4][k2][k3] = v4;
                    });
                });
            });
        }
    });
    return temp;
}
