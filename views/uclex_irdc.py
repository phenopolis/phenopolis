from views import *
from lookups import *

@app.route('/irdc/summary')
def irdc_summary():
    if session['user'] == 'demo':
        return 'Error, you do not have access to this page'
    hpo_freq = get_hpo_size_freq('hpo_freq.tsv')
    hpo_dot_file = os.path.join('dot','irdc_hpo.json')
    hpo_dot_inf = open(hpo_dot_file,'r')
    hpo_dot = json.load(hpo_dot_inf)
    return render_template('irdc_summary.html', 
            hpo_dot = json.dumps(hpo_dot),
            hpo_freq = json.dumps(hpo_freq))

