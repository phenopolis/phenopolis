from views import *

@app.route('/pubmedbatch/', methods=['GET', 'POST'])
@requires_auth
def pubmedbatch_main():
    # this is the main page of pubmedBatch
    # It allows a user to create a "folder" to hold results in the mangodb. Post request will handle the 'folder' creation
    '''
    db.results = {
            'user_id':'Jing',
            'folder':{
                foo':{
                    'test':{JSON}
                    'test(1):{JSON}
            },
        },
    }
    '''

    user = session.get('user')
    db = get_db('pubmedbatch')
    user_db = db.results.find_one({'user_id':user},{'_id':0})
    print 'user_id:%s:' % user_db
    if request.method == 'POST':
        # time to create a folder
        # first check if the folder already exists. Yes? pass. No? create
        folder = request.form['create-folder']
        
        user_folders =  [d for d in user_db['folder']]
        if folder not in user_folders:
            # set dot.notation
            user_folder = 'folder.' + folder
            db.results.update({'user_id':user}, {'$set':{user_folder:{} }})
    else:
        # get folders. If user not exists in pubmed db, create it
        #print(res)
        if not user_db:
            # A new user dropping by...
            print('A "user" is being made in pubmedDB!')
            db.results.insert({'user_id':user, 'folder':{}})
    # let's render the template  
    user_db = db.results.find_one({'user_id':user},{'_id':0})
    print user_db
    return render_template( 'pubmedbatch_main.html', 
            user = user,
            folders = [d for d in user_db['folder']]
    )

'''
delete a folder in pubmedbatch
'''
@app.route('/pubmedbatch_folderdel/<folder>', methods=['POST'])
@requires_auth
def pubmedbatch_delfolder(folder):
    user = session.get('user')
    db = get_db('pubmedbatch')
    user_db = db.results.find_one({'user_id':user},{'_id':0})
    folders = [d for d in user_db['folder']]
    if folder in folders:
        folder = 'folder.' + folder
        db.results.update({'user_id':user}, {'$unset': {folder:''}})
        return redirect('/pubmedbatch/')
    else:
        return 'Folder: <b>%s</b> does not exist!' % folder

'''
main business
'''
@app.route('/pubmedbatch/<folder>', methods=['GET', 'POST'])
@requires_auth
def pubmedbatch(folder):
    # This is the main business
    user = session.get('user') or app.config['DEFAULT_USER']
    db = get_db('pubmedbatch')
    if request.method == 'POST':
        # post. get form data, return JSON
        ##########################################################
        # get form data
        #########################################################
        column = int(request.form['column']) - 1
        Entrez.email = request.form['email']
        AND_term = request.form['AND']
        OR_term = request.form['OR']
        verbose = request.form.get('verbose','') 
        # read RetNet file
        retnet_f = 'retnet.json'
        RETNET = json.load(open(retnet_f, 'r')) 
        # set session
        session['AND'] = AND_term
        session['OR'] = OR_term
        session['EMAIL'] = request.form['email']
        session['COLUMN'] = column + 1 
        csv_file = request.files['csv_upload']
        file_name = csv_file.filename
        known_genes = request.files['known_genes'] or ''
        mask_genes = request.files['mask_genes'] or ''
        #########################################################
        # parse known genes and mask genes
        known_genes = known_genes.read().split() if known_genes else ''
        mask_genes = mask_genes.read().split() if mask_genes else ''
        #########################################################
        # read csv. has to make the csv string looking like a filehandle
        csv_content = csv_file.read() 
        '''
        progress bar initiation
        '''
        # how many lines? minus header
        csv_lines = len(re.findall(r'\n', csv_content)) - 1
        # a random id
        progress_id = user + request.form['rid']
        bad = init_progress_bar(progress_id, csv_lines)
        if bad:
            print 'weird, progress id conflict!'
        '''
        '''
        csvreader = csv.reader(StringIO.StringIO(csv_content), delimiter=',', quotechar='"')
        # life time from config
        life = app.config['PUBMEDBATCH_LIFE']
        # number of lines?
        line_num = len(re.findall('\n', csv_content))
        # format terms
        OR = OR_term.split()
        OR.sort()
        # make reg term for later use (calculating pubmed score
        reg = '\\b|\\b'.join(OR)
        reg = '\\b' + reg + '\\b'
        reg = re.compile(reg, re.IGNORECASE) 
        AND = AND_term.split()
        AND.sort()
        smashed_OR = ['"' + t + '"' + '[Title/Abstract]' for t in OR]
        smashed_AND = ['"' + t + '"' + '[Title/Abstract]' for t in AND]
        smashed_OR = ' OR '.join(smashed_OR)
        smashed_AND = ' AND '.join(smashed_AND)
        smashed_term = ' AND (' + smashed_OR + ')'
        if smashed_AND:
            smashed_term += ' AND ' + smashed_AND
        ###########################################################
        # it's time to read the csv file
        ###########################################################
        row_num = -1
        header = [] # header
        output = [] # all the results get pushed here
        for row in csvreader:
            row_num += 1
            # read header
            if not row_num:
                header = row
                # replace . with _ in header
                header = [re.sub('\.','_', h) for h in header]
                # add 2 columns after HUGO
                header[column+1:column+1] = ['ref(pubmedID)', 'pubmed_score']
                # add a pred score at the beginning
                header[0:0] = ['pred_score']
                continue
            # read in real data
            gene_name = row[column]
            if gene_name == 'NA':
                continue
            # get rid of any parentheses and their content
            gene_name = re.sub(r'\([^)]*\)?','',gene_name) 
            # update progress bar
            update_progress_bar({'id':progress_id, 'message': '%s;' % gene_name}) 
            genes={} # storing 
            print gene_name
            smashed_all = gene_name + smashed_term
            ####################################################
            # first to see if masked, pass
            # else
            #   db.cache
            #       when in db.cache and up-to-date, get the pubmed result
            #       when in db.cache and out-of-date, add new pubmed result
            #       when not in db.cache, search pubmed and store
            ####################################################
            # db.cache's structure:
            #  {['key': '_'.join([gene_name.upper(), ','.join(OR).lower(), ','.join(AND).lower()]), 
            #   'result': {pubmed_result},
            #   'date': now]}
            if gene_name in mask_genes:
                # masked. don't search, just return the row
                if not verbose:
                    continue
                row[column+1:column+1] = [{'total_score': 0, 'results': ['masked']}, 0]
                # add a placeholder for pred_score
                row[0:0] = [0]
                ha = {}
                for col in range(len(header)):
                    ha[header[col]] = row[col]
                ha['pred_score'] = get_pred_score(ha)
                output.append(ha)
            else:
                # not masked
                now = time.mktime(time.localtime()) #get time in seconds
                lag = 0 # use it as a flag of how to search. 0 = search; now-saved['date'] = update; 
                term = '_'.join([gene_name.upper(), ','.join(OR).lower(), ','.join(AND).lower()])
                # now check if the result in the db is uptodate
                saved = db.cache.find_one({'key': term})
                if saved:
                    lag = now - saved['date']
                    # has record. let's see if it is out of date
                    if lag  <= life:
                        # not out of date. let's use it
                        genes[gene_name] = saved['result'] 
                if gene_name not in genes:
                    # need to search
                    # handle = Entrez.einfo()
                    # record = Entrez.read(handle)
                    # handle.close()
                    results = scrutinise({'lag':lag, 'smashed_all': smashed_all, 'reg':reg}) 
                    # update the database, and maybe results
                    if lag:
                        # need to update the database
                        # old result = saved
                        results['total_score'] = results['total_score'] + saved['result']['total_score']
                        results['results'].extend(saved['result']['results'])
                        # update database now
                        db.cache.update({'key': term}, {'$set': {'result': results, 'date': now}})
                    else:
                        # write to the database
                        db.cache.insert({'key': term, 'result': results, 'date': now}) 
                    genes[gene_name] = results 
                # verbose and no score?
                if not (verbose or genes[gene_name]['total_score']):
                    continue
                else:
                    # known genes?
                    if gene_name in known_genes:
                        genes[gene_name]['known'] = 1
                        # give the rest a minimal score to keep them on the list
                        genes[gene_name]['total_score'] = max(1, genes[gene_name]['total_score'])
                    else:
                        genes[gene_name]['known'] = 0
                    # give the rest a minimal score to keep them on the list
                    genes[gene_name]['total_score'] = max(1, genes[gene_name]['total_score'])
                    # Retnet?
                    if gene_name in RETNET:
                        genes[gene_name]['disease'] = RETNET[gene_name]['disease']
                        genes[gene_name]['omim'] = RETNET[gene_name]['omim']
                        genes[gene_name]['mode'] = RETNET[gene_name]['mode']
                        # reassign total score according to mode
                        if RETNET[gene_name]['mode'] == 'd' or RETNET[gene_name]['mode'] == 'x':
                            genes[gene_name]['total_score'] = max(100, genes[gene_name]['total_score'])
                        elif re.search(r'/(?!dominant)/', AND_term):
                            # not searching doimant, also assign others to 100
                            genes[gene_name]['total_score'] = max(100, genes[gene_name]['total_score'])
                        else:
                            # give the rest a minimal score to keep them on the list
                            genes[gene_name]['total_score'] = max(1, genes[gene_name]['total_score']) 
                    # add pubmed result
                    row[column+1:column+1] = [genes[gene_name], genes[gene_name]['total_score']]
                    # add a placeholder for pred_score
                    row[0:0] = [0]
                    ha = {}
                    for col in range(len(header)):
                        ha[header[col]] = row[col]
                    ha['pred_score'] = get_pred_score(ha)
                    output.append(ha) 
        # save results to folder
        # if the file name already exists?
        #folders = db.results.find_one({'user_id':user},{'_id':0})['folder']
        the_folder = db.results.find_one({
            'user_id': user},{'_id': 0})['folder'][folder]
        files = [ f for f in the_folder]
        file_name = os.path.splitext(file_name)[0]
        num = 0
        while file_name in files:
            num += 1
            file_name = re.sub(r'(\(\d+\))?$', '(%s)' % num, file_name) 
        # get the search term, to display
        search_term = 'AND[ <b>%s</b> ]; OR[ <b>%s</b> ]' % (', '.join(AND), ', '.join(OR))
        # write to database and win
        final_result = [header, output, search_term]
        the_file = 'folder.' + folder + '.' + file_name
        db.results.update({'user_id': user}, {
            "$set": {the_file: final_result}
        })
        return json.dumps([header, output, search_term, file_name])
    else:
        # get. display page
        # First see if folder exists. if not, return error
        user_folders = db.results.find_one({
            'user_id': user},{'_id': 0})['folder']
        # get AND an OR field
        AND = session.get('AND') or ''
        OR = session.get('OR') or app.config['PUBMEDBATCH_OR']
        EMAIL = session.get('EMAIL') or app.config['PUBMED_EMAIL']
        COLUMN = session.get('COLUMN') or '4'
        #user_folders = user_folders['folder']
        if folder not in user_folders:
            return "Error: " + folder + " does not exist!"
        # get the files in the folder to display
        files = [e for e in user_folders[folder]]
        return render_template('pubmedbatch.html',
            home_pubmedbatch = home_pubmedbatch,
            files = files,
            user = user, 
            folder = folder,
            AND = AND,
            OR = OR,
            EMAIL = EMAIL,
            COLUMN = COLUMN
        )



'''
Rename a file in pubmedbatch
'''
@app.route('/pubmedbatch_rename/<folder>', methods=['POST'])
@requires_auth
def pubmedbatch_rename(folder):
    user = session.get('user') 
    new_name = request.form['new-name']
    old_name = request.form['old-name'] 
    # check if new-name already exists
    db = get_db('pubmedbatch')
    the_folder = db.results.find_one({'user_id':user}, {'_id':0})['folder'][folder]
    files = [ f for f in the_folder]
    # sanity check
    if new_name in files:
        # throw an error
        return '<span style="color:orange">Oops, the new name: <b>%s</b> already exists in the folder</span>' % new_name
    if old_name not in files: # shouldn't happen though
        return '<span style="color:orange">Oops, the old name: <b>%s</b> does not exist in the folder, somehow</span>' % old_name
    # rename
    # mongodb doesn't like multiple '$' for the time being, so need to find hard index
    old_name = 'folder.' + folder + '.' + old_name
    new_name = 'folder.' + folder + '.' + new_name
    db.results.update({'user_id': user}, {'$rename':{old_name : new_name}}) 
    return 'The file has been successfully renamed!'


'''
Delete a file in pubmedbatch
'''
@app.route('/pubmedbatch_del/<path:path>', methods=['POST'])
@requires_auth
def pubmedbatch_del(path):
    user = session.get('user') 
    # get folder and file from route
    m = re.search(r'([^/]+)/([^/]+)', path)
    (folder_name, file_name) = m.groups()
    file_name = re.sub('%20', ' ', file_name) 
    db = get_db('pubmedbatch') 
    print db.results.update({'user_id':user}, {'$unset':{'folder.%s.%s' % (folder_name, file_name):''}}) 
    return '1'

'''
get a file in pubmedbatch
'''
@app.route('/pubmedbatch/<path:path>', methods=['POST'])
@requires_auth
def pubmedbatch_getfile(path):
    user = session.get('user')
    
    # get folder and file from route
    m = re.search(r'([^/]+)/([^/]+)', path)
    (folder_name, file_name) = m.groups()
    file_name = re.sub('%20', ' ', file_name)
    # get file from database
    db = get_db('pubmedbatch')
    file = db.results.find_one({'user_id': user})['folder'][folder_name][file_name]
    return json.dumps(file)


@app.route('/pubmedbatch_progress_bar/<id>')
def pubmedbatch_progress(id):
    '''
    progress bar query
    '''
    user = session.get('user') or app.config['DEFAULT_USER']
    progress_id = user + id
    return json.dumps(PROGRESS_BAR[progress_id])


