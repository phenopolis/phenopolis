
from flask import jsonify
import lookups
from os import listdir, chdir
from os.path import isfile, join
import pymongo
from collections import defaultdict, Counter

class User(object):
    def __init__(self, user, user_db, groups=[], email='', affiliation=''):
        User.db=user_db
        data=User.db.users.find_one({'user':user},{'_id':False})
        if data:
            self.__dict__.update(data)
            self.status={ 'message':'User account exists already.'%user, 'http_code':401}
            return
        data=User.db.new_users.find_one({'user':user},{'_id':False})
        if data:
            self.__dict__.update(data)
            self.status={ 'message':'User account %s request already created, still unapproved.'%user, 'http_code':401}
            return
        self.user=user
        self.email=email
        self.affiliation=affiliation
        self.groups=groups
        self.status={ 'message':'User account request created for %s.'%user, 'http_code':200}
        # add to unapproved table
        User.db.new_users.ensure_index('user',unique=True)
        User.db.new_users.insert_one( self.__dict__ )
    def __getattribute__(self, key):
        "Emulate type_getattro() in Objects/typeobject.c"
        v = object.__getattribute__(self, key)
        if hasattr(v, '__get__'): return v.__get__(None, self)
        return v
    def json(self):
        if '_id' in self.__dict__: del self.__dict__['_id']
        return jsonify(result=self.__dict__)
    def save(self):
        print('writing', self.external_id, 'to database')
        return Patient.db.user.update({'user':self.user},self.__dict__,upsert=True)
    @property
    def password(self):
        pass
    @property
    def external_ids(self):
        pass
    @property
    def individuals(self):
        pass
    @property
    def approved(self):
        pass
