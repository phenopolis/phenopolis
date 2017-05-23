
# Uncomment to run this module directly.
#import sys, os
#sys.path.append(os.path.join(os.path.dirname(__file__), '..'))
# End of uncomment.

import load_data
import unittest
import runserver


class LoginTestCase(unittest.TestCase):

    def setUp(self):
        load_data.load_user_data()
        runserver.app.config['DB_NAME_USERS'] = 'test_users'
        runserver.app.config['TESTING'] = True
        self.app = runserver.app.test_client()

    def tearDown(self):
        pass

    # Smoke test to check that unittest is running ok.
    def test_smoke_test(self):
       pass 

    def login(self, username, password):
        return self.app.post('/login', data=dict(
            name=username,
            password=password
        ), follow_redirects=True)

    def logout(self):
        return self.app.get('/logout', follow_redirects=True)

    def change_password(self, username, password, new_pass_1):
        return self.app.post('/change_password', data=dict(
            change_pwd_name=username,
            current_password=password,
            new_password_1=new_pass_1,
        ), follow_redirects=True)

    def test_login_logout(self):
        rv = self.login('demox', 'demo123')
        assert rv.status_code == 401
        assert 'Invalid Credentials. Please try again.' in rv.data
        rv = self.login('demo', 'demo123x')
        assert rv.status_code == 401
        assert 'Invalid Credentials. Please try again' in rv.data
        rv = self.login('demo', 'demo123')
        assert rv.status_code == 200
        assert 'Authenticated' in rv.data
        rv = self.logout()
        assert rv.status_code == 200
        assert 'Please login' and 'username' and 'password' in rv.data

    def test_change_password(self):
        rv = self.login('test', 'test123')
        assert rv.status_code == 200
        assert 'Authenticated' in rv.data

        rv = self.change_password('test', 'test123', 'test456')
        assert rv.status_code == 200
        print(rv.data)
        assert 'Password for username \'test\' changed' in rv.data

        rv = self.login('test', 'test456')
        assert rv.status_code == 200
        assert 'Authenticated' in rv.data

        rv = self.login('test', 'test123')
        assert rv.status_code == 401

        rv = self.change_password('test', 'test456', 'test123')
        assert rv.status_code == 200

        rv = self.change_password('demo', 'demo123', 'demo456')
        assert rv.status_code == 401

        rv = self.change_password('test', 'x', 'test456')
        assert rv.status_code == 401

        rv = self.change_password('x', 'test123', 'test456')
        assert rv.status_code == 401

if __name__ == '__main__':
    unittest.main()
