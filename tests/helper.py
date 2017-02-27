

def login(app):
    return app.post('/login', data=dict(
        username='demo',
        password='demo123'
    ), follow_redirects=True)
