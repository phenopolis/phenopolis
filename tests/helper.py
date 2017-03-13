

def login(app):
    return app.post('/login', data=dict(
        name='demo',
        password='demo123'
    ), follow_redirects=True)
