import bs4
import functools
from htmlmin.minify import html_minify

def prettify(route_function):
    @functools.wraps(route_function)
    def wrapped(*args, **kwargs):
        yielded_html = route_function(*args, **kwargs)
        soup = bs4.BeautifulSoup(yielded_html, 'lxml')
        return soup.prettify()

    return wrapped

def uglify(route_function):
    @functools.wraps(route_function)
    def wrapped(*args, **kwargs):
        yielded_html = route_function(*args, **kwargs)
        minified_html = html_minify(yielded_html)
        return minified_html

    return wrapped
