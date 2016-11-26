import jinja2
from jinja2.utils import soft_unicode
from jinja2.utils import Markup

def get_attributes(value, *args):
    """
    """
    return (value[a] for a in args)


def map_format(value, pattern):
    """
    Apply python string formatting on an object:
        .. sourcecode:: jinja
    {{ "%s - %s"|format("Hello?", "Foo!") }}
        -> Hello? - Foo!
    """
    value=tuple((value))
    s=soft_unicode(pattern).format(*value)
    return s


def unique(value):
    """
    """
    return list(set((value)))


def href(value,link):
    return Markup(soft_unicode("<br>".join(["<a href=/"+link+"/"+x+" target=_blank>"+x+"</a>" for x in value])))


class FilterModule(object):
    ''' jinja2 filters '''
    def filters(self): return { 'map_format': map_format, 'get_attributes': get_attributes, 'unique':unique, 'href':href }

jinja2.filters.FILTERS['map_format'] = map_format
jinja2.filters.FILTERS['get_attributes'] = get_attributes
jinja2.filters.FILTERS['unique'] = unique
jinja2.filters.FILTERS['href'] = href

