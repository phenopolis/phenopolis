import jinja2
from jinja2.utils import soft_unicode

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

class FilterModule(object):
    ''' jinja2 filters '''
    def filters(self): return { 'map_format': map_format, 'get_attributes': get_attributes, }

jinja2.filters.FILTERS['map_format'] = map_format
jinja2.filters.FILTERS['get_attributes'] = get_attributes

