import doctest
import inspect
import types
import collections
import random
import functools
import numpy as np

from interval import interval, inf

#from http://svn.colorstudy.com/home/ianb/ruby_blocks.py
def magic_set(obj):
    """
    Adds a function/method to an object.  Uses the name of the first
    argument as a hint about whether it is a method (``self``), class
    method (``cls`` or ``klass``), or static method (anything else).
    Works on both instances and classes.

        >>> class color:
        ...     def __init__(self, r, g, b):
        ...         self.r, self.g, self.b = r, g, b
        >>> c = color(0, 1, 0)
        >>> c      # doctest: +ELLIPSIS
        <__main__.color instance at ...>
        >>> @magic_set(color)
        ... def __repr__(self):
        ...     return '<color %s %s %s>' % (self.r, self.g, self.b)
        >>> c
        <color 0 1 0>
        >>> @magic_set(color)
        ... def red(cls):
        ...     return cls(1, 0, 0)
        >>> color.red()
        <color 1 0 0>
        >>> c.red()
        <color 1 0 0>
        >>> @magic_set(color)
        ... def name():
        ...     return 'color'
        >>> color.name()
        'color'
        >>> @magic_set(c)
        ... def name(self):
        ...     return 'red'
        >>> c.name()
        'red'
        >>> @magic_set(c)
        ... def name(cls):
        ...     return cls.__name__
        >>> c.name()
        'color'
        >>> @magic_set(c)
        ... def pr(obj):
        ...     print obj
        >>> c.pr(1)
        1
    """
    def decorator(func):
        is_class = (isinstance(obj, type)
                    or isinstance(obj, types.ClassType))
        args, varargs, varkw, defaults = inspect.getargspec(func)
        if not args or args[0] not in ('self', 'cls', 'klass'):
            # Static function/method
            if is_class:
                replacement = staticmethod(func)
            else:
                replacement = func
        elif args[0] == 'self':
            if is_class:
                replacement = func
            else:
                def replacement(*args, **kw):
                    return func(obj, *args, **kw)
                try:
                    replacement.func_name = func.func_name
                except:
                    pass
        else:
            if is_class:
                replacement = classmethod(func)
            else:
                def replacement(*args, **kw):
                    return func(obj.__class__, *args, **kw)
                try:
                    replacement.func_name = func.func_name
                except:
                    pass
        setattr(obj, func.func_name, replacement)
        return replacement
    return decorator

class memoized(object):
    """Decorator that caches a function's return value each time it is called.
    If called later with the same arguments, the cached value is returned, and
    not re-evaluated.
    """
    def __init__(self, func):
        self.func = func
        self.cache = {}
    def __call__(self, *args):
        try:
            return self.cache[args]
        except KeyError:
            value = self.func(*args)
            self.cache[args] = value
            return value
        except TypeError:
            # uncachable -- for instance, passing a list as an argument.
            # Better to not cache than to blow up entirely.
            return self.func(*args)
    def __repr__(self):
        """Return the function's docstring."""
        return self.func.__doc__
    def __get__(self, obj, objtype):
        """Support instance methods."""
        return functools.partial(self.__call__, obj)

def range_dist(m, range):
    '''given a value 'm' and a range of values as a tuple 'range', give the
       absolute distance of m from the middle of this range.
    '''
    r_mean = sum(range) / 2
    return abs(m - r_mean)

def mod_str(s, slice, repl):
    '''given a slice of the string to modify, the original string,
       and a replacement, return a new string with the area in the slice region
       replaced
    '''
    return s[:slice.start] + repl + s[slice.stop:]

def find_all(str, substr):
    '''given a string and a substring, return all occurrences of the substring
       inside the string as a tuple of ints, corresponding to the starting
       position of the substr
    '''
    starts = ()
    last_start = str.find(substr)
    while last_start != -1:
        starts += (last_start,)
        last_start = str.find(substr, last_start + 1)
    return starts

def generate_nmers(seq, bounds):
    '''generate a list of n-mers at each position in the seq (a string), relative
    to each position i, and remove any substrings that are too short
    '''

    #if one bound, start at 0 rel to i
    if len(bounds) == 1:
        bounds = (0,) + bounds
    elif len(bounds) != 2:
        raise BoundsException

    motif_len = bounds[1] - bounds[0]
    str_list = []
    loc_list = []

    for i in range(len(seq)):

        i_bounds = (bounds[0] + i, bounds[1] + i)
        i_slice = slice(*i_bounds)
        #generate n-mers from bounds relative to position i
        if len(seq[i_slice]) >= motif_len:
            str_list.append(str(seq[i_slice]))
            loc_list.append(i_bounds)

    return (str_list, loc_list)

def to_raw(plain_str):
    '''
    This simple function converts a string into a raw string literal so that,
    for instance, regular expressions do not mistake a literal '.' for a 'match
    anything' character.
    '''

    escape_dict = {
        '.' : '\.',
        '(' : '\(',
        ')' : '\)',
    }

    temp_str = "%r" % plain_str
    new_raw_str = temp_str[1:-1]
    return "".join([escape_dict.get(char,char) for char in new_raw_str])

def iter_len(iterable):
    return sum(1 for _ in iterable)

def max_none(iterable):
    try:
        return max(iterable)
    except ValueError:
        return None

def mean_none(iterable):
    try:
        return np.mean(iterable)
    except FloatingPointError:
        return None

@magic_set(interval)
def invert(self):
    '''
    separates all the components and takes the intersection of their individual
    inverses
    '''

    ivls = []

    inverse_func = lambda cmp: interval([-inf, cmp[0][0]], [cmp[0][1], inf])
    inverses = map(inverse_func, self.components)
    return reduce(lambda x, y: x & y, inverses, interval((-inf, inf)))

@magic_set(interval)
def to_tuple(self):
    '''
    return start and end of the first (and hopefully only component) to ivl
    '''
    return (int(self[0][0]), int(self[0][1]))

def irandomize(iterable, bufsize=1000, seed=1):
    random.seed(seed)
    #print "Here's your random seed:\t{}".format(seed)
    iterable = iter(iterable)
    items = []
    try:
        while True:
            for i in xrange(random.randint(1, bufsize)):
                items.append(iterable.next())
            random.shuffle(items)
            for i in xrange(random.randint(1, bufsize)):
                if items:
                    yield items.pop()
                else:
                    break
    except StopIteration:
        random.shuffle(items)
        for item in items:
            yield item
        raise StopIteration

@magic_set(interval)
def overlaps(self, ivl2):
    '''
    A & B returns true if there is one point the same, i,e between
    (a,b) and (b,c), but we only want to return true if these intervals have
    some non-empty space in common, not just one point.
    '''
    isct = self & ivl2
    return bool((isct and len(isct) >= 1 and (isct[0][1] - isct[0][0]) > 0))

def str_diff(str1, str2, case_sensitive=True):
    '''
    for two strings of the same length, tells you positions that are different
    case sensitive by default!
    '''

    if not case_sensitive:
        str1 = str1.upper()
        str2 = str2.upper()

    diffs = ()

    if len(str1) != len(str2): return ValueError
    if str1 == str2: return ()
    for i in range(len(str1)):
        if str1[i] != str2[i]: diffs += (i,)
    return diffs

def shuffle_and_return(sequence):
    random.shuffle(sequence)
    return sequence

@magic_set(interval)
def sum_len(self):
    ''' this gives the summed length of every component of the interval. '''
    ivl_len = 0
    for comp in self:
        ivl_len += comp[1] - comp[0]
    return ivl_len

