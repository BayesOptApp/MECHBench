### THIS IS A MODULE TO SET-UP THE CASE FROM THE START

# Module Properties
__author__ = "Ivan Olarte Rodriguez"


# Module Imports
import warnings
from collections import namedtuple, OrderedDict



def setup_default_options_(
        open_radioss_main_path='/home/ivanolar/Documents/OpenRadioss2/OpenRadioss_linux64/OpenRadioss/Tools/openradioss_gui'\
                        '# or any path which may exist',
        write_vtk='0 # or 1 if you need to generate the vtk files of the simulation',
        h_level = '1 # Set integers given a mesh refinement level',
        nt = '1 # Set integers as number of threads in the case',
        np = '1 # Set integers as number of processors in the case',
        gmsh_verbosity = '0 # Set integers as verbosity level of gmsh',
        save_mesh_vtk = '0 # or 1 if you want GMSH/mesher to save an extra file with the generated mesh'+\
                        'as a vtk file'
        
        
):
    r"""
    This is just a function defined in order to define default options on how to set-up cases in
    the OpenRadioss Environment.
    """

    return dict(locals())


_default_options:dict = setup_default_options_() # This will be later reassigned
_allowed_options_keys = dict([s.lower(), s] for s in _default_options) # Get the keys from the defaults

def safe_str(s, known_words:dict=None):
    """return ``s`` as `str` safe to `eval` or raise an exception.

    Strings in the `dict` `known_words` are replaced by their values
    surrounded with a space, which the caller considers safe to evaluate
    with `eval` afterwards.

    """
    safe_chars = ' 0123456789.,+-*()[]e<>=/\\'
    if s != str(s):
        return str(s)
    if not known_words:
        known_words = {}
    stest = s[:]  # test this string
    sret = s[:]  # return this string
    for word in sorted(known_words.keys(), key=len, reverse=True):
        stest = stest.replace(word, '  ')
        sret = sret.replace(word, " %s " % known_words[word])
    for c in stest:
        if c not in safe_chars:
            raise ValueError('"%s" is not a safe string'
                             ' (known words are %s)' % (s, str(known_words)))
    return sret

class RunnerOptions(dict):
    r"""
    This is a dictionary with the available options and default values to
    run each of the test-cases/benchmarks included in this repository
    """

    @staticmethod
    def defaults()->dict:
        """return a dictionary with default option values and description"""
        return _default_options
    
    def check(self, options=None):
        """check for ambiguous keys and move attributes into dict"""
        self.check_values(options)
        self.check_attributes(options)
        self.check_values(options)
        return self
    
    def check_values(self, options=None):
        corrected_key = RunnerOptions().corrected_key  # caveat: infinite recursion
        validated_keys = []
        original_keys = []
        if options is None:
            options = self
        for key in options:
            correct_key = corrected_key(key)
            if correct_key is None:
                raise ValueError('%s is not a valid option.\n'
                                 'Similar valid options are %s\n'
                                 'Valid options are %s' %
                                (key, str(list(_default_options(key))),
                                 str(list(_default_options))))
            if correct_key in validated_keys:
                if key == correct_key:
                    key = original_keys[validated_keys.index(key)]
                raise ValueError("%s was not a unique key for %s option"
                    % (key, correct_key))
            validated_keys.append(correct_key)
            original_keys.append(key)
        return options
    
    def check_attributes(self, opts=None):
        """check for attributes and moves them into the dictionary"""
        if opts is None:
            opts = self
        if 11 < 3:
            if hasattr(opts, '__dict__'):
                for key in opts.__dict__:
                    if key not in self._attributes:
                        raise ValueError("""
                        Assign options with ``opts['%s']``
                        instead of ``opts.%s``
                        """ % (opts.__dict__.keys()[0],
                               opts.__dict__.keys()[0]))
            return self
        else:
        # the problem with merge is that ``opts['ftarget'] = new_value``
        # would be overwritten by the old ``opts.ftarget``.
        # The solution here is to empty opts.__dict__ after the merge
            if hasattr(opts, '__dict__'):
                for key in list(opts.__dict__):
                    if key in self._attributes:
                        continue
                    warnings.warn(
                        """
        An option attribute has been merged into the dictionary,
        thereby possibly overwriting the dictionary value, and the
        attribute has been removed. Assign options with

            ``opts['%s'] = value``  # dictionary assignment

        or use

            ``opts.set('%s', value)  # here isinstance(opts, RunnerOptions)

        instead of

            ``opts.%s = value``  # attribute assignment
                        """ % (key, key, key), 'check', 'RunnerOptions')

                    opts[key] = opts.__dict__[key]  # getattr(opts, key)
                    delattr(opts, key)  
            return opts
    

    def __init__(self, s=None, **kwargs):
        """return a `RunnerOptions` instance.

        Return default options if ``s is None and not kwargs``,
        or all options whose name or description contains `s`, if
        `s` is a (search) string (case is disregarded in the match),
        or with entries from dictionary `s` as options,
        or with kwargs as options if ``s is None``,
        in any of the latter cases not complemented with default options
        or settings.

        Returns: see above.

        Details: as several options start with ``'s'``, ``s=value`` is
        not valid as an option setting.

        """
        # if not RunnerOptions.defaults:  # this is different from self.defaults!!!
        #     RunnerOptions.defaults = fmin([],[])
        if s is None and not kwargs:
            super(RunnerOptions, self).__init__(RunnerOptions.defaults())  # dict.__init__(self, RunnerOptions.defaults()) should be the same
            # self = RunnerOptions.defaults()
            s = 'nocheck'
        elif isinstance(s,str) and not s.startswith('unchecked'):
            super(RunnerOptions, self).__init__(RunnerOptions().match(s))
            # we could return here
            s = 'nocheck'
        elif isinstance(s, dict):
            if kwargs:
                raise ValueError('Dictionary argument must be the only argument')
            super(RunnerOptions, self).__init__(s)
        elif kwargs and (s is None or s.startswith('unchecked')):
            super(RunnerOptions, self).__init__(kwargs)
        else:
            raise ValueError('The first argument must be a string or a dict or a keyword argument or `None`')
        if not isinstance(s,str) or not s.startswith(('unchecked', 'nocheck')):
            # was main offender
            self.check()  # caveat: infinite recursion
            for key in list(self.keys()):
                correct_key = self.corrected_key(key)
                if correct_key not in RunnerOptions.defaults():
                    warnings.warn('invalid key ``' + str(key) +
                                   '`` removed', '__init__', 'RunnerOptions')
                    self.pop(key)
                elif key != correct_key:
                    self[correct_key] = self.pop(key)
        # self.evaluated = False  # would become an option entry
        self._lock_setting = False
        self._attributes = self.__dict__.copy()  # are not valid keys
        self._attributes['_attributes'] = len(self._attributes)

    
    def corrected_key(self, key):
        """return the matching valid key, if ``key.lower()`` is a unique
        starting sequence to identify the valid key, ``else None``

        """
        if key.startswith('_'):
            return key
        matching_keys = []
        key = key.lower()  # this was somewhat slow, so it is speed optimized now
        if key in _allowed_options_keys:
            return _allowed_options_keys[key]
        for allowed_key in _allowed_options_keys:
            if allowed_key.startswith(key):
                if len(matching_keys) > 0:
                    return None
                matching_keys.append(allowed_key)
        return matching_keys[0] if len(matching_keys) == 1 else None
    
    def complement(self):
        """add all missing options with their default values"""

        # add meta-parameters, given options have priority
        self.check()
        for key in RunnerOptions.defaults():
            if key not in self:
                self[key] = RunnerOptions.defaults()[key]
        return self
    
    def eval(self, key, default=None, loc=None, correct_key=True):
        """Evaluates and sets the specified option value in
        environment `loc`. Many options need ``N`` to be defined in
        `loc`, some need `popsize`.

        Details
        -------
        Keys that contain 'filename' are not evaluated.
        For `loc` is None, the self-dict is used as environment

        :See: `evalall()`, `__call__`

        """
        # TODO: try: loc['dim'] = loc['N'] etc
        if correct_key:
            # in_key = key  # for debugging only
            key = self.corrected_key(key)
        self[key] = self(key, default, loc)
        return self[key]

    def evalall(self, loc=None, defaults=None):
        """Evaluates all option values in environment `loc`.

        :See: `eval()`

        """
        self.check()
        if defaults is None:
            defaults = _default_options
        # TODO: this needs rather the parameter N instead of loc
        # if 'N' in loc:  # TODO: __init__ of CMA can be simplified
        #     popsize = self('popsize', defaults['popsize'], loc)
        #     for k in list(self.keys()):
        #         k = self.corrected_key(k)
        #         if k.startswith('_'):
        #             continue
        #         self.eval(k, defaults[k],
        #                   {'N':loc['N'], 'popsize':popsize})
        self._lock_setting = True
        return self

    def match(self, s=''):
        """return all options that match, in the name or the description,
        with string `s`, case is disregarded.

        """
        match = s.lower()
        res = {}
        for k in sorted(self):
            s = str(k) + '=\'' + str(self[k]) + '\''
            if match in s.lower():
                res[k] = self[k]
        return RunnerOptions(res)
    
    def amend_integer_options(self, inopts:dict):
        if inopts.get('h_level') in (None, _default_options['h_level']):
            self['h_level'] = 1
        else:
            self['h_level'] = int(inopts.get('h_level'))
        
        if inopts.get('nt') in (None, _default_options['nt']):
            self['nt'] = 1
        else:
            self['nt'] = int(inopts.get('nt'))
        
        if inopts.get('np') in (None, _default_options['np']):
            self['np'] = 1
        else:
            self['np'] = int(inopts.get('np'))

    def __call__(self, key, default=None, loc=None):
        """evaluate and return the value of option `key` on the fly, or
        return those options whose name or description contains `key`,
        case disregarded.

        Details
        -------
        Keys that contain `filename` are not evaluated.
        For ``loc==None``, `self` is used as environment
        but this does not define ``N``.

        :See: `eval()`, `evalall()`

        """
        global_env = dict(globals())

        try:
            val = self[key]
        except Exception:
            return self.match(key)

        if loc is None:
            loc = self  # TODO: this hack is not so useful: popsize could be there, but N is missing
        try:
            if isinstance(val,str):
                val = val.split('#')[0].strip()  # remove comments
                if (key.find('filename') or key.find('path')) < 0 and not (key == 'seed' and val.startswith('time')):
                        # and key.find('mindx') < 0:
                    val = eval(safe_str(val), global_env, loc)
            # invoke default
            # TODO: val in ... fails with array type, because it is applied element wise!
            # elif val in (None,(),[],{}) and default is not None:
            elif val is None and default is not None:
                val = eval(safe_str(default), global_env, loc)
        except Exception as e:
            if not str(e).startswith('"initial"'):
                warnings.warn(str(e))
            pass  # slighly optimistic: the previous is bug-free
        return val
    
    def from_namedtuple(self, t:OrderedDict):
        """update options from a `collections.namedtuple`.
        :See also: `to_namedtuple`
        """
        return self.update(t._asdict())
    
    def pprint(self, linebreak=80):
        for i in sorted(self.items()):
            s = str(i[0]) + "='" + str(i[1]) + "'"
            a = s.split(' ')

            # print s in chunks
            line = ''  # start entire to the left
            while a:
                while a and len(line) + len(a[0]) < linebreak:
                    line += ' ' + a.pop(0)
                print(line)
                line = '        '  # tab for subsequent lines


default_options = RunnerOptions(_default_options)