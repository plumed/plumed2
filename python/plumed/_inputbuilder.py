# /* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#   Copyright (c) 2011-2016 The plumed team
#   (see the PEOPLE file at the root of the distribution for a list of names)
#
#   See http://www.plumed.org for more information.
#
#   This file is part of plumed, version 2.
#
#   plumed is free software: you can redistribute it and/or modify
#   it under the terms of the GNU Lesser General Public License as published by
#   the Free Software Foundation, either version 3 of the License, or
#   (at your option) any later version.
#
#   plumed is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU Lesser General Public License for more details.
#
#   You should have received a copy of the GNU Lesser General Public License
#   along with plumed.  If not, see <http://www.gnu.org/licenses/>.
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
#
"""Input construction and metadata access through the public kernel interface."""

import array
import os
import re
import sys
import types

from ._plumed import Plumed


def _guessplumedroot(kernel=None):
    """Guess plumed root.

    kernel: path to the plumed kernel

    In case the Plumed object cannot be created, try to launch a `plumed` executable and obtain the root
    dir from there.
    """
    try:
        from tempfile import TemporaryDirectory

        # not a NamedTemporaryFile: on Windows the kernel could not open it
        # for writing while it is open here
        with TemporaryDirectory() as tmpdir:
            log = os.path.join(tmpdir, "plumed.log")
            with Plumed(kernel) as p:
                p.cmd("setLogFile", log)
                p.cmd("init")
            root = None
            # the kernel writes UTF-8 (paths included)
            with open(log, encoding="utf-8", errors="replace") as fin:
                for line in fin:
                    if re.match("PLUMED: Root: ", line):
                        root = re.sub("PLUMED: Root: ", "", line).rstrip("\n")
                        break
            if root and len(root) > 0:
                return root
    except Exception:
        # Ordinary lookup failures retain the legacy fallback; interrupts propagate.
        pass
    # alternative solution, search for a plumed executable in the path
    import subprocess

    return (
        subprocess.check_output(["plumed", "--no-mpi", "info", "--root"])
        .decode("utf-8")
        .rstrip()
    )


def _readkernelmetadata(kernel=None, instance=None):
    """Read live metadata from a selected kernel or an existing PLUMED instance.

    An existing instance takes precedence and is borrowed without finalizing
    or initializing it. Return None for older kernels so callers can use
    their installed Vim files.
    """
    import json
    from contextlib import nullcontext

    with nullcontext(instance) if instance is not None else Plumed(kernel) as p:
        api = array.array("i", [0])
        p.cmd("getApiVersion", api)
        if api[0] < 11:
            return None
        length = array.array("i", [0])
        p.cmd("getInputBuilderJSONLength", length)
        if length[0] <= 0:
            raise ValueError("PLUMED returned an invalid InputBuilder metadata length")
        buffer = array.array("i", [0]) * length[0]
        p.cmd("getInputBuilderJSON", buffer)
    # the kernel stores one byte per int: convert the values, not the raw memory of the array
    syntax = json.loads(array.array("B", buffer).tobytes().decode("utf-8"))
    dictionary = {}
    doc = {}
    for action, entry in syntax.items():
        if not isinstance(entry, dict) or not isinstance(entry.get("syntax"), dict):
            continue
        dictionary[action] = {}
        doc[action] = (
            "Create action "
            + action
            + "\n"
            + entry.get("help", entry.get("description", ""))
            + "\n"
        )
        doc[action] = doc[action].replace("\\f$", "").replace("\\", "")
        for name, keyword in entry["syntax"].items():
            # LABEL is supplied by InputBuilder itself. All other hidden
            # keywords remain usable, as they are in the legacy Vim dictionary.
            if name in ("LABEL", "output"):
                continue
            if not isinstance(keyword, dict) or "type" not in keyword:
                continue
            if "numbered" not in keyword:
                raise ValueError(
                    "PLUMED metadata lacks precise numbering for " + action + "." + name
                )
            numbered = keyword["numbered"]
            if isinstance(numbered, dict):
                numbered = numbered.get("enabled")
            if not isinstance(numbered, bool):
                raise ValueError(
                    "PLUMED metadata has invalid numbering for " + action + "." + name
                )
            if keyword["type"] == "flag":
                menu = "(flag)"
            elif numbered:
                menu = "(numbered)"
            else:
                menu = "(option)"
            dictionary[action][name] = menu
        # Match the Vim completion-word order, including the '=' of options,
        # so generated signatures and keyword handling retain their legacy order.
        dictionary[action] = dict(
            sorted(
                dictionary[action].items(),
                key=lambda item: item[0] + ("=" if item[1] == "(option)" else ""),
            )
        )
    return dictionary, doc


def _readvimdict(plumedroot=None, kernel=None):
    """Read VIM dictionary given the path to PLUMED root and return (dictionary,doc).

    If plumedroot is not given, it is guessed by launching `plumed` executable.
    The dictionary is structured as follows:
    - The keys are the names of the actions (e.g. "RESTRAINT").
    - The values are dictionaries structured as follows:
      - The keys are the available options.
      - The value is a string describing the option type.
    For instance `dictionary["RESTRAINT"]["NUMERICAL_DERIVATIVES"]=="(flag)"`.

    The doc is a dictionary structured as follows:
    - The keys are the names of the actions (e.g. "RESTRAINT").
    - The values are docstrings for the corresponding actions.
    """
    if plumedroot is None:
        plumedroot = _guessplumedroot(kernel)
    syntax_file = plumedroot + "/vim/syntax/plumed.vim"
    help_dir = plumedroot + "/vim/help"
    # local dictionary, read from VIM
    plumedDictionary = {}
    pattern = re.compile(r"^let b:plumedDictionary\[.*$")
    # the Vim files are UTF-8 (not the locale encoding, e.g. on Windows)
    with open(syntax_file, encoding="utf-8") as fin:
        for line in fin:
            if pattern.match(line):
                line = re.sub("^let b:", "", line)
                exec(
                    line, {"__builtins__": None}, {"plumedDictionary": plumedDictionary}
                )
    ret = {}
    doc = {}
    for action in plumedDictionary:
        ret[action] = {}
        doc[action] = {}
        # read documentation
        with open(help_dir + "/" + action + ".txt", encoding="utf-8") as fin:
            thisdoc = ""
            for line in fin:
                thisdoc += line
                if line.rstrip("\n") == "****************************************":
                    thisdoc = "Create action " + action + "\n"
            doc[action] = thisdoc
            # remove LaTex stuff
            doc[action] = doc[action].replace("\\f$", "").replace("\\", "")
        # read dictionary
        for opt in plumedDictionary[action]:
            # skip label (it is added automatically)
            if opt["menu"] != "(label)":
                ret[action][re.sub("=$", "", opt["word"])] = opt["menu"]
    return ret, doc


def _create_functions(dictionary, *, doc=None, append_underscores=False):
    """Create functions given dictionary and, optionally, documentation.

    Given a dictionary a doc produced with _readvimdict, it returns a dictionary
    containing the functions corresponding to each action. The functions are stored
    in string that should be then evaluated with exec. For each action (say, "RESTRAINT")
    we define both a function `def RESTRAINT and a docstring `RESTRAINT.__doc__`.
    These functions are only using builtins and functions from the _format_tools dictionary.
    """
    functions = {}
    for action in dictionary:
        # skip actions with incorrect name
        if re.match(".*-.*", action):
            continue
        fname = action
        if append_underscores:
            fname += "__"
        string = ""
        string += "def " + fname
        string += '(self,LABEL="",verbatim=None'
        if len(dictionary[action]) > 0:
            string += ",*"
        for w in dictionary[action]:
            # skip arguments with incorrect name
            if re.match(".*-.*", w):
                continue
            if dictionary[action][w] == "(flag)":
                string += "," + w + "=False"
            else:
                string += "," + w + "=None"
        string += ",**kwargs):\n"
        string += '  ret=""\n'
        string += "  ret+=_format_label(self,LABEL)\n"
        string += '  ret+="' + action + '"\n'
        string += "  retlist=[]\n"
        for w in dictionary[action]:
            # skip arguments with incorrect name
            if re.match(".*-.*", w):
                continue
            t = dictionary[action][w]
            if t == "(flag)":
                string += '  retlist.append(_format_flag(self,"' + w + '",' + w + "))\n"
            elif t == "(option)":
                string += '  retlist.append(_format_opt(self,"' + w + '",' + w + "))\n"
            elif t == "(numbered)":
                string += (
                    '  retlist.append(_format_numbered(self,"' + w + '",' + w + "))\n"
                )
            else:
                raise TypeError("error parsing dictionary, unknown type " + t)
        # now process kwargs to allow numbered arguments
        string += "  for arg in kwargs:\n"
        string += "      import re\n"
        string += "      allowed=[]\n"
        for x in dictionary[action]:
            string += '      allowed.append("' + x + '")\n'
        string += '      if not re.sub("[0-9]*$","",arg) in allowed:\n'
        string += '         raise TypeError("unknown arg " + arg)\n'
        string += "      retlist.append(_format_anything(self,arg,kwargs[arg]))\n"
        # sorting is necessary to make the line reproducible in regtests
        string += "  retlist.sort()\n"
        string += "  for x in retlist:\n"
        string += "      if(len(x)>0):\n"
        string += "          ret+=x\n"
        string += "  ret+=_format_verbatim(self,verbatim)\n"
        string += "  return _format_return(self,ret," + repr(action) + ")\n"
        # set docstring
        if doc is not None:
            string += fname + ".__doc__ = " + repr("\n" + doc[action]) + "\n"
        functions[fname] = string
    return functions


# formatting tools


class _numbered:
    def __init__(self, *args):
        # this is experimental:
        # it allows calling _numbered(arg1,arg2,arg3)
        # it is however on purpose not implemented in the numbered() method below
        if len(args) == 1:
            if isinstance(args[0], dict):
                self.arg = args[0]
            elif hasattr(args[0], "__iter__"):
                self.arg = {}
                i = 0
                for x in args[0]:
                    self.arg[i] = x
                    i += 1
            else:
                raise TypeError(
                    "when calling numbered with 1 argument, it should be a list/tuple/dictionary"
                )
        else:
            self.arg = {}
            i = 0
            for x in args:
                self.arg[i] = x
                i += 1


class _replicas:
    def __init__(self, *args):
        # this is experimental:
        # it allows calling _replicas(arg1,arg2,arg3)
        # it is however on purpose not implemented in the replicas() method below
        if len(args) == 1:
            if hasattr(args[0], "__iter__"):
                self.arg = args[0]
            else:
                raise TypeError(
                    "when calling replicas with 1 argument, it should be a list/tuple"
                )
        else:
            self.arg = args


## tool to format at strings (optional)


def _format_at_one_residue(builder, name, residue, chain):
    if isinstance(chain, int):
        chain = str(chain)
    digit = False
    for i in chain:
        if i.isdigit():
            digit = True
    if digit:
        chain = chain + "_"
    if isinstance(residue, int) or isinstance(residue, str):
        return "@" + name + "-" + chain + str(residue)
    else:
        assert False


def _format_at_one_chain(builder, name, residue, chain):
    res = ""
    if hasattr(residue, "__iter__") and not isinstance(residue, str):
        for x in residue:
            res += builder._separator + _format_at_one_residue(builder, name, x, chain)
    else:
        res += builder._separator + _format_at_one_residue(
            builder, name, residue, chain
        )

    return res


def _format_at(builder, name, residue, chain=""):
    res = ""
    if hasattr(chain, "__iter__") and not isinstance(chain, str):
        for x in chain:
            res += _format_at_one_chain(builder, name, residue, x)
    else:
        res += _format_at_one_chain(builder, name, residue, chain)
    return res[len(builder._separator) :]


class _at:
    def __init__(self, builder):
        import weakref

        self._builder = weakref.ref(builder)
        _at_global = [
            "mdatoms",
            "allatoms",
            "water",
            "nucleic",
            "protein",
            "water",
            "ions",
            "hydrogens",
            "nonhydrogens",
        ]
        for x in _at_global:
            exec("self." + x + '="@' + x + '"', None, {"self": self})
        _at_residue = [
            "phi",
            "psi",
            "omega",
            "chi1",
            "alpha",
            "beta",
            "gamma",
            "delta",
            "epsilon",
            "zeta",
            "v0",
            "v1",
            "v2",
            "v3",
            "v4",
            "chi",
            "back",
            "sugar",
            "base",
            "lcs",
        ]
        for x in _at_residue:
            ldict = {}
            exec(
                "def "
                + x
                + '(self,residue,chain=""):\n          return _format_at(self._builder(),"'
                + x
                + '",residue,chain)\n',
                {"_format_at": _format_at},
                ldict,
            )
            exec(
                "self." + x + " = types.MethodType( ldict['" + x + "'],self)",
                None,
                {"self": self, "types": types, "ldict": ldict},
            )

    def __call__(self, name, residue, chain=""):
        return _format_at(self._builder(), name, residue, chain)


## end of tool to format at strings


def _fix_braces(builder, arg, comma):
    """Fix braces.

    If comma is true, consider comma as a separator. Otherwise, only consider space as a separator.
    """
    always = not builder._minimize_braces
    # recursively remove braces to find non-matching ones
    tmp = arg
    go = 1
    while go > 0:
        (tmp, go) = re.subn("{[^{}]*}", "", tmp)
    if re.match(".*[{}]", tmp):
        raise ValueError('option "' + arg + '" contains nonmatching braces')
    # @replicas: beginning strings are not treated.
    # the reason is that these strings are expected to be already fixed by possibly adding braces after the :
    if re.match("@replicas:", arg):
        return arg
    if (
        always
        or arg == ""
        or (not comma and re.match(".*[ {}]", arg))
        or (comma and re.match(".*[ {},]", arg))
    ):
        return "{" + arg + "}"
    return arg


def _format_single(builder, arg, level=0):
    """Format a single argument.

    Level keeps track of recursions.
    """
    import numbers

    if builder._pre_format is not None:
        arg = builder._pre_format(arg)
    # only import if already loaded
    # this is to avoid slow startup times.
    if builder._enable_mda_groups and "MDAnalysis" in sys.modules:
        import MDAnalysis

        if isinstance(arg, MDAnalysis.core.groups.AtomGroup):
            arg = arg.indices + 1
    if isinstance(arg, str):
        return re.sub("[\n\t]", " ", arg)
    if isinstance(arg, numbers.Number):
        return str(arg)
    if isinstance(arg, dict):
        raise TypeError("options cannot be a dictionary")
    if isinstance(arg, _replicas):
        if level > 1:
            raise TypeError("@replica syntax only allowed for scalar or rank 1 vectors")
        return "@replicas:" + _fix_braces(
            builder, _format_single(builder, arg.arg, level + 1), comma=True
        )
    if hasattr(arg, "__iter__"):
        string = ""
        for x in arg:
            string += builder._separator + _fix_braces(
                builder, _format_single(builder, x, level + 1), comma=True
            )
        return string[len(builder._separator) :]
    raise TypeError("options should be string, number or iterable")


def _format_numbered(builder, name, arg):
    """Format a numbered argument."""
    if isinstance(arg, _numbered):
        arg = arg.arg
        ret = ""
        for x in arg:
            if not isinstance(x, int):
                raise TypeError("numbered types should have integer keys")
            ret += _format_opt(builder, name + str(x), arg[x])
        return ret
    return _format_opt(builder, name, arg)


def _format_flag(builder, name, arg):
    """Format a flag."""
    if arg is None:
        return ""
    if isinstance(arg, bool):
        if arg:
            return " " + name
        else:
            return ""
    raise TypeError(name + " should be of bool type")


def _format_label(builder, label):
    """Format a label."""
    if label is None:
        return ""
    if isinstance(label, str):
        if len(label) == 0:
            return ""
        return label + ": "
    raise TypeError("label should be of str type")


def _format_opt(builder, name, arg):
    """Format an option."""
    if arg is None:
        return ""
    if isinstance(arg, _replicas):
        string = "@replicas:" + _fix_braces(
            builder, _format_single(builder, arg.arg), comma=False
        )
    else:
        string = _fix_braces(builder, _format_single(builder, arg), comma=False)
    return " " + name + "=" + string


def _format_return(builder, ret, action=None):
    """Format the return statement"""
    if builder._toplumed is not None:
        try:
            builder._toplumed.cmd("readInputLine", ret)
        except:
            builder.last_failure = ret
            raise
    if builder._tofile is not None:
        builder._tofile.write(ret + "\n")
    builder.history.append(ret)
    # LOAD changes the attached instance's registry only after the command succeeds.
    if builder._toplumed is not None and builder._has_dict:
        # Verbatim input also reaches this function without a known action name.
        if action == "LOAD" or (
            action is None and re.match(r"^\s*(?:[^\s:]+:\s*)?LOAD(?:\s|$)", ret)
        ):
            builder.refresh()
    if builder._post_format is not None:
        return builder._post_format(ret)
    return ret + "\n"


def _format_verbatim(builder, verbatim):
    if verbatim is not None:
        if len(verbatim) > 0:
            return " " + re.sub("[\n\t]", " ", verbatim)
    return ""


def _format_anything(builder, name, arg):
    """Choose format based on arg type"""
    ret = ""
    if name == "verbatim":
        ret += _format_verbatim(builder, arg)
    elif isinstance(arg, bool):
        ret += _format_flag(builder, name, arg)
    elif isinstance(arg, _numbered):
        ret += _format_numbered(builder, name, arg)
    else:
        ret += _format_opt(builder, name, arg)
    return ret


# this functions are made available for the exec commands
_format_tools = {
    "_format_numbered": _format_numbered,
    "_format_flag": _format_flag,
    "_format_label": _format_label,
    "_format_opt": _format_opt,
    "_format_return": _format_return,
    "_format_verbatim": _format_verbatim,
    "_format_anything": _format_anything,
}


class InputBuilder:
    """Object used to construct plumed input files.

    An instance of this object can be used to construct plumed input files.
    Check the constructor to see all the available options.
    """

    def __init__(
        self,
        *,
        tofile=None,
        toplumed=None,
        kernel=None,
        append_underscores=False,
        comma_separator=False,
        minimize_braces=True,
        enable_at=True,
        load_dict=True,
        enable_mda_groups=True,
        post_format=None,
        pre_format=None,
    ):
        """Constructor.

        Parameters
        ----------
        tofile:
          PLUMED input is also forwarded to this file. Useful for preparing input files.
        toplumed:
          PLUMED input is also forwarded to this PLUMED object. This is useful for
          interactively testing an input. Its live metadata is used when supported,
          including previously loaded actions. A successful LOAD through this builder
          refreshes its methods. The supplied object is borrowed and remains usable.
        kernel:
          Path to PLUMED kernel used when no toplumed object is supplied, or for legacy
          Vim-file lookup. By default, the default kernel (PLUMED_KERNEL) is used.
        append_underscores:
          Append the two underscores to method names to avoid autocompletion problems in older ipython versions.
        comma_separator:
          Use comma as a separator rather than a space.
        minimize_braces:
          Minimize the number of braces. Setting this to False triggers a bug in the PLUMED parser
          (fixed in 2.4.5 and 2.5.1).
        enable_at:
          Enable `ib.at.chi(1,"A")` syntax.
        load_dict:
          Load full dictionary. Set to false to make initialization faster, at the price of loosing autocompletion.
        enable_mda_groups:
          Enable MDAnalysis groups. Notice that MDAnalysis is not explicitly imported. Only set to false if you
          have a non-working MDAnalysis module installed.
        """
        self._toplumed = toplumed
        self._tofile = tofile
        if comma_separator:
            self._separator = ","
        else:
            self._separator = " "
        self._minimize_braces = minimize_braces
        self._append_underscores = append_underscores
        # history of plumed lines
        self.history = []
        # last line leading to a failure in a Plumed object
        self.last_failure = ""
        self._has_dict = False
        self._enable_mda_groups = enable_mda_groups
        self._pre_format = pre_format
        self._post_format = post_format
        if load_dict:
            metadata = _readkernelmetadata(kernel=kernel, instance=toplumed)
            if metadata is None:
                metadata = _readvimdict(kernel=kernel)
            self._set_dictionary(*metadata)

        if enable_at:
            self.at = _at(self)

    def _set_dictionary(self, dictionary, doc):
        """Replace generated methods, syntax and help with one metadata snapshot."""
        functions = _create_functions(
            dictionary, doc=doc, append_underscores=self._append_underscores
        )
        methods = {}
        for name, source in functions.items():
            namespace = {}
            exec(source, _format_tools, namespace)
            methods[name] = types.MethodType(namespace[name], self)
        # Compile everything first so a malformed snapshot leaves existing methods intact.
        for name in self.__dict__.get("_functions", {}):
            if name not in functions:
                self.__dict__.pop(name, None)
        self.__dict__.update(methods)
        self._vimdict, self._doc = dictionary, doc
        self._functions = functions
        self._has_dict = True

    def refresh(self):
        """Refresh syntax, help and methods from the attached PLUMED instance.

        Use after loading an action directly through that instance. LOAD sent
        through this builder refreshes automatically. Return True on refresh,
        or False for an unattached builder, disabled dictionary loading, or
        an older kernel using static Vim metadata. The instance is not finalized.
        """
        if self._toplumed is None or not self._has_dict:
            return False
        metadata = _readkernelmetadata(instance=self._toplumed)
        if metadata is None:
            return False
        self._set_dictionary(*metadata)
        return True

    def __call__(self, action, LABEL="", verbatim=None, **kwargs):
        ret = ""
        ret += _format_label(self, LABEL)
        ret += action
        retlist = []
        for arg in sorted(kwargs):
            retlist.append(_format_anything(self, arg, kwargs[arg]))
        retlist.sort()
        for x in retlist:
            if len(x) > 0:
                ret += x
        ret += _format_verbatim(self, verbatim)
        return _format_return(self, ret, action)

    def __getattr__(self, name):
        if self._has_dict:

            class _callme:
                def __init__(self, builder, name):
                    self._builder = builder
                    self._name = name
                    self.__doc__ = getattr(self._builder, name).__doc__

                def __call__(self, LABEL="", verbatim=None, **kwargs):
                    func = getattr(self._builder, self._name)
                    return func(LABEL, verbatim, **kwargs)

            if self._append_underscores and not re.match(".*__$", name):
                name__ = name + "__"
                if name__ in self._functions.keys():
                    return _callme(self, name__)
            if not self._append_underscores and re.match(".*__$", name):
                name__ = re.sub("__$", "", name)
                if name__ in self._functions.keys():
                    return _callme(self, name__)
            raise AttributeError("unknown method " + name)
        else:

            class _callme:
                def __init__(self, builder, name):
                    self._builder = builder
                    self._name = name
                    self.__doc__ = "dynamic method printing " + name

                def __call__(self, LABEL="", verbatim=None, **kwargs):
                    return self._builder(self._name, LABEL, verbatim, **kwargs)

            if re.match(".*[^_]__$", name):
                name = name[:-2]
            return _callme(self, name)

    def verbatim(self, line):
        """Create an arbitrary line."""
        return _format_return(self, re.sub("[\n\t]", " ", line))

    def numbered(self, arg):
        """Shortcut for numbered syntax.

        Accepts either a list/tuple or a dictionary with integer keys.
        """
        return _numbered(arg)

    def replicas(self, arg):
        """Shortcut for replica syntax.

        Accepts a list/tuple.
        """
        return _replicas(arg)
