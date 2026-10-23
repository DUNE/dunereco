#!/usr/bin/env python3
'''
wcls-summary.py - summarize a Wire-Cell Toolkit / LArSoft (WCLS) job
configuration as a human-oriented PDF document.

The script:

1. Parses a FHiCL file (enough of the language to follow #include,
   PROLOG, @local::, @table::, @sequence:: and dotted overrides) and
   finds a table with a "wcls_main" WCLS tool configuration.

2. Converts wcls_main.params to Jsonnet std.extVar() strings
   (--ext-str) and wcls_main.structs to std.extVar() code (--ext-code)
   in the same way that larwirecell's WCLS tool does.

3. Compiles the named Wire-Cell Jsonnet configuration with the
   "jsonnet" command line program.

4. Inspects the resulting list of WCT component configurations,
   determines which job types are present (simulation, noise
   filtering, signal processing, ...) and writes a LaTeX document with
   tables and plots: detector geometry, physics constants, field,
   electronics and total responses, noise spectra, deconvolution
   filters, etc.  Each table names the parameters, their values and
   the WCT component "type:name" instances that provide them.

5. Runs pdflatex to produce the PDF.

Example, running from a dunereco source checkout (not installed):

  ./wcls-summary.py \\
      --dunewirecell . \\
      --fcl-tree /path/to/dunecore \\
      -J /path/to/wire-cell-toolkit/cfg \\
      -P /path/to/wire-cell-data \\
      -o summary.pdf \\
      wirecell_dune.fcl  dune10kt-vd/wcls-sim-drift-simchannel-nf-sp.jsonnet

The FHiCL table can be selected with --block (dotted path, eg
"physics.producers.tpcrawdecoder" or "dune10kt_vertdrift_sim_nfsp").
If not given, the first table whose wcls_main.configs matches the
given Jsonnet file is used (or the only WCLS table if no Jsonnet file
is given).  Use --list to see all WCLS tables found.

Values that the FHiCL parser can not resolve (eg @local:: references
to files not found) may be supplied with --ext-str / --ext-code which
always override values from FHiCL.

Requirements: python3 with numpy and matplotlib, a "jsonnet" program
(go-jsonnet or C++ jsonnet), pdflatex and optionally graphviz "dot".
'''

import argparse
import bz2
import collections
import copy
import datetime
import glob
import gzip
import json
import math
import os
import re
import shutil
import subprocess
import sys
import tempfile

import numpy as np

os.environ.setdefault('MPLCONFIGDIR', tempfile.mkdtemp(prefix='wcls-summary-mpl-'))
import matplotlib                       # noqa: E402
matplotlib.use('Agg')
import matplotlib.pyplot as plt         # noqa: E402
from matplotlib import colors as mcolors  # noqa: E402
from matplotlib.patches import Rectangle  # noqa: E402


# ----------------------------------------------------------------------
# WCT system of units (see WireCellUtil/Units.h)
# ----------------------------------------------------------------------

class U:
    mm = 1.0
    um = 1e-3 * mm
    cm = 10.0 * mm
    m = 1000.0 * mm
    ns = 1.0
    us = 1e3 * ns
    ms = 1e6 * ns
    s = 1e9 * ns
    eplus = 1.0
    coulomb = eplus / 1.602176634e-19
    fC = 1e-15 * coulomb
    MeV = 1.0
    joule = MeV / 1.602176634e-13
    kg = joule * s ** 2 / m ** 2
    gram = 1e-3 * kg
    Bq = 1.0 / s
    volt = 1e-6
    mV = 1e-3 * volt
    Hz = 1.0 / s
    kHz = 1e3 * Hz
    MHz = 1e6 * Hz
    kelvin = 1.0


# ----------------------------------------------------------------------
# Logging
# ----------------------------------------------------------------------

WARNINGS = []


def warn(msg):
    if msg in WARNINGS:
        return
    WARNINGS.append(msg)
    sys.stderr.write('wcls-summary: warning: %s\n' % msg)


def info(msg):
    sys.stderr.write('wcls-summary: %s\n' % msg)


# ----------------------------------------------------------------------
# FHiCL parsing
# ----------------------------------------------------------------------

class FhiclError(Exception):
    pass


class Atom:
    '''A FHiCL atomic value.

    kind is one of: str (quoted), word (unquoted), num, bool, nil,
    complex.  src is "file:line" of the literal.  via lists the
    @local:: (etc) references through which the value was obtained.'''
    __slots__ = ('value', 'text', 'kind', 'src', 'via')

    def __init__(self, value, text, kind, src, via=()):
        self.value = value
        self.text = text
        self.kind = kind
        self.src = src
        self.via = tuple(via)

    def __repr__(self):
        return 'Atom(%r)' % (self.text,)


class Unresolved:
    '''A reference that could not be resolved.'''
    __slots__ = ('ref', 'src', 'via')

    def __init__(self, ref, src, via=()):
        self.ref = ref
        self.src = src
        self.via = tuple(via)

    def __repr__(self):
        return 'Unresolved(%s)' % self.ref


class Table(collections.OrderedDict):
    pass


ERASE = object()

_TOKRE = re.compile(r'''
 (?P<ws>\s+)
|(?P<include>\#include\s+"(?P<incfile>[^"]*)")
|(?P<comment>(?:\#|//)[^\n]*)
|(?P<dstr>"(?:[^"\\]|\\.)*")
|(?P<sstr>'[^']*')
|(?P<protect>@protect_ignore:|@protect_error:)
|(?P<ref>@(?:local|table|sequence|id|db)::[A-Za-z0-9_.\[\]]+)
|(?P<erase>@erase)
|(?P<nil>@nil)
|(?P<num>[-+]?(?:\d+\.\d*|\.\d+|\d+)(?:[eE][-+]?\d+)?(?![A-Za-z_]))
|(?P<name>[-+]?[A-Za-z_][A-Za-z0-9_]*(?:\.[A-Za-z_][A-Za-z0-9_]*|\[\d+\])*)
|(?P<punct>[{}\[\]:,()])
''', re.X)

_PATHSEG = re.compile(r'([A-Za-z_][A-Za-z0-9_]*)|\[(\d+)\]')


def split_path(path):
    segs = []
    pos = 0
    path = path.strip()
    while pos < len(path):
        if path[pos] == '.':
            pos += 1
            continue
        m = _PATHSEG.match(path, pos)
        if not m:
            raise FhiclError('bad FHiCL name path: %r' % path)
        segs.append(m.group(1) if m.group(1) is not None else int(m.group(2)))
        pos = m.end()
    return segs


def copy_val(v, ref):
    '''Deep copy a value obtained via reference ref.'''
    if isinstance(v, Atom):
        return Atom(v.value, v.text, v.kind, v.src, (ref,) + v.via)
    if isinstance(v, Unresolved):
        return Unresolved(v.ref, v.src, (ref,) + v.via)
    if isinstance(v, Table):
        t = Table()
        for k, vv in v.items():
            t[k] = copy_val(vv, ref)
        return t
    if isinstance(v, list):
        return [copy_val(vv, ref) for vv in v]
    return copy.deepcopy(v)


def fhicl_lookup(root, path):
    '''Return value at dotted path or raise KeyError.'''
    cur = root
    for seg in split_path(path):
        if isinstance(seg, int):
            if not isinstance(cur, list) or seg >= len(cur):
                raise KeyError(path)
            cur = cur[seg]
        else:
            if not isinstance(cur, Table) or seg not in cur:
                raise KeyError(path)
            cur = cur[seg]
    return cur


def fhicl_assign(root, path, value):
    segs = split_path(path)
    cur = root
    for i, seg in enumerate(segs):
        last = (i == len(segs) - 1)
        if isinstance(seg, int):
            if not isinstance(cur, list):
                raise FhiclError('index into non-sequence in %r' % path)
            while len(cur) <= seg:
                cur.append(Atom(None, 'nil', 'nil', '?'))
            if last:
                cur[seg] = value
            else:
                if not isinstance(cur[seg], (Table, list)):
                    cur[seg] = Table()
                cur = cur[seg]
        else:
            if not isinstance(cur, Table):
                raise FhiclError('member of non-table in %r' % path)
            if last:
                if value is ERASE:
                    cur.pop(seg, None)
                else:
                    if seg in cur:
                        cur.pop(seg)  # re-insert to keep "last definition" order
                    cur[seg] = value
            else:
                nxt = cur.get(seg)
                if not isinstance(nxt, (Table, list)):
                    nxt = Table()
                    cur[seg] = nxt
                cur = nxt


class FhiclDocument:
    '''Parse a FHiCL file and its includes into one namespace.'''

    def __init__(self, include_dirs=(), trees=()):
        self.include_dirs = list(include_dirs)
        env = os.environ.get('FHICL_FILE_PATH', '')
        self.include_dirs += [d for d in env.split(':') if d]
        self.tree_index = {}
        for tree in trees:
            for dp, dns, fns in os.walk(tree):
                dns[:] = [d for d in dns if not d.startswith('.')]
                for fn in fns:
                    if fn.endswith('.fcl'):
                        self.tree_index.setdefault(fn, os.path.join(dp, fn))
        self.root = Table()
        self.prolog_keys = set()
        self.files = []         # all files loaded
        self.missing = []       # (include, from)
        self._stack = []

    def find_include(self, name, fromdir):
        if os.path.isabs(name):
            return name if os.path.exists(name) else None
        for d in [fromdir] + self.include_dirs:
            p = os.path.join(d, name)
            if os.path.exists(p):
                return p
        base = os.path.basename(name)
        return self.tree_index.get(base)

    def load(self, path):
        path = os.path.abspath(path)
        if path in self._stack:
            warn('recursive include of %s ignored' % path)
            return
        self._stack.append(path)
        self.files.append(path)
        try:
            text = open(path, encoding='utf-8', errors='replace').read()
            _FhiclParser(self, text, path).parse_document()
        except FhiclError as err:
            warn('FHiCL parse error, rest of file ignored: %s' % err)
        finally:
            self._stack.pop()

    def resolve(self, ref, src):
        kind, path = ref[1:].split('::', 1)
        try:
            val = fhicl_lookup(self.root, path)
        except KeyError:
            return Unresolved(ref, src)
        return copy_val(val, ref)


class _FhiclParser:

    def __init__(self, doc, text, path):
        self.doc = doc
        self.path = path
        self.dir = os.path.dirname(path)
        self.toks = []
        self.pos = 0
        self.in_prolog = False
        linestarts = [0] + [m.end() for m in re.finditer('\n', text)]
        import bisect
        pos = 0
        while pos < len(text):
            m = _TOKRE.match(text, pos)
            if not m:
                line = bisect.bisect_right(linestarts, pos)
                raise FhiclError('%s:%d: can not tokenize near %r'
                                 % (path, line, text[pos:pos + 30]))
            kind = m.lastgroup
            if kind == 'incfile':
                kind = 'include'
            if kind not in ('ws', 'comment'):
                line = bisect.bisect_right(linestarts, pos)
                val = m.group('incfile') if kind == 'include' else m.group(kind)
                self.toks.append((kind, val, line))
            pos = m.end()

    # token helpers
    def peek(self, off=0):
        i = self.pos + off
        return self.toks[i] if i < len(self.toks) else ('eof', None, -1)

    def next(self):
        t = self.peek()
        self.pos += 1
        return t

    def expect(self, kind, val=None):
        t = self.next()
        if t[0] != kind or (val is not None and t[1] != val):
            raise FhiclError('%s:%d: expected %s %s got %s %r'
                             % (self.path, t[2], kind, val or '', t[0], t[1]))
        return t

    def src(self, tok):
        return '%s:%d' % (os.path.basename(self.path), tok[2])

    def parse_document(self):
        root = self.doc.root
        while True:
            tok = self.peek()
            kind, val, line = tok
            if kind == 'eof':
                return
            if kind == 'include':
                self.next()
                inc = self.doc.find_include(val, self.dir)
                if inc is None:
                    self.doc.missing.append((val, self.src(tok)))
                    warn('FHiCL include "%s" not found (from %s)' % (val, self.src(tok)))
                else:
                    self.doc.load(inc)
                continue
            if kind == 'name' and val in ('BEGIN_PROLOG', 'END_PROLOG'):
                self.next()
                self.in_prolog = (val == 'BEGIN_PROLOG')
                continue
            if kind == 'ref' and val.startswith('@table::'):
                self.next()
                self.splice_table(root, self.doc.resolve(val, self.src(tok)), tok)
                continue
            self.parse_binding(root)

    def splice_table(self, table, val, tok):
        if isinstance(val, Table):
            for k, v in val.items():
                table.pop(k, None)
                table[k] = v
        else:
            warn('%s: can not splice unresolved %s' % (self.src(tok), tok[1]))

    def parse_binding(self, table):
        tok = self.next()
        if tok[0] != 'name':
            raise FhiclError('%s:%d: expected a name, got %r' % (self.path, tok[2], tok[1]))
        name = tok[1]
        if self.peek()[0] == 'protect':
            self.next()
        else:
            self.expect('punct', ':')
            if self.peek()[0] == 'protect':  # eg "a: @protect_ignore: 1" (lenient)
                self.next()
        value = self.parse_value()
        if table is self.doc.root and self.in_prolog:
            self.doc.prolog_keys.add(split_path(name)[0])
        fhicl_assign(table, name, value)

    def parse_value(self):
        tok = self.next()
        kind, val, line = tok
        src = self.src(tok)
        if kind == 'punct' and val == '{':
            table = Table()
            while True:
                t = self.peek()
                if t[0] == 'punct' and t[1] == '}':
                    self.next()
                    return table
                if t[0] == 'eof':
                    raise FhiclError('%s: unterminated table' % src)
                if t[0] == 'ref' and t[1].startswith('@table::'):
                    self.next()
                    self.splice_table(table, self.doc.resolve(t[1], self.src(t)), t)
                    continue
                if t[0] == 'punct' and t[1] == ',':  # lenient
                    self.next()
                    continue
                self.parse_binding(table)
        if kind == 'punct' and val == '[':
            seq = []
            while True:
                t = self.peek()
                if t[0] == 'punct' and t[1] == ']':
                    self.next()
                    return seq
                if t[0] == 'punct' and t[1] == ',':
                    self.next()
                    continue
                if t[0] == 'eof':
                    raise FhiclError('%s: unterminated sequence' % src)
                if t[0] == 'ref' and t[1].startswith('@sequence::'):
                    self.next()
                    v = self.doc.resolve(t[1], self.src(t))
                    if isinstance(v, list):
                        seq.extend(v)
                    else:
                        seq.append(v)
                    continue
                seq.append(self.parse_value())
        if kind == 'punct' and val == '(':
            re_ = self.next()
            self.expect('punct', ',')
            im_ = self.next()
            self.expect('punct', ')')
            text = '(%s,%s)' % (re_[1], im_[1])
            return Atom(text, text, 'complex', src)
        if kind == 'ref':
            if val.startswith('@id::'):
                return Atom(val, val, 'word', src)
            return self.doc.resolve(val, src)
        if kind == 'erase':
            return ERASE
        if kind == 'nil':
            return Atom(None, 'nil', 'nil', src)
        if kind == 'dstr':
            s = val[1:-1]
            try:
                s = json.loads(val)
            except ValueError:
                pass
            return Atom(s, val, 'str', src)
        if kind == 'sstr':
            return Atom(val[1:-1], val, 'str', src)
        if kind == 'num':
            try:
                v = int(val)
            except ValueError:
                v = float(val)
            return Atom(v, val, 'num', src)
        if kind == 'name':
            if val in ('true', 'false'):
                return Atom(val == 'true', val, 'bool', src)
            if val == 'nil':
                return Atom(None, val, 'nil', src)
            if val.lstrip('+-') == 'infinity':
                return Atom(float(val.replace('infinity', 'inf')), val, 'num', src)
            return Atom(val, val, 'word', src)
        raise FhiclError('%s: unexpected token %r' % (src, val))


def fhicl_to_python(v):
    if isinstance(v, Table):
        return {k: fhicl_to_python(vv) for k, vv in v.items()}
    if isinstance(v, list):
        return [fhicl_to_python(vv) for vv in v]
    if isinstance(v, Atom):
        return v.value
    if isinstance(v, Unresolved):
        return '<unresolved %s>' % v.ref
    return v


def find_wcls_blocks(doc):
    '''Return list of (dotted path, table) of tables holding a wcls_main.'''
    found = []

    def walk(tab, path):
        for k, v in tab.items():
            if not isinstance(v, Table):
                continue
            p = path + [k]
            wm = v.get('wcls_main')
            if isinstance(wm, Table) and 'configs' in wm:
                found.append(('.'.join(p), v))
            walk(v, p)
    walk(doc.root, [])
    # tables outside of prolog first
    found.sort(key=lambda pv: pv[0].split('.')[0] in doc.prolog_keys)
    return found


def fhicl_strings(v):
    if isinstance(v, list):
        return [a.value for a in v if isinstance(a, Atom)]
    if isinstance(v, Atom):
        return [v.value]
    return []


class ExtVar:
    def __init__(self, name, kind, value, atom=None, override=False):
        self.name = name
        self.kind = kind        # 'str' or 'code'
        self.value = value
        self.atom = atom
        self.override = override

    @property
    def origin(self):
        if self.override:
            return 'command line'
        a = self.atom
        if a is None:
            return ''
        via = ' via ' + ', '.join(a.via) if a.via else ''
        return a.src + via


def make_extvars(wcls_main, overrides_str, overrides_code):
    '''Turn wcls_main params/structs into ext vars like WCLS_tool does.'''
    evs = collections.OrderedDict()
    unresolved = []
    params = wcls_main.get('params')
    if isinstance(params, Table):
        for k, a in params.items():
            if isinstance(a, Unresolved):
                unresolved.append((k, a))
                evs[k] = ExtVar(k, 'str', None, a)
                continue
            if isinstance(a, Atom):
                val = a.text if a.kind in ('num', 'bool') else ('' if a.value is None else str(a.value))
                evs[k] = ExtVar(k, 'str', val, a)
            else:
                warn('wcls_main.params.%s is not an atom, passing as JSON' % k)
                evs[k] = ExtVar(k, 'str', json.dumps(fhicl_to_python(a)), None)
    structs = wcls_main.get('structs')
    if isinstance(structs, Table):
        for k, a in structs.items():
            if isinstance(a, Unresolved):
                unresolved.append((k, a))
                evs[k] = ExtVar(k, 'code', None, a)
                continue
            if isinstance(a, Atom):
                if a.kind == 'bool':
                    code = 'true' if a.value else 'false'
                elif a.kind == 'nil':
                    code = 'null'
                elif a.kind == 'num':
                    code = a.text.lstrip('+')
                    if code.lower().endswith('infinity'):
                        warn('struct %s is infinite, not representable in Jsonnet' % k)
                else:
                    code = str(a.value)  # WCLS passes the string content as code
                evs[k] = ExtVar(k, 'code', code, a)
            else:
                warn('wcls_main.structs.%s is not an atom, passing as JSON code' % k)
                evs[k] = ExtVar(k, 'code', json.dumps(fhicl_to_python(a)), None)
    for kind, ovs in (('str', overrides_str), ('code', overrides_code)):
        for ov in ovs:
            if '=' not in ov:
                raise SystemExit('bad override, need key=value: %s' % ov)
            k, v = ov.split('=', 1)
            evs[k] = ExtVar(k, kind, v, None, override=True)
            unresolved = [u for u in unresolved if u[0] != k]
    return evs, unresolved


# ----------------------------------------------------------------------
# Jsonnet and data file handling
# ----------------------------------------------------------------------

class Paths:
    '''Search paths for Jsonnet files and WCT data files.'''

    def __init__(self, args, tmpdir):
        self.jpaths = []
        # Synthetic tree of symlinks for --map, highest priority.
        maps = list(args.map or [])
        if args.dunewirecell:
            d = os.path.abspath(args.dunewirecell)
            maps.append('pgrapher/experiment=' + d)
            for sub in ('common', 'dune', 'ui'):
                if os.path.isdir(os.path.join(d, sub)):
                    maps.append('pgrapher/%s=%s' % (sub, os.path.join(d, sub)))
        if maps:
            mroot = os.path.join(tmpdir, 'jmap')
            for m in maps:
                if '=' not in m:
                    raise SystemExit('bad --map, need prefix=dir: %s' % m)
                pre, d = m.split('=', 1)
                link = os.path.join(mroot, pre)
                os.makedirs(os.path.dirname(link), exist_ok=True)
                if os.path.lexists(link):
                    os.unlink(link)
                os.symlink(os.path.abspath(d), link)
            self.jpaths.append(mroot)
        self.jpaths += [os.path.abspath(j) for j in (args.jpath or [])]
        self.wcpaths = [p for p in os.environ.get('WIRECELL_PATH', '').split(':') if p]
        self.datapaths = [os.path.abspath(p) for p in (args.data_path or [])]
        self.data_index = {}
        self.first_index = {}
        if args.dunewirecell:
            self._index(args.dunewirecell, self.first_index)
        for tree in args.data_tree or []:
            self._index(tree, self.data_index)

    @staticmethod
    def _index(tree, index):
        '''Index WCT data files under tree by base name.'''
        for dp, dns, fns in os.walk(tree):
            dns[:] = [x for x in dns if not x.startswith('.') and x != 'pgrapher']
            for fn in fns:
                if re.search(r'\.json(\.bz2|\.gz)?$', fn):
                    index.setdefault(fn, os.path.join(dp, fn))

    def all_jpaths(self, extra=()):
        return self.jpaths + list(extra) + self.wcpaths

    def find_jsonnet(self, name, extra=()):
        if os.path.exists(name):
            return os.path.abspath(name)
        for d in self.all_jpaths(extra):
            p = os.path.join(d, name)
            if os.path.exists(p):
                return p
        return None

    def find_data(self, name):
        if not name:
            return None
        if os.path.isabs(name) and os.path.exists(name):
            return name
        if os.path.basename(name) in self.first_index:
            return self.first_index[os.path.basename(name)]
        for d in self.datapaths + self.jpaths + self.wcpaths:
            p = os.path.join(d, name)
            if os.path.exists(p):
                return p
        return self.data_index.get(os.path.basename(name))


_JSON_CACHE = {}


def load_json(path):
    if path in _JSON_CACHE:
        return _JSON_CACHE[path]
    if path.endswith('.bz2'):
        fp = bz2.open(path, 'rt')
    elif path.endswith('.gz'):
        fp = gzip.open(path, 'rt')
    else:
        fp = open(path)
    with fp:
        dat = json.load(fp)
    _JSON_CACHE[path] = dat
    return dat


def compile_jsonnet(jsonnet_exe, jfile, jpaths, extvars):
    '''Return (list of component configs, command line).'''
    cmd = [jsonnet_exe]
    # jsonnet: the right-most -J has precedence, WCT: first path wins.
    for d in reversed(jpaths):
        cmd += ['-J', d]
    for ev in extvars.values():
        if ev.value is None:
            continue
        cmd += ['--ext-str' if ev.kind == 'str' else '--ext-code',
                '%s=%s' % (ev.name, ev.value)]
    cmd.append(jfile)
    info('compiling %s' % jfile)
    proc = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    if proc.returncode != 0:
        sys.stderr.write(proc.stderr)
        raise SystemExit('jsonnet compilation failed')
    dat = json.loads(proc.stdout)
    if isinstance(dat, dict):
        dat = [dat]
    return dat, cmd


# ----------------------------------------------------------------------
# WCT configuration model
# ----------------------------------------------------------------------

def tn(comp):
    name = comp.get('name') or ''
    return comp['type'] + (':' + name if name else '')


class WCConfig:

    def __init__(self, comps):
        self.comps = [c for c in comps if isinstance(c, dict) and 'type' in c]
        self.by_tn = {}
        for c in self.comps:
            self.by_tn[tn(c)] = c

    def get(self, ref):
        if not isinstance(ref, str) or not ref:
            return None
        c = self.by_tn.get(ref)
        if c is None and ref.endswith(':'):
            c = self.by_tn.get(ref[:-1])
        if c is None and ':' not in ref:
            c = self.by_tn.get(ref + ':')
        return c

    def of_type(self, *types):
        return [c for c in self.comps if c['type'] in types]

    def types(self):
        return collections.Counter(c['type'] for c in self.comps)

    def apps(self):
        return [c for c in self.comps if isinstance(c.get('data'), dict)
                and 'edges' in c['data']]


JOB_TYPES = [
    ('sim', 'Simulation',
     {'Drifter', 'DepoSetDrifter', 'DepoTransform', 'Ductor', 'MultiDuctor', 'DepoSplat',
      'DepoZipper', 'AddNoise', 'Digitizer', 'DepoBagger', 'TrackDepos', 'BlipSource',
      'wclsSimDepoSource', 'wclsSimDepoSetSource', 'wclsDepoFluxWriter', 'wclsSimChannelSink'}),
    ('nf', 'Noise filtering', {'OmnibusNoiseFilter', 'OmniChannelNoiseDB'}),
    ('sp', 'Signal processing', {'OmnibusSigProc', 'L1SPFilter', 'L1SPFilterPD'}),
    ('dnn', 'DNN ROI finding', {'DNNROIFinding', 'TorchService', 'TorchScript', 'TritonService'}),
    ('img', 'Imaging', {'MaskSlices', 'SumSlices', 'GridTiling', 'BlobSolving', 'BlobClustering',
                        'BlobGrouping', 'ChargeSolving', 'LocalGeomClustering'}),
    ('clus', 'Clustering', {'MultiAlgBlobClustering', 'PointTreeBuilding', 'UbooneClusterSource'}),
]


def detect_jobs(cfg):
    have = set(cfg.types())
    return [(key, title, sorted(have & types)) for key, title, types in JOB_TYPES if have & types]


# ----------------------------------------------------------------------
# Formatting helpers
# ----------------------------------------------------------------------

def tex(s):
    '''Escape text for LaTeX and allow line breaks at punctuation.'''
    s = str(s)
    rep = {'\\': r'\textbackslash{}', '&': r'\&', '%': r'\%', '$': r'\$', '#': r'\#',
           '_': r'\_', '{': r'\{', '}': r'\}', '~': r'\textasciitilde{}',
           '^': r'\textasciicircum{}', '<': r'\textless{}', '>': r'\textgreater{}'}
    out = []
    for ch in s:
        out.append(rep.get(ch, ch))
        if ch in '/_:.,-=' or ch == '[':
            out.append(r'\allowbreak{}')
    return ''.join(out)


def cell(s):
    '''Escape a table cell unless it already is LaTeX (contains a backslash).'''
    s = str(s)
    if '\\' in s:
        return s
    return tex(s).replace(r'\allowbreak{}', '')


def ttex(s):
    return r'\texttt{%s}' % tex(s)


def num(v, prec=6):
    if isinstance(v, bool):
        return 'true' if v else 'false'
    if isinstance(v, int):
        return str(v)
    if isinstance(v, float):
        if v == 0:
            return '0'
        return '%.*g' % (prec, v)
    return str(v)


def _time_unit(v):
    a = abs(v)
    if a == 0:
        return 'us'
    if a >= 0.1 * U.ms:
        return 'ms'
    if a >= 0.1 * U.us:
        return 'us'
    return 'ns'


_TUNITS = {'ns': (U.ns, 'ns'), 'us': (U.us, r'$\mu$s'), 'ms': (U.ms, 'ms')}

_TIME_KEYS = re.compile(
    r'^(tick|period|shaping|start_time|readout_time|time_offset|toffset|ftoffset|ctoffset|'
    r'start|time_smear|window_start|window_duration|reference_time|g4_ref_time|gate|'
    r'overall_short_padding|long_padding|width|drift_dt|tstart)$')


def phys(key, v, ctype=None, cdata=None):
    '''Return a LaTeX string with a physical-units interpretation of
    value v of parameter key (last path element) or None.'''
    if isinstance(v, bool) or not isinstance(v, (int, float)):
        return None
    k = re.sub(r'\[\d+\]$', '', key.split('.')[-1])
    if k in ('DL', 'DT'):
        return r'%s cm$^2$/s' % num(v / (U.cm ** 2 / U.s), 4)
    if k == 'density':
        return r'%s g/cm$^3$' % num(v / (U.gram / U.cm ** 3), 4)
    if k == 'ar39activity':
        return '%s Bq/kg' % num(v / (U.Bq / U.kg), 4)
    if k == 'lifetime':
        return '%s ms' % num(v / U.ms, 4)
    if k in ('drift_speed', 'speed'):
        return r'%s mm/$\mu$s' % num(v / (U.mm / U.us), 5)
    if k in ('anode', 'response', 'cathode') and ctype == 'AnodePlane':
        return '%s cm' % num(v / U.cm, 5)
    if k == 'xregions' or (k in ('anode', 'response', 'cathode') and ctype in ('Drifter',)):
        return '%s cm' % num(v / U.cm, 5)
    if k == 'gain' and v != 0 and abs(v) < 1e-6:
        return '%s mV/fC' % num(v / (U.mV / U.fC), 4)
    if k in ('baselines', 'fullscale'):
        return '%s mV' % num(v / U.mV, 5)
    if k == 'ADC_mV':
        return '%s ADC/mV' % num(v * U.mV, 5)
    if ctype in ('HfFilter', 'LfFilter') and k in ('sigma', 'tau', 'max_freq'):
        mf = (cdata or {}).get('max_freq', 1)
        if mf != 1:
            return '%s MHz' % num(v / U.MHz, 4)
        return None
    if k in ('pitch', 'location', 'origin'):
        return '%s mm' % num(v / U.mm, 5)
    if k == 'wire_length_scale':
        return '%s cm' % num(v / U.cm, 4)
    if _TIME_KEYS.match(k):
        if k == 'start' and ctype not in (None, 'ColdElecResponse', 'JsonElecResponse',
                                          'ResponseSys', 'RCResponse', 'WarmElecResponse'):
            return None
        if k == 'width' and ctype not in ('RCResponse',):
            return None
        fac, lab = _TUNITS[_time_unit(v)]
        return '%s %s' % (num(v / fac, 6), lab)
    return None


def fmt_value(key, v, ctype=None, cdata=None, maxlen=160):
    '''Format a configuration value as LaTeX, with physical units if known.'''
    if isinstance(v, (list, tuple)) and v and all(isinstance(x, (int, float)) and not isinstance(x, bool) for x in v) and len(v) <= 12:
        ph = [phys(key, x, ctype, cdata) for x in v]
        raw = '[' + ', '.join(num(x) for x in v) + ']'
        s = ttex(raw)
        if all(p is not None for p in ph):
            s += r' $\rightarrow$ [' + ', '.join(ph) + ']'
        return s
    if isinstance(v, (dict, list)):
        txt = json.dumps(v, separators=(',', ':'))
        if len(txt) > maxlen:
            n = len(v)
            txt = txt[:maxlen] + '... (%d %s)' % (n, 'keys' if isinstance(v, dict) else 'elements')
        return ttex(txt)
    if isinstance(v, str):
        if len(v) > maxlen:
            v = v[:maxlen] + '...'
        return ttex(repr(v) if v == '' else v)
    s = ttex(num(v))
    p = phys(key, v, ctype, cdata)
    if p:
        s += r' = ' + p
    return s


def _ranges(nums):
    nums = sorted(set(nums))
    rng = []
    start = prev = nums[0]
    for x in nums[1:] + [None]:
        if x is not None and x == prev + 1:
            prev = x
            continue
        rng.append(str(start) if start == prev else '%d-%d' % (start, prev))
        if x is not None:
            start = prev = x
    return ','.join(rng)


def compress_names(names):
    '''Compress a list of names with integers into ranges.

    "a1", "a2", "a3" -> "a{1-3}".  Names with several numbers compress to
    a cartesian product when possible, eg "crm{0-159}/face{0-1}", or, if
    all numbers in each name are equal, "crm#sigproc# #={0-319}".'''
    templs = collections.OrderedDict()
    other = []
    for n in names:
        parts = re.split(r'(\d+)', n)
        nums = tuple(int(x) for x in parts[1::2])
        if not nums:
            other.append(n)
            continue
        templs.setdefault(tuple(parts[0::2]), []).append(nums)
    out = []

    def fill(texts, fields):
        return ''.join(t + (fields[i] if i < len(fields) else '') for i, t in enumerate(texts))

    for texts, tuples in templs.items():
        tuples = sorted(set(tuples))
        k = len(tuples[0])
        if len(tuples) == 1:
            out.append(fill(texts, [str(x) for x in tuples[0]]))
            continue
        sets = [sorted(set(t[i] for t in tuples)) for i in range(k)]
        if k > 1 and all(len(set(t)) == 1 for t in tuples):
            out.append('%s #={%s}' % (fill(texts, ['#'] * k), _ranges([t[0] for t in tuples])))
            continue
        if np.prod([len(x) for x in sets]) == len(tuples):
            fields = [str(x[0]) if len(x) == 1 else '{%s}' % _ranges(x) for x in sets]
            out.append(fill(texts, fields))
            continue
        # fall back: vary the last number
        byhead = collections.OrderedDict()
        for t in tuples:
            byhead.setdefault(t[:-1], []).append(t[-1])
        for head, lasts in byhead.items():
            out.append(fill(texts, [str(x) for x in head] + ['{%s}' % _ranges(lasts) if len(lasts) > 1 else str(lasts[0])]))
    return out + other


def names_tex(names, maxn=6):
    c = compress_names(names)
    s = ', '.join(c[:maxn])
    if len(c) > maxn:
        s += ', ... (%d more)' % (len(c) - maxn)
    return ttex(s)


def flatten(d, prefix='', depth=0):
    out = collections.OrderedDict()
    if isinstance(d, dict) and depth < 3:
        for k, v in d.items():
            key = '%s.%s' % (prefix, k) if prefix else k
            if isinstance(v, dict) and v and depth < 2 and len(v) <= 30:
                out.update(flatten(v, key, depth + 1))
            elif isinstance(v, list) and v and len(v) <= 4 and all(isinstance(x, dict) for x in v):
                for i, x in enumerate(v):
                    out.update(flatten(x, '%s[%d]' % (key, i), depth + 1))
            else:
                out[key] = v
    else:
        out[prefix or 'value'] = d
    return out


def hashable(v):
    return json.dumps(v, sort_keys=True)


def group_params(comps, keys=None, skip=()):
    '''Return list of (key, [(value, [tn,...]),...]) over components.'''
    flats = [(tn(c), flatten(c.get('data', {}) or {})) for c in comps]
    allkeys = []
    for _, f in flats:
        for k in f:
            if k not in allkeys:
                allkeys.append(k)
    if keys is not None:
        allkeys = [k for k in keys if k in allkeys] + [k for k in allkeys if k not in keys and False]
    rows = []
    for k in allkeys:
        if any(re.match(s, k) for s in skip):
            continue
        vals = collections.OrderedDict()
        for name, f in flats:
            if k not in f:
                continue
            vals.setdefault(hashable(f[k]), (f[k], []))[1].append(name)
        rows.append((k, list(vals.values())))
    return rows


# ----------------------------------------------------------------------
# LaTeX document builder
# ----------------------------------------------------------------------

class Doc:

    def __init__(self, figdir):
        self.parts = []
        self.figdir = figdir
        self.nfig = 0

    def add(self, s):
        self.parts.append(s)

    def section(self, title, label=None):
        self.add(r'\clearpage' + '\n' + r'\section{%s}' % title)

    def subsection(self, title):
        self.add(r'\subsection{%s}' % title)

    def subsubsection(self, title):
        self.add(r'\subsubsection{%s}' % title)

    def para(self, s):
        self.add(s + '\n')

    def table(self, header, rows, widths, caption=None, fontsize=r'\footnotesize'):
        cols = ''.join('p{%s}' % w for w in widths)
        out = [r'{%s' % fontsize, r'\begin{longtable}{%s}' % cols]
        if caption:
            out.append(r'\caption{%s}\\' % caption)
        out.append(r'\toprule')
        out.append(' & '.join(r'\textbf{%s}' % cell(h) for h in header) + r'\\')
        out.append(r'\midrule\endhead')
        out.append(r'\bottomrule\endfoot')
        for r in rows:
            out.append(' & '.join(cell(c) for c in r) + r'\\')
        out.append(r'\end{longtable}}')
        self.add('\n'.join(out))

    def kv_table(self, rows, caption=None):
        '''rows: (param, value-tex, source-tex)'''
        self.table(['Parameter', 'Value', 'Source'], rows,
                   ['0.24\\linewidth', '0.46\\linewidth', '0.24\\linewidth'], caption)

    def group_table(self, comps, caption=None, skip=(), maxvals=6):
        '''Table of parameters of a set of components of one type.'''
        if not comps:
            return
        ctype = comps[0]['type']
        cdata0 = comps[0].get('data', {})
        rows = []
        allnames = [tn(c) for c in comps]
        for key, vals in group_params(comps, skip=skip):
            if len(vals) == 1:
                val, names = vals[0]
                src = 'all %d' % len(names) if len(names) == len(allnames) and len(names) > 1 \
                    else names_tex(names)
                rows.append((ttex(key), fmt_value(key, val, ctype, cdata0), src))
                continue
            def as_str(v):
                if isinstance(v, list) and len(v) == 1 and isinstance(v[0], str):
                    return v[0]
                return v
            if len(vals) > maxvals and all(isinstance(as_str(v), str) for v, _ in vals):
                allv = [as_str(v) for v, _ in vals]
                txt = ', '.join(compress_names(allv))
                if len(txt) > 300:
                    txt = txt[:300] + '...'
                rows.append((ttex(key), ttex(txt) + r' \emph{(%d distinct values)}' % len(vals),
                             'one per instance' if len(vals) == len(allnames) else names_tex(allnames, 3)))
                continue
            for i, (val, names) in enumerate(vals[:maxvals]):
                rows.append((ttex(key) if i == 0 else '', fmt_value(key, val, ctype, cdata0),
                             names_tex(names, 4)))
            if len(vals) > maxvals:
                rows.append(('', r'\emph{... %d more distinct values}' % (len(vals) - maxvals), ''))
        if caption is None:
            caption = r'Parameters of %d %s component(s): %s' % (
                len(comps), ttex(ctype), names_tex(allnames, 4))
        if not rows:
            self.para(r'%s: \emph{no configuration parameters.}' % names_tex(allnames, 4))
            return
        self.kv_table(rows, caption)

    def figure(self, fig, caption, width=r'\linewidth'):
        self.nfig += 1
        fname = 'fig%03d.pdf' % self.nfig
        fig.savefig(os.path.join(self.figdir, fname), bbox_inches='tight')
        plt.close(fig)
        self.add(r'\begin{figure}[H]\centering\includegraphics[width=%s,height=0.8\textheight,keepaspectratio]{%s}'
                 r'\caption{%s}\end{figure}' % (width, fname, caption))

    def include_pdf(self, fname, caption):
        self.add(r'\begin{figure}[H]\centering\includegraphics[width=\linewidth,height=0.85\textheight,keepaspectratio]{%s}'
                 r'\caption{%s}\end{figure}' % (fname, caption))

    def verbatim(self, text):
        self.add(r'\begin{lstlisting}' + '\n' + text + '\n' + r'\end{lstlisting}')

    def render(self, title, subtitle):
        head = r'''\documentclass[10pt]{article}
\usepackage[a4paper,margin=1.8cm]{geometry}
\usepackage[T1]{fontenc}
\usepackage{lmodern}
\usepackage{graphicx}
\usepackage{longtable}
\usepackage{booktabs}
\usepackage{float}
\usepackage{xcolor}
\usepackage{listings}
\usepackage[hidelinks]{hyperref}
\lstset{basicstyle=\ttfamily\scriptsize,breaklines=true,columns=fullflexible,
        breakatwhitespace=false,postbreak=\mbox{\textcolor{gray}{$\hookrightarrow$}\space}}
\setlength{\LTcapwidth}{\linewidth}
\renewcommand{\arraystretch}{1.1}
\title{%s}
\author{%s}
\date{%s}
\begin{document}
\maketitle
\tableofcontents
''' % (title, subtitle, datetime.datetime.now().strftime('%Y-%m-%d %H:%M'))
        return head + '\n\n'.join(self.parts) + '\n\\end{document}\n'


# ----------------------------------------------------------------------
# Physics computations
# ----------------------------------------------------------------------

def coldelec(t, gain, shaping):
    '''Vectorized WireCell::Response::coldelec().'''
    t = np.asarray(t, dtype=float)
    out = np.zeros_like(t)
    ok = (t > 0) & (t < 10 * U.us)
    r = t[ok] / shaping
    g = gain * 10 * 1.012
    e1 = np.exp(-2.94809 * r)
    e2 = np.exp(-2.82833 * r)
    e3 = np.exp(-2.40318 * r)
    c1, s1 = np.cos(1.19361 * r), np.sin(1.19361 * r)
    c2, s2 = np.cos(2.38722 * r), np.sin(2.38722 * r)
    c3, s3 = np.cos(2.5928 * r), np.sin(2.5928 * r)
    c4, s4 = np.cos(5.18561 * r), np.sin(5.18561 * r)
    v = (4.31054 * e1 - 2.6202 * e2 * c1 - 2.6202 * e2 * c1 * c2 + 0.464924 * e3 * c3
         + 0.464924 * e3 * c3 * c4 + 0.762456 * e2 * s1 - 0.762456 * e2 * c2 * s1
         + 0.762456 * e2 * c1 * s2 - 2.620200 * e2 * s1 * s2 - 0.327684 * e3 * s3
         + 0.327684 * e3 * c4 * s3 - 0.327684 * e3 * c3 * s4 + 0.464924 * e3 * s3 * s4)
    out[ok] = v * g
    return out


class ElecResponse:
    '''Evaluate an electronics response component as a function of time.'''

    def __init__(self, comp, paths):
        self.comp = comp
        self.type = comp['type']
        self.data = comp.get('data', {})
        self.ok = True
        self.note = ''
        d = self.data
        self.postgain = d.get('postgain', 1.0)
        if self.type == 'ColdElecResponse':
            self.gain = d.get('gain', 14 * U.mV / U.fC)
            self.shaping = d.get('shaping', 2 * U.us)
        elif self.type == 'JsonElecResponse':
            fn = d.get('filename')
            p = paths.find_data(fn)
            if not p:
                self.ok = False
                self.note = 'file %s not found' % fn
                warn('JsonElecResponse file not found: %s' % fn)
            else:
                dat = load_json(p)
                self.times = np.array(dat['times'], dtype=float)
                self.amps = np.array(dat['amplitudes'], dtype=float)
                self.path = p
        else:
            self.ok = False
            self.note = 'response type %s not implemented in this script' % self.type

    def __call__(self, t):
        t = np.asarray(t, dtype=float)
        if self.type == 'ColdElecResponse':
            return coldelec(t, self.gain, self.shaping) * self.postgain
        if self.type == 'JsonElecResponse':
            v = np.interp(t, self.times, self.amps, left=0.0, right=0.0)
            return v * self.postgain
        return np.zeros_like(t)


class FieldResponse:

    def __init__(self, path):
        dat = load_json(path)
        fr = dat.get('FieldResponse', dat)
        self.path = path
        self.axis = fr.get('axis')
        self.origin = fr.get('origin')
        self.tstart = fr.get('tstart', 0.0)
        self.period = fr.get('period')
        self.speed = fr.get('speed')
        self.planes = []
        for p in fr['planes']:
            p = p.get('PlaneResponse', p)
            paths = []
            for pr in p['paths']:
                pr = pr.get('PathResponse', pr)
                cur = pr['current']
                if isinstance(cur, dict):
                    cur = cur.get('array', cur).get('elements')
                paths.append((pr['pitchpos'], np.array(cur, dtype=float)))
            paths.sort(key=lambda x: x[0])
            self.planes.append(dict(planeid=p['planeid'], location=p.get('location'),
                                    pitch=p['pitch'], paths=paths))

    def symmetric_paths(self, iplane):
        paths = self.planes[iplane]['paths']
        pp = [pr[0] for pr in paths]
        if min(pp) >= -1e-6:  # half-plane: mirror it
            mirrored = [(-x, c) for x, c in paths if x > 1e-6]
            paths = sorted(mirrored + paths, key=lambda x: x[0])
        return paths

    def wire_average(self, iplane, wire=0):
        '''Average current over impact positions in the region of a wire.'''
        pl = self.planes[iplane]
        pitch = pl['pitch']
        lo, hi = (wire - 0.5) * pitch, (wire + 0.5) * pitch
        tot = None
        wsum = 0.0
        for x, cur in self.symmetric_paths(iplane):
            if x < lo - 1e-6 or x > hi + 1e-6:
                continue
            w = 0.5 if (abs(x - lo) < 1e-6 or abs(x - hi) < 1e-6) else 1.0
            tot = cur * w if tot is None else tot + cur * w
            wsum += w
        if tot is None:
            return None
        return tot / wsum

    def nsamples(self):
        return len(self.planes[0]['paths'][0][1])


def rc_kernel(width, dt, n):
    t = np.arange(n) * dt
    h = -dt / width * np.exp(-t / width)
    h[0] += 1.0
    return h


def fft_convolve(a, b, n=None):
    n = n or (len(a) + len(b) - 1)
    nfft = 1 << int(math.ceil(math.log2(max(n, 2))))
    out = np.fft.irfft(np.fft.rfft(a, nfft) * np.fft.rfft(b, nfft), nfft)
    return out[:n]


def total_response(fr, iplane, er, rcs=(), long_padding=1.5 * U.ms, er_window=40 * U.us):
    '''FR (wire-0 average) convolved with electronics and RC responses.

    Returns (times, volts per electron) sampled at the FR period.'''
    dt = fr.period
    cur = fr.wire_average(iplane)
    if cur is None:
        return None, None
    ner = int(er_window / dt)
    ert = (np.arange(ner) + 0.5) * dt
    erw = er(ert)
    v = fft_convolve(cur * dt, erw)
    if rcs:
        n = len(v) + int(long_padding / dt)
        for width in rcs:
            v = fft_convolve(v, rc_kernel(width, dt, n), n)
    t = fr.tstart + np.arange(len(v)) * dt
    return t, v


def hf_filter(f, sigma, power, flag):
    f = np.asarray(f, dtype=float)
    with np.errstate(divide='ignore', invalid='ignore', over='ignore'):
        out = np.exp(-0.5 * np.power(np.abs(f) / sigma, power)) if sigma > 0 else np.ones_like(f)
    if flag:
        out = np.where(f == 0, 0.0, out)
    return out


def lf_filter(f, tau):
    f = np.asarray(f, dtype=float)
    if tau <= 0:
        return np.ones_like(f)
    return 1 - np.exp(-np.power(f / tau, 2))


PLANE_NAMES = {0: 'U', 1: 'V', 2: 'W'}


# ----------------------------------------------------------------------
# Wire geometry
# ----------------------------------------------------------------------

class WireGeom:

    def __init__(self, path):
        st = load_json(path)['Store']
        self.path = path
        P = np.array([[p['Point']['x'], p['Point']['y'], p['Point']['z']] for p in st['points']])
        self.anodes = {}
        wires = st['wires']
        for an in st['anodes']:
            an = an['Anode']
            faces = []
            for fi in an['faces']:
                face = st['faces'][fi]['Face']
                planes = []
                for pi in face['planes']:
                    pl = st['planes'][pi]['Plane']
                    ws = [wires[w]['Wire'] for w in pl['wires']]
                    tails = P[[w['tail'] for w in ws]]
                    heads = P[[w['head'] for w in ws]]
                    chans = np.array([w['channel'] for w in ws])
                    planes.append(self._plane_info(pl['ident'], tails, heads, chans))
                faces.append(dict(ident=face['ident'], planes=planes))
            self.anodes[an['ident']] = faces

    @staticmethod
    def _plane_info(ident, tails, heads, chans):
        info_ = dict(ident=ident, nwires=len(tails), tails=tails, heads=heads,
                     chmin=int(chans.min()) if len(chans) else -1,
                     chmax=int(chans.max()) if len(chans) else -1,
                     nchan=len(set(chans.tolist())))
        pts = np.vstack([tails, heads]) if len(tails) else np.zeros((1, 3))
        info_['x'] = float(np.mean(pts[:, 0]))
        info_['ymin'], info_['ymax'] = float(pts[:, 1].min()), float(pts[:, 1].max())
        info_['zmin'], info_['zmax'] = float(pts[:, 2].min()), float(pts[:, 2].max())
        pitch = angle = float('nan')
        n = len(tails)
        if n >= 3:
            ctr = 0.5 * (tails + heads)
            d = heads - tails
            dl = np.linalg.norm(d, axis=1)
            mid = n // 2
            u = d[mid] / dl[mid]
            pitches = []
            for i in range(max(0, mid - 5), min(n - 1, mid + 5)):
                dv = ctr[i + 1] - ctr[i]
                perp = dv - np.dot(dv, u) * u
                pitches.append(np.linalg.norm(perp))
            pitch = float(np.median(pitches))
            angle = math.degrees(math.atan2(u[1], u[2]))
            if angle > 90:
                angle -= 180
            if angle <= -90:
                angle += 180
        info_['pitch'] = pitch
        info_['angle'] = angle
        return info_


# ----------------------------------------------------------------------
# The summarizer
# ----------------------------------------------------------------------

class Summary:

    def __init__(self, args, paths, cfg, doc, figdir):
        self.args = args
        self.paths = paths
        self.cfg = cfg
        self.doc = doc
        self.figdir = figdir
        self._fr = {}
        self._wires = {}
        self.covered = set(self.COVERED)

    # -- data access

    def field(self, comp):
        if comp is None:
            return None
        fn = comp.get('data', {}).get('filename')
        if fn in self._fr:
            return self._fr[fn]
        p = self.paths.find_data(fn)
        fr = None
        if not p:
            warn('field response file not found: %s' % fn)
        else:
            try:
                fr = FieldResponse(p)
            except Exception as err:   # noqa: BLE001
                warn('failed to load field response %s: %s' % (p, err))
        self._fr[fn] = fr
        return fr

    def wires(self, fn):
        if fn in self._wires:
            return self._wires[fn]
        p = self.paths.find_data(fn)
        wg = None
        if not p:
            warn('wire schema file not found: %s' % fn)
        else:
            info('loading wires %s' % p)
            try:
                wg = WireGeom(p)
            except Exception as err:   # noqa: BLE001
                warn('failed to load wires %s: %s' % (p, err))
        self._wires[fn] = wg
        return wg

    def anode_of(self, comp):
        d = comp.get('data', {})
        for k in ('anode', 'anode_tn'):
            if isinstance(d.get(k), str):
                return self.cfg.get(d[k])
        return None

    # -- sections

    def run(self, jobinfo):
        sections = [
            ('Job', self.sec_job, jobinfo),
            ('Detector', self.sec_detector, None),
            ('LAr', self.sec_lar, None),
            ('Field', self.sec_fields, None),
            ('Elec', self.sec_elec, None),
            ('Sim', self.sec_sim, None),
            ('NF', self.sec_nf, None),
            ('SP', self.sec_sp, None),
            ('Other', self.sec_other, None),
            ('Appendix', self.sec_appendix, None),
        ]
        for name, meth, arg in sections:
            try:
                meth(arg) if arg is not None else meth()
            except Exception as err:   # noqa: BLE001
                import traceback
                traceback.print_exc()
                warn('section %s failed: %s' % (name, err))
                self.doc.para(r'\emph{Section %s failed: %s}' % (tex(name), tex(err)))
        self.sec_warnings()

    # Job ------------------------------------------------------------

    def sec_job(self, jobinfo):
        d = self.doc
        d.section('Job')
        d.para('This document summarizes the Wire-Cell Toolkit (WCT) configuration produced '
               'by compiling a Jsonnet file with the external variables given by a WCLS '
               'FHiCL configuration.  All quantities are given in WCT system of units '
               r'(mm, ns, MeV, $e^+$) followed, where known, by a conversion to more '
               'common units.')
        d.subsection('FHiCL source')
        rows = [
            ('FHiCL file', ttex(jobinfo['fcl']), ''),
            ('WCLS table', ttex(jobinfo['block']), 'FHiCL'),
        ]
        if jobinfo.get('alternatives'):
            rows.append(('other matching tables', ttex(', '.join(jobinfo['alternatives'])), 'FHiCL'))
        wm = jobinfo['wcls_main']
        for key in ('tool_type', 'apps', 'plugins', 'configs', 'inputers', 'outputers',
                    'logsinks', 'loglevels', 'paths'):
            if key in wm:
                v = fhicl_to_python(wm[key])
                rows.append((ttex('wcls_main.' + key), ttex(json.dumps(v)), 'FHiCL'))
        if jobinfo['block_table'].get('module_type') is not None:
            rows.insert(2, ('module_type', ttex(fhicl_to_python(jobinfo['block_table']['module_type'])), 'FHiCL'))
        d.kv_table(rows, 'WCLS configuration table found in FHiCL.')
        if jobinfo['missing_includes']:
            d.para(r'{\raggedright\textbf{FHiCL includes not found} (values from these are unresolved '
                   r'unless overridden): ' + ttex(', '.join(sorted(set(m[0] for m in jobinfo['missing_includes']))))
                   + r'\par}')

        d.subsection('Jsonnet external variables')
        d.para(r'The WCLS tool passes each \texttt{wcls\_main.params} entry as a Jsonnet '
               r'external \emph{string} and each \texttt{wcls\_main.structs} entry as '
               r'external \emph{code}.  The ``origin'' column gives where the value was '
               r'defined in FHiCL and the references through which it was obtained.')
        rows = []
        for ev in jobinfo['extvars'].values():
            rows.append((ttex(ev.name), ev.kind,
                         ttex(ev.value) if ev.value is not None else r'\textcolor{red}{unresolved}',
                         tex(ev.origin)))
        d.table(['Name', 'Kind', 'Value', 'Origin'], rows,
                ['0.2\\linewidth', '0.06\\linewidth', '0.3\\linewidth', '0.36\\linewidth'],
                'Jsonnet std.extVar() values.')
        d.subsection('Jsonnet compilation')
        d.para('Main Jsonnet file(s): ' + ', '.join(ttex(f) for f in jobinfo['jsonnet_files']))
        d.para('Command:')
        d.verbatim(' '.join(jobinfo['command']))

        d.subsection('Job types')
        jobs = detect_jobs(self.cfg)
        rows = [(title, ttex(', '.join(types))) for key, title, types in jobs]
        d.table(['Job type', 'Detected by component types'], rows,
                ['0.3\\linewidth', '0.64\\linewidth'],
                'Job types detected from the component types present in the configuration.')
        types = self.cfg.types()
        rows = [(ttex(t), str(n)) for t, n in sorted(types.items())]
        d.table(['Component type', 'Count'], rows, ['0.6\\linewidth', '0.2\\linewidth'],
                'All %d configured components by type.' % len(self.cfg.comps))
        self.dataflow()

    def dataflow(self):
        d = self.doc
        apps = self.cfg.apps()
        if not apps:
            return
        for app in apps:
            edges = app['data']['edges']
            d.subsection('Data flow graph: %s' % ttex(tn(app)))
            nodes = set()
            tedges = collections.Counter()
            for e in edges:
                t, h = e['tail']['node'], e['head']['node']
                nodes.add(t)
                nodes.add(h)
                tedges[(t.split(':')[0], h.split(':')[0])] += 1
            ntypes = collections.Counter(n.split(':')[0] for n in nodes)
            d.para('The application %s has %d edges between %d nodes of %d types.  The figure '
                   'collapses all nodes of one type into one box (with the number of '
                   'instances); edge labels count the collapsed edges.'
                   % (ttex(tn(app)), len(edges), len(nodes), len(ntypes)))
            if not shutil.which('dot'):
                d.para(r'\emph{graphviz dot not found, graph not drawn.}')
                continue
            lines = ['digraph G {', 'rankdir=TB;', 'node [shape=box, fontsize=10, style=rounded];',
                     'edge [fontsize=8];']
            ids = {t: 'n%d' % i for i, t in enumerate(sorted(ntypes))}
            for t, n in ntypes.items():
                lines.append('%s [label="%s\\n×%d"];' % (ids[t], t, n))
            for (a, b), n in tedges.items():
                lines.append('%s -> %s [label="%d"];' % (ids[a], ids[b], n))
            lines.append('}')
            self.doc.nfig += 1
            base = os.path.join(self.figdir, 'graph%03d' % self.doc.nfig)
            open(base + '.dot', 'w').write('\n'.join(lines))
            r = subprocess.run(['dot', '-Tpdf', '-o', base + '.pdf', base + '.dot'],
                               stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
            if r.returncode == 0:
                d.include_pdf(os.path.basename(base) + '.pdf',
                              'Type-collapsed data flow graph of %s.' % ttex(tn(app)))
            else:
                warn('dot failed: %s' % r.stderr)

    # Detector --------------------------------------------------------

    def sec_detector(self):
        d = self.doc
        anodes = self.cfg.of_type('AnodePlane')
        if not anodes:
            return
        d.section('Detector')
        d.para(r'The detector is described by %d \texttt{AnodePlane} components.  Each '
               r'anode has up to two faces.  For each face the configuration gives the '
               r'$x$ position of the anode plane (\texttt{anode}), of the plane where '
               r'the field response calculation starts (\texttt{response}) and of the '
               r'cathode (\texttt{cathode}).  Wire (strip) geometry comes from the '
               r'wire schema file of the \texttt{WireSchemaFile} component.' % len(anodes))
        d.para(r'\textbf{Coordinate system:} WCT uses the global coordinates of the wire '
               r'schema file.  The $x$ axis is the drift axis (field responses are '
               r'calculated along their \texttt{axis}, normally $(1,0,0)$) and wires lie '
               r'in $y$--$z$ planes.  In the DUNE/LArSoft convention $y$ is vertical (up) '
               r'and $z$ is along the beam.  Electrons drift from the cathode toward the anode.')

        schemas = collections.OrderedDict()
        for a in anodes:
            ws = self.cfg.get(a['data'].get('wire_schema'))
            fn = ws['data'].get('filename') if ws else None
            schemas.setdefault(fn, []).append(a)
        rows = []
        for fn, ans in schemas.items():
            ws = self.cfg.get(ans[0]['data'].get('wire_schema'))
            rows.append((ttex('filename'), ttex(fn), ttex(tn(ws)) if ws else '?'))
            p = self.paths.find_data(fn)
            rows.append((ttex('(found at)'), ttex(p or 'NOT FOUND'), ''))
        rows += [(ttex(k), fmt_value(k, v, 'AnodePlane'), ttex('AnodePlane') + ' (%s)' % src)
                 for k, vals in group_params(anodes, skip=(r'^faces', r'^ident$', r'^wire_schema$'))
                 for v, names in vals for src in [('all' if len(names) == len(anodes) else names_tex(names))]]
        d.kv_table(rows, 'Wire schema and common anode parameters.')

        # per anode geometry
        geo = []   # (anode comp, face index, face cfg, face geom or None)
        for fn, ans in schemas.items():
            wg = self.wires(fn) if fn else None
            for a in ans:
                ident = a['data'].get('ident')
                faces = a['data'].get('faces', [])
                fgeo = wg.anodes.get(ident) if wg else None
                for fi, fc in enumerate(faces):
                    if fc is None:
                        continue
                    fg = fgeo[fi] if fgeo and fi < len(fgeo) else None
                    geo.append((a, fi, fc, fg))
        self.geo = geo
        self.det_totals(geo, anodes)
        self.det_checks(geo)
        self.det_figures(geo)
        self.det_assignments(anodes)
        self.det_table(geo)

    def det_totals(self, geo, anodes):
        d = self.doc
        nwires = nchan = 0
        allpts = []
        for a, fi, fc, fg in geo:
            if fg:
                for pl in fg['planes']:
                    nwires += pl['nwires']
                    allpts.append([pl['ymin'], pl['ymax'], pl['zmin'], pl['zmax']])
        chans = set()
        for fn, wg in self._wires.items():
            if not wg:
                continue
            for faces in wg.anodes.values():
                for f in faces:
                    for pl in f['planes']:
                        chans.update(range(pl['chmin'], pl['chmax'] + 1))
        nchan = len(chans)
        xs = [fc[k] for a, fi, fc, fg in geo for k in ('anode', 'cathode')]
        rows = [('number of anodes (TPCs)', str(len(anodes)), ttex('AnodePlane')),
                ('number of configured faces', str(len(geo)), ttex('AnodePlane.data.faces'))]
        if allpts:
            arr = np.array(allpts)
            rows += [('number of wires (all faces)', str(nwires), 'wire schema'),
                     ('channel range span', '%d (upper bound on channel count)' % nchan, 'wire schema'),
                     ('y extent', '%.1f to %.1f cm (%.1f cm)' % (arr[:, 0].min() / U.cm, arr[:, 1].max() / U.cm,
                                                               (arr[:, 1].max() - arr[:, 0].min()) / U.cm), 'wire schema'),
                     ('z extent', '%.1f to %.1f cm (%.1f cm)' % (arr[:, 2].min() / U.cm, arr[:, 3].max() / U.cm,
                                                               (arr[:, 3].max() - arr[:, 2].min()) / U.cm), 'wire schema')]
        if xs:
            rows.append(('x extent (anode/cathode)', '%.1f to %.1f cm (%.1f cm)'
                         % (min(xs) / U.cm, max(xs) / U.cm, (max(xs) - min(xs)) / U.cm), ttex('AnodePlane.data.faces')))
        d.subsection('Overall dimensions')
        d.kv_table(rows, 'Overall detector dimensions.')

    def det_checks(self, geo):
        '''Consistency of face x positions and field response origin.'''
        d = self.doc
        # field responses used by each anode (sim and sp)
        fr_of = collections.defaultdict(set)
        for c in self.cfg.of_type('DepoTransform', 'Ductor', 'DepoSplat', 'DepoZipper'):
            a = self.anode_of(c)
            if not a:
                continue
            for p in c['data'].get('pirs', []):
                pir = self.cfg.get(p)
                if pir:
                    fr_of[tn(a)].add(pir['data'].get('field_response'))
        for c in self.cfg.of_type('OmnibusSigProc'):
            a = self.anode_of(c)
            if a:
                fr_of[tn(a)].add(c['data'].get('field_response'))
        groups = collections.OrderedDict()
        for a, fi, fc, fg in geo:
            if not fg or not fg['planes']:
                continue
            xcol = fg['planes'][-1]['x']
            key = (round(fc['anode'] - xcol, 1), round(fc['response'] - xcol, 1),
                   round(fc['cathode'] - fc['anode'], 1), tuple(sorted(x for x in fr_of.get(tn(a), []) if x)))
            groups.setdefault(key, []).append('%s/face%d' % (tn(a), fi))
        if not groups:
            return
        d.subsection('Consistency checks')
        d.para(r'For each distinct face configuration: $x$ distances relative to the last '
               r'(collection) wire plane of the face, the drift length (cathode--anode) and, '
               r'for comparison, the \texttt{origin} of the field response(s) used by the '
               r'anode.  The magnitude of the response plane distance should equal the '
               r'field response origin.')
        rows = []
        for (da, dr, dl, frs), names in groups.items():
            frtxt = []
            for f in frs:
                fr = self.field(self.cfg.get(f))
                if fr:
                    delta = abs(dr) - fr.origin
                    ok = abs(delta) < 0.05 * U.mm
                    color = 'orange' if abs(delta) < 2 * U.mm else 'red'
                    frtxt.append(r'%s origin=%.2f cm: %s' % (
                        tex(f), fr.origin / U.cm, r'\textcolor{green!50!black}{OK}' if ok else
                        r'\textcolor{%s}{differs by %+.2f mm}' % (color, delta / U.mm)))
            rows.append((names_tex(names, 3), '%.2f' % (da / U.cm), '%.2f' % (dr / U.cm), '%.2f' % (dl / U.cm),
                         '; '.join(frtxt) or '-'))
        d.table(['Faces', 'anode-coll. [cm]', 'response-coll. [cm]', 'drift length [cm]', 'Field response check'],
                rows, ['0.24\\linewidth', '0.1\\linewidth', '0.11\\linewidth', '0.1\\linewidth', '0.33\\linewidth'],
                'Face geometry consistency.')

    def det_figures(self, geo):
        d = self.doc
        withg = [(a, fi, fc, fg) for a, fi, fc, fg in geo if fg and fg['planes']]
        if not withg:
            return
        d.subsection('Detector drawings')
        # group faces into panels by anode x and drift direction
        panels = collections.OrderedDict()
        for a, fi, fc, fg in withg:
            sign = 1 if fc['cathode'] > fc['anode'] else -1
            key = (round(fc['anode'] / (10 * U.cm)), sign)
            panels.setdefault(key, []).append((a, fi, fc, fg))
        for (xk, sign), faces in panels.items():
            self.det_panel(faces, sign)

        # x-z and x-y views
        fig, axes = plt.subplots(1, 2, figsize=(14, 6))
        for ax, (j, lab) in zip(axes, ((2, 'z'), (1, 'y'))):
            for a, fi, fc, fg in withg:
                pl = fg['planes']
                lo = min(p[lab + 'min'] for p in pl) / U.m
                hi = max(p[lab + 'max'] for p in pl) / U.m
                xa, xr, xc = fc['anode'] / U.m, fc['response'] / U.m, fc['cathode'] / U.m
                ax.add_patch(Rectangle((min(xa, xc), lo), abs(xc - xa), hi - lo, color='C0',
                                       alpha=min(0.1, 1.0 / len(withg)), lw=0))
                ax.plot([xa, xa], [lo, hi], 'C0-', lw=1.2)
                ax.plot([xr, xr], [lo, hi], 'C2--', lw=0.6)
                ax.plot([xc, xc], [lo, hi], 'C3:', lw=1.0)
            ax.plot([], [], 'C0-', label='anode')
            ax.plot([], [], 'C2--', label='response plane')
            ax.plot([], [], 'C3:', label='cathode')
            ax.set_xlabel('x [m] (drift axis)')
            ax.set_ylabel('%s [m]' % lab)
            ax.legend(fontsize=8, loc='best')
            ax.set_title('%s vs x' % lab)
            ax.autoscale_view()
        d.figure(fig, r'Drift volumes seen from above ($z$ vs $x$) and from the side ($y$ vs $x$).  Shaded: '
                 r'volume between anode and cathode of each face.  Source: \texttt{AnodePlane.data.faces}.')

        # wire detail for first anode
        a, fi, fc, fg = withg[0]
        pl = fg['planes']
        fig, axes = plt.subplots(1, len(pl), figsize=(5 * len(pl), 5), squeeze=False)
        for ax, p in zip(axes[0], pl):
            n = p['nwires']
            step = max(1, n // 40)
            for i in range(0, n, step):
                t, h = p['tails'][i] / U.cm, p['heads'][i] / U.cm
                ax.plot([t[2], h[2]], [t[1], h[1]], lw=0.6, color='C%d' % (p['ident'] % 10))
            ax.set_aspect('equal')
            ax.set_xlabel('z [cm]')
            ax.set_ylabel('y [cm]')
            ax.set_title('plane %d: %d wires, pitch %.3f mm, angle %.1f°\nch %d-%d, x=%.2f cm (every %d. wire drawn)'
                         % (p['ident'], n, p['pitch'] / U.mm, p['angle'], p['chmin'], p['chmax'], p['x'] / U.cm, step),
                         fontsize=8)
        d.figure(fig, r'Wires of each plane of face %d of %s (angle measured from the $z$ axis toward $y$).  '
                 r'Source: wire schema, \texttt{AnodePlane} ident %s.' % (fi, ttex(tn(a)), a['data'].get('ident')))

        # plane geometry table of distinct plane layouts
        layouts = collections.OrderedDict()
        for a, fi, fc, fg in withg:
            key = tuple((p['ident'], p['nwires'], round(p['pitch'], 2), round(p['angle'], 1)) for p in fg['planes'])
            layouts.setdefault(key, []).append('%s/face%d' % (tn(a), fi))
        rows = []
        for key, names in layouts.items():
            for i, (pid, nw, pitch, ang) in enumerate(key):
                rows.append((names_tex(names, 3) if i == 0 else '', str(pid), str(nw),
                             '%.3f' % pitch, '%.1f' % ang))
        d.table(['Faces', 'Plane', 'Wires', 'Pitch [mm]', 'Angle [deg]'], rows,
                ['0.4\\linewidth', '0.08\\linewidth', '0.1\\linewidth', '0.12\\linewidth', '0.12\\linewidth'],
                'Distinct wire plane layouts (from the wire schema).')

    def det_panel(self, faces, sign):
        '''Draw anode faces of one panel in the y-z plane, labeled by ident.'''
        d = self.doc
        rects = []
        anodes = collections.OrderedDict()
        for a, fi, fc, fg in faces:
            pl = fg['planes']
            r = (min(p['zmin'] for p in pl) / U.m, max(p['zmax'] for p in pl) / U.m,
                 min(p['ymin'] for p in pl) / U.m, max(p['ymax'] for p in pl) / U.m)
            rects.append(r)
            u = anodes.get(tn(a))
            anodes[tn(a)] = (a, r if u is None else (min(u[1][0], r[0]), max(u[1][1], r[1]),
                                                      min(u[1][2], r[2]), max(u[1][3], r[3])))
        z0 = min(r[0] for r in rects)
        z1 = max(r[1] for r in rects)
        y0 = min(r[2] for r in rects)
        y1 = max(r[3] for r in rects)
        zspan, yspan = z1 - z0, max(y1 - y0, 1e-3)
        # split long detectors into several rows of aspect ratio ~2:1
        nrow = max(1, min(6, int(round(zspan / (2.5 * yspan)))))
        seg = zspan / nrow
        width_in = 10.0
        row_h = width_in * yspan / seg
        if row_h * nrow > 8.5:  # tall detector: shrink the width to fit the page
            scale = 8.5 / (row_h * nrow)
            width_in *= scale
            row_h *= scale
        fig, axes = plt.subplots(nrow, 1, figsize=(max(width_in, 3.0), row_h * nrow + 0.8), squeeze=False)
        # font size from the typical anode box width
        bw = np.median([u[1][1] - u[1][0] for u in anodes.values()]) / seg * width_in * 72
        bh = np.median([u[1][3] - u[1][2] for u in anodes.values()]) / yspan * row_h * 72
        nchar = max(len(str(u[0]['data'].get('ident'))) for u in anodes.values())
        fsize = max(3, min(12, bw / (0.65 * nchar), bh / 1.5))
        for irow, ax in enumerate(axes[:, 0]):
            for r in rects:
                ax.add_patch(Rectangle((r[0], r[2]), r[1] - r[0], r[3] - r[2], fill=False, lw=0.3, color='gray'))
            for a, u in anodes.values():
                ax.add_patch(Rectangle((u[0], u[2]), u[1] - u[0], u[3] - u[2], fill=False, lw=0.8, color='k'))
                zc = 0.5 * (u[0] + u[1])
                if z0 + irow * seg - 1e-6 <= zc <= z0 + (irow + 1) * seg + 1e-6:
                    ax.text(zc, 0.5 * (u[2] + u[3]), str(a['data'].get('ident')),
                            ha='center', va='center', fontsize=fsize, color='C3')
            ax.set_xlim(z0 + irow * seg - 0.01 * seg, z0 + (irow + 1) * seg + 0.01 * seg)
            ax.set_ylim(y0 - 0.02 * yspan, y1 + 0.02 * yspan)
            ax.set_aspect('equal')
            ax.set_ylabel('y [m]')
            ax.tick_params(labelsize=7)
        axes[-1, 0].set_xlabel('z [m]')
        xs = faces[0][2]['anode'] / U.cm
        axes[0, 0].set_title('%d faces of %d anodes with anode plane at x=%.1f cm, electrons drift toward %sx; '
                             'labels: AnodePlane ident' % (len(faces), len(anodes), xs, '-' if sign > 0 else '+'),
                             fontsize=8)
        d.figure(fig, r'Anodes seen along the drift ($x$) axis with anode plane at $x=%.1f$\,cm.  Black: '
                 r'$y$--$z$ extent of all wires of each anode (labeled with the \texttt{AnodePlane} '
                 r'\texttt{ident}, the per-TPC number); gray: individual faces.  Overall extent '
                 r'$%.2f\times%.2f$\,m$^2$ ($z\times y$).  Source: \texttt{AnodePlane} components and wire schema.'
                 % (xs, zspan, yspan))

    def det_assignments(self, anodes):
        '''Which responses, noise and digitizer settings apply to each anode.'''
        d = self.doc
        per = collections.defaultdict(dict)
        for c in self.cfg.of_type('DepoTransform', 'Ductor', 'DepoSplat', 'DepoZipper'):
            a = self.anode_of(c)
            if not a:
                continue
            frs, ers = set(), set()
            for p in c['data'].get('pirs', []):
                pir = self.cfg.get(p)
                if pir:
                    frs.add(pir['data'].get('field_response'))
                    ers.update(pir['data'].get('short_responses', []))
                    nrc = len(pir['data'].get('long_responses', []))
            per[tn(a)]['sim field response'] = ', '.join(sorted(x for x in frs if x))
            per[tn(a)]['sim short responses'] = ', '.join(sorted(ers))
            per[tn(a)]['sim RC responses'] = str(nrc)
        for c in self.cfg.of_type('EmpiricalNoiseModel'):
            a = self.anode_of(c)
            if a:
                per[tn(a)]['noise spectra'] = c['data'].get('spectra_file', '')
        for c in self.cfg.of_type('Digitizer'):
            a = self.anode_of(c)
            if a:
                dd = c['data']
                per[tn(a)]['digitizer'] = '%s bit, baselines %s mV, fullscale %s mV' % (
                    dd.get('resolution'), [num(x / U.mV, 5) for x in dd.get('baselines', [])],
                    [num(x / U.mV, 5) for x in dd.get('fullscale', [])])
        for c in self.cfg.of_type('OmnibusSigProc'):
            a = self.anode_of(c)
            if a:
                per[tn(a)]['SP field response'] = c['data'].get('field_response', '')
                per[tn(a)]['SP elec response'] = c['data'].get('elecresponse', '')
        if not per:
            return
        groups = collections.OrderedDict()
        for a in anodes:
            key = tuple(sorted(per.get(tn(a), {}).items()))
            groups.setdefault(key, []).append(tn(a))
        d.subsection('Per-anode response and noise assignment')
        d.para('Anodes grouped by the responses, noise spectra and digitizer settings that '
               'the configuration attaches to them (via the \\texttt{anode} parameter of '
               'the simulation, noise, digitizer and signal processing components).')
        rows = []
        for key, names in groups.items():
            for i, (k, v) in enumerate(key):
                rows.append((names_tex(names, 3) if i == 0 else '', tex(k), ttex(v)))
            rows.append(('', '', ''))
        d.table(['Anodes', 'Aspect', 'Value'], rows,
                ['0.28\\linewidth', '0.2\\linewidth', '0.46\\linewidth'],
                'Per-anode assignment of responses, noise and digitization.')

    def det_table(self, geo):
        '''One row per anode.  Faces with identical x positions are merged.'''
        d = self.doc
        per = collections.OrderedDict()
        for a, fi, fc, fg in geo:
            per.setdefault(tn(a), (a, []))[1].append((fi, fc, fg))
        rows = []
        maxrows = self.args.max_anode_rows
        for a, faces in list(per.values())[:maxrows]:
            pls = [p for fi, fc, fg in faces if fg for p in fg['planes']]
            if pls:
                yr = '%.1f / %.1f' % (min(p['ymin'] for p in pls) / U.cm, max(p['ymax'] for p in pls) / U.cm)
                zr = '%.1f / %.1f' % (min(p['zmin'] for p in pls) / U.cm, max(p['zmax'] for p in pls) / U.cm)
                ch = '%d-%d' % (min(p['chmin'] for p in pls), max(p['chmax'] for p in pls))
            else:
                yr = zr = ch = '-'
            xs = collections.OrderedDict()
            for fi, fc, fg in faces:
                key = tuple('%.2f' % (fc[k] / U.cm) for k in ('anode', 'response', 'cathode'))
                xs.setdefault(key, []).append(str(fi))
            xcols = [' / '.join(k[i] for k in xs) for i in range(3)]
            rows.append((str(a['data'].get('ident')), ttex(a.get('name', '')), ','.join(str(f[0]) for f in faces),
                         xcols[0], xcols[1], xcols[2], yr, zr, ch))
        cap = ('Per anode geometry: x positions of anode plane, response plane and cathode from '
               'AnodePlane faces (one value per distinct face), y/z ranges and channel range from the wire schema.')
        if len(per) > maxrows:
            cap += '  Only the first %d of %d anodes are listed (see --max-anode-rows).' % (maxrows, len(per))
        d.subsection('Anodes')
        d.table(['ident', 'name', 'faces', 'x anode [cm]', 'x resp. [cm]', 'x cath. [cm]',
                 'y min/max [cm]', 'z min/max [cm]', 'channels'], rows,
                ['0.05\\linewidth', '0.1\\linewidth', '0.04\\linewidth', '0.08\\linewidth', '0.08\\linewidth',
                 '0.08\\linewidth', '0.14\\linewidth', '0.14\\linewidth', '0.12\\linewidth'], cap,
                fontsize=r'\scriptsize')

    # LAr / physics constants ----------------------------------------

    def sec_lar(self):
        d = self.doc
        drifters = self.cfg.of_type('Drifter')
        others = self.cfg.of_type('DepoSetDrifter', 'DepoTransform', 'Ductor', 'DepoSplat', 'DepoZipper')
        if not drifters and not others:
            return
        d.section('Physics constants')
        d.para('Liquid argon transport parameters as used by the simulation.  '
               r'The \texttt{Drifter} applies diffusion (DL, DT), attachment (lifetime) and '
               r'fluctuations while drifting depositions to the response plane.  The '
               r'drift speed is also used by the components that convert drifted depositions '
               r'into signals.')
        rows = []
        for key in ('DL', 'DT', 'lifetime', 'drift_speed', 'fluctuate', 'time_offset', 'charge_scale'):
            for comps in (drifters, others):
                for k, vals in group_params(comps):
                    if k != key:
                        continue
                    for v, names in vals:
                        rows.append((ttex(key), fmt_value(key, v, comps[0]['type']), names_tex(names, 3)))
        d.kv_table(rows, 'Physics constants.')
        frs = {}
        for c in self.cfg.of_type('FieldResponse'):
            fr = self.field(c)
            if fr:
                frs[tn(c)] = fr
        if frs:
            rows = [(ttex(k), r'%s mm/$\mu$s' % num(fr.speed / (U.mm / U.us), 5), ttex(os.path.basename(fr.path)))
                    for k, fr in frs.items()]
            d.table(['Field response', 'Drift speed in FR calculation', 'File'], rows,
                    ['0.3\\linewidth', '0.25\\linewidth', '0.39\\linewidth'],
                    'Drift speed assumed by field response calculations (for comparison).')
        for c in drifters:
            xr = c['data'].get('xregions', [])
            if xr:
                rows = [(str(i), '%.2f' % (r.get('anode', float('nan')) / U.cm),
                         '%.2f' % (r.get('response', float('nan')) / U.cm),
                         '%.2f' % (r.get('cathode', float('nan')) / U.cm)) for i, r in enumerate(xr) if r]
                d.table(['Region', 'anode x [cm]', 'response x [cm]', 'cathode x [cm]'], rows,
                        ['0.1\\linewidth'] + ['0.2\\linewidth'] * 3,
                        r'Drift regions (\texttt{xregions}) of %s.  Depositions outside are dropped; '
                        r'those inside are drifted to the response plane.' % ttex(tn(c)))

    # Field responses -------------------------------------------------

    def sec_fields(self):
        d = self.doc
        frcs = self.cfg.of_type('FieldResponse')
        if not frcs:
            return
        d.section('Field responses')
        d.para(r'Each \texttt{FieldResponse} component loads a file of induced currents '
               r'on the wires of each plane for electrons drifting from the response plane '
               r'(at distance \texttt{origin} from the last plane) along paths at a set of '
               r'transverse impact positions (\texttt{pitchpos}).')
        seen = {}
        for c in frcs:
            fn = c['data'].get('filename')
            seen.setdefault(fn, []).append(tn(c))
        for fn, names in seen.items():
            fr = self.field(self.cfg.get(names[0]))
            d.subsection('Field response %s' % ttex(fn))
            if not fr:
                d.para(r'\emph{File not found.}')
                continue
            rows = [('filename', ttex(fn), names_tex(names)), ('path', ttex(fr.path), '')]
            for k in ('axis', 'origin', 'tstart', 'period', 'speed'):
                rows.append((k, fmt_value(k, getattr(fr, k)), 'file'))
            rows.append(('samples per path', str(fr.nsamples()), 'file'))
            for i, p in enumerate(fr.planes):
                pp = [x for x, _ in p['paths']]
                ch = fr.wire_average(i)
                q = ch.sum() * fr.period if ch is not None else float('nan')
                rows.append(('plane %d' % i, 'planeid=%s, location=%s, pitch=%s, %d paths, pitchpos %.2f to %.2f mm, '
                             'wire-0 integrated charge %.3f e' % (p['planeid'], fmt_value('location', p['location']),
                                                                  fmt_value('pitch', p['pitch']), len(pp), min(pp), max(pp), q),
                             'file'))
            d.kv_table(rows, 'Field response %s.' % ttex(fn))
            self.plot_field(fr, fn, names)

    def plot_field(self, fr, fn, names):
        d = self.doc
        np_ = len(fr.planes)
        fig, axes = plt.subplots(1, np_, figsize=(5 * np_, 4), squeeze=False)
        for i, ax in enumerate(axes[0]):
            t = (fr.tstart + np.arange(fr.nsamples()) * fr.period) / U.us
            for w in (0, 1, 2):
                cur = fr.wire_average(i, w)
                if cur is None:
                    continue
                ax.plot(t, cur * U.us, label='wire %+d' % w if w else 'wire 0', lw=1.2 if w == 0 else 0.8)
            ax.set_title('plane %d (%s)' % (fr.planes[i]['planeid'], PLANE_NAMES.get(fr.planes[i]['planeid'], '?')))
            ax.set_xlabel(r'time [$\mu$s]')
            ax.set_ylabel(r'current [$e^+$/$\mu$s] per unit charge')
            ax.legend(fontsize=8)
            # zoom on the interesting part
            cur0 = fr.wire_average(i, 0)
            if cur0 is not None:
                big = np.where(np.abs(cur0) > 0.01 * np.abs(cur0).max())[0]
                if len(big):
                    ax.set_xlim(t[max(0, big[0] - 20)], t[min(len(t) - 1, big[-1] + 40)])
        d.figure(fig, r'Field response %s: induced current, as stored, averaged over impact positions within the '
                 r'region of the central wire (wire 0) and its neighbors.  WCT multiplies this by the charge of '
                 r'each deposition (negative for ionization electrons).  Source: %s.'
                 % (ttex(fn), names_tex(names, 3)))
        fig, axes = plt.subplots(1, np_, figsize=(5 * np_, 4), squeeze=False)
        for i, ax in enumerate(axes[0]):
            paths = fr.symmetric_paths(i)
            pp = np.array([x for x, _ in paths])
            img = np.array([c for _, c in paths]) * U.us
            t = (fr.tstart + np.arange(img.shape[1]) * fr.period) / U.us
            vmax = np.abs(img).max()
            norm = mcolors.SymLogNorm(linthresh=vmax * 1e-3, vmin=-vmax, vmax=vmax)
            im = ax.pcolormesh(t, pp, img, cmap='RdBu_r', norm=norm, shading='auto', rasterized=True)
            cur0 = fr.wire_average(i, 0)
            big = np.where(np.abs(cur0) > 0.01 * np.abs(cur0).max())[0] if cur0 is not None else []
            if len(big):
                ax.set_xlim(t[max(0, big[0] - 20)], t[min(len(t) - 1, big[-1] + 40)])
            ax.set_xlabel(r'time [$\mu$s]')
            ax.set_ylabel('impact position [mm]')
            ax.set_title('plane %d' % fr.planes[i]['planeid'])
            for w in range(-12, 13):
                y = (w + 0.5) * fr.planes[i]['pitch']
                if pp.min() <= y <= pp.max():
                    ax.axhline(y, color='gray', lw=0.2)
            fig.colorbar(im, ax=ax, label=r'e/$\mu$s')
        d.figure(fig, r'Field response %s: current vs time and impact position for all paths '
                 r'(symmetric log color scale; gray lines are wire region boundaries).' % ttex(fn))

    # Electronics / total responses ---------------------------------

    def sec_elec(self):
        d = self.doc
        etypes = ('ColdElecResponse', 'JsonElecResponse', 'WarmElecResponse', 'ParamsColdElecResponse')
        ercs = self.cfg.of_type(*etypes)
        rcs = self.cfg.of_type('RCResponse')
        pirs = self.cfg.of_type('PlaneImpactResponse')
        if not ercs and not pirs:
            return
        d.section('Electronics and total responses')
        d.para(r'Electronics response components give the response of the front end '
               r'electronics to a unit charge.  The simulation combines a field response '
               r'with ``short'' (electronics) and ``long'' (RC) responses in '
               r'\texttt{PlaneImpactResponse} components; signal processing uses the field '
               r'and electronics responses named in its configuration.')
        for c in ercs:
            d.group_table([c])
        if rcs:
            d.group_table(rcs)
        if ercs:
            fig, ax = plt.subplots(figsize=(8, 4.5))
            t = np.linspace(0, 20 * U.us, 2001)
            for c in ercs:
                er = ElecResponse(c, self.paths)
                if not er.ok:
                    d.para(r'\emph{%s: %s}' % (ttex(tn(c)), tex(er.note)))
                    continue
                y = er(t) / (U.mV / U.fC)
                ax.plot(t / U.us, y, label='%s (peak %.2f mV/fC)' % (tn(c), y.max()))
            ax.set_xlabel(r'time [$\mu$s]')
            ax.set_ylabel('response [mV/fC] (incl. postgain)')
            ax.legend(fontsize=8)
            ax.grid(alpha=0.3)
            d.figure(fig, 'Electronics responses (including postgain).  Source: %s.'
                     % names_tex([tn(c) for c in ercs], 6))
        if rcs:
            fig, ax = plt.subplots(figsize=(8, 3.5))
            for c in rcs:
                w = c['data'].get('width', 1 * U.ms)
                tt = np.linspace(0, 5 * w, 500)
                ax.plot(tt / U.ms, -np.exp(-tt / w) / (w / U.us), label='%s: width %.3g ms (tail per $\\mu$s)' % (tn(c), w / U.ms))
            ax.set_xlabel('time [ms]')
            ax.set_ylabel(r'tail amplitude [1/$\mu$s]')
            ax.legend(fontsize=8)
            nuse = collections.Counter()
            for p in pirs:
                nuse[len(p['data'].get('long_responses', []))] += 1
            use = ', '.join('%d PIRs apply %d RC' % (n, k) for k, n in sorted(nuse.items()))
            d.figure(fig, r'RC responses: a delta function followed by a negative exponential tail. '
                     r'The number of RC responses applied by each \texttt{PlaneImpactResponse} is given by '
                     r'its \texttt{long\_responses}: %s.' % tex(use or 'no PIRs'))
        if pirs:
            d.subsection('Plane impact responses (simulation)')
            d.group_table(pirs)
            self.plot_total_sim(pirs)

    def adc_per_volt_for(self, anode_tns):
        for c in self.cfg.of_type('Digitizer'):
            a = self.anode_of(c)
            if a and tn(a) in anode_tns:
                dd = c['data']
                fs = dd.get('fullscale', [0, 2 * U.volt])
                res = dd.get('resolution', 12)
                return dd.get('gain', 1.0) * ((1 << int(res)) - 1) / (fs[1] - fs[0]), tn(c)
        return None, None

    def plot_total_sim(self, pirs):
        d = self.doc
        # group PIR trios by (field, shorts, longs)
        groups = collections.OrderedDict()
        for p in pirs:
            dd = p['data']
            key = (dd.get('field_response'), tuple(dd.get('short_responses', [])),
                   tuple(dd.get('long_responses', [])))
            groups.setdefault(key, []).append(p)
        # anodes using each group
        users = collections.defaultdict(list)
        for c in self.cfg.of_type('DepoTransform', 'Ductor', 'DepoSplat', 'DepoZipper'):
            a = self.anode_of(c)
            for p in c['data'].get('pirs', []):
                pc = self.cfg.get(p)
                if pc is None or a is None:
                    continue
                dd = pc['data']
                key = (dd.get('field_response'), tuple(dd.get('short_responses', [])),
                       tuple(dd.get('long_responses', [])))
                users[key].append(tn(a))
        nelec = 1e4
        for key, ps in groups.items():
            frn, shorts, longs = key
            fr = self.field(self.cfg.get(frn))
            if not fr:
                continue
            ers = []
            for s in shorts:
                sc = self.cfg.get(s)
                if sc and sc['type'] in ('ColdElecResponse', 'JsonElecResponse'):
                    ers.append(ElecResponse(sc, self.paths))
            if len(ers) != 1 or not ers[0].ok:
                d.para(r'\emph{Total response for %s not drawn (need exactly one supported electronics short response).}'
                       % names_tex([tn(p) for p in ps]))
                continue
            widths = [self.cfg.get(l)['data'].get('width', U.ms) for l in longs if self.cfg.get(l)]
            anode_tns = sorted(set(users.get(key, [])))
            adcpv, digi = self.adc_per_volt_for(set(anode_tns))
            tick = ps[0]['data'].get('tick', 0.5 * U.us)
            fig, axes = plt.subplots(1, len(fr.planes), figsize=(5 * len(fr.planes), 4), squeeze=False)
            for i, ax in enumerate(axes[0]):
                t, v = total_response(fr, i, ers[0], widths)
                if t is None:
                    continue
                y = -v * nelec  # ionization electrons carry negative charge
                ylab = 'mV'
                if adcpv:
                    y = y * adcpv
                    ylab = 'ADC'
                else:
                    y = y / U.mV
                ax.plot(t / U.us, y, lw=1, label='FR period sampling')
                step = max(1, int(round(tick / fr.period)))
                ax.plot(t[::step] / U.us, y[::step], '.', ms=3, label='tick %.2f $\\mu$s' % (tick / U.us))
                big = np.where(np.abs(y) > 0.01 * np.abs(y).max())[0]
                if len(big):
                    ax.set_xlim(t[big[0]] / U.us - 5, t[min(len(t) - 1, big[0] + int(40 * U.us / fr.period))] / U.us)
                ax.set_title('plane %d: peak %+.1f / %+.1f %s' % (fr.planes[i]['planeid'], y.max(), y.min(), ylab), fontsize=9)
                ax.set_xlabel(r'time [$\mu$s]')
                ax.set_ylabel('%s per %d electrons' % (ylab, nelec))
                ax.legend(fontsize=7)
                ax.grid(alpha=0.3)
            cap = (r'Total response (field response wire-0 average $\otimes$ %s%s) for %d ionization electrons '
                   r'(charge $-%d\,e$) at the response plane, without diffusion.  '
                   r'PIRs: %s.  Field: %s.  '
                   % (ttex(shorts[0]), (r' $\otimes$ %d RC' % len(widths)) if widths else '', nelec, nelec,
                      names_tex([tn(p) for p in ps], 3), ttex(frn)))
            if adcpv:
                cap += r'Converted to ADC with %s.  ' % ttex(digi)
            if anode_tns:
                cap += 'Used by anodes: %s.' % names_tex(anode_tns, 4)
            d.figure(fig, cap)

    # Simulation ------------------------------------------------------

    def sec_sim(self):
        d = self.doc
        jobs = dict((k, t) for k, _, t in detect_jobs(self.cfg))
        if 'sim' not in jobs:
            return
        d.section('Simulation')
        d.para('Components of the simulation chain.  Parameters are grouped over all instances of '
               'each component type; where values differ between instances each distinct value '
               'is listed with the instances that use it.')
        order = ['wclsSimDepoSource', 'wclsSimDepoSetSource', 'TrackDepos', 'BlipSource', 'DepoBagger',
                 'DepoSetDrifter', 'Drifter', 'DepoSetFilter', 'DepoTransform', 'Ductor', 'MultiDuctor',
                 'DepoSplat', 'DepoZipper', 'Reframer', 'AddNoise', 'EmpiricalNoiseModel', 'Digitizer',
                 'wclsDepoFluxWriter', 'wclsSimChannelSink', 'wclsFrameSaver']
        for t in order:
            comps = self.cfg.of_type(t)
            if not comps:
                continue
            d.subsection(ttex(t))
            d.group_table(comps)
        self.plot_noise()

    def plot_noise(self):
        d = self.doc
        models = self.cfg.of_type('EmpiricalNoiseModel')
        if not models:
            return
        files = collections.OrderedDict()
        for m in models:
            files.setdefault(m['data'].get('spectra_file'), []).append(m)
        d.subsection('Noise spectra')
        d.para(r'\texttt{EmpiricalNoiseModel} components provide, per plane, noise amplitude '
               r'spectra measured or modeled for several wire lengths; the model interpolates '
               r'in wire length and \texttt{AddNoise} adds random noise with this spectrum to '
               r'the (voltage) waveforms before digitization.  Amplitudes are shown as stored '
               r'(WCT voltage units converted to mV).')
        for fn, ms in files.items():
            p = self.paths.find_data(fn)
            if not p:
                d.para(r'\emph{Noise file %s not found.}' % ttex(fn))
                continue
            spectra = load_json(p)
            anodes = sorted(set(tn(self.anode_of(m)) for m in ms if self.anode_of(m)))
            planes = sorted(set(s['plane'] for s in spectra))
            fig, axes = plt.subplots(1, len(planes), figsize=(5 * len(planes), 4), squeeze=False)
            rows = []
            for ax, pl in zip(axes[0], planes):
                ss = sorted([s for s in spectra if s['plane'] == pl], key=lambda s: s.get('wirelen', 0))
                cmap = plt.get_cmap('viridis')
                for i, s in enumerate(ss):
                    f = np.array(s['freqs']) / U.MHz
                    a = np.array(s['amps']) / U.mV
                    ax.plot(f, a, color=cmap(i / max(1, len(ss) - 1)), lw=0.8,
                            label='%.0f cm' % (s.get('wirelen', 0) / U.cm))
                    rows.append((str(pl), '%.1f' % (s.get('wirelen', 0) / U.cm), num(s.get('gain', 0) / (U.mV / U.fC), 4),
                                 num(s.get('shaping', 0) / U.us, 4), str(s.get('nsamples')), num(s.get('period', 0) / U.us, 4),
                                 num(s.get('const', 0) / U.mV, 4), str(len(s['freqs']))))
                ax.set_xlabel('frequency [MHz]')
                ax.set_ylabel('amplitude [mV]')
                ax.set_title('plane %d' % pl)
                ax.legend(fontsize=6, title='wire length', title_fontsize=7)
                ax.grid(alpha=0.3)
            d.figure(fig, r'Noise spectra in %s, used by %s (%s).' % (
                ttex(fn), names_tex([tn(m) for m in ms], 3),
                'anodes ' + names_tex(anodes, 4) if anodes else 'no anode'))
            d.table(['Plane', 'Wire length [cm]', 'Gain [mV/fC]', r'Shaping [$\mu$s]', 'nsamples',
                     r'Period [$\mu$s]', 'Const [mV]', 'nfreqs'], rows,
                    ['0.06\\linewidth', '0.12\\linewidth', '0.12\\linewidth', '0.12\\linewidth',
                     '0.1\\linewidth', '0.12\\linewidth', '0.1\\linewidth', '0.08\\linewidth'],
                    'Entries of noise spectra file %s.' % ttex(fn), fontsize=r'\scriptsize')

    # Noise filtering -------------------------------------------------

    def sec_nf(self):
        d = self.doc
        nfs = self.cfg.of_type('OmnibusNoiseFilter')
        dbs = self.cfg.of_type('OmniChannelNoiseDB')
        if not nfs and not dbs:
            return
        d.section('Noise filtering')
        if nfs:
            d.subsection(ttex('OmnibusNoiseFilter'))
            d.group_table(nfs)
        filt_types = set()
        for n in nfs:
            for k in ('channel_filters', 'grouped_filters', 'channel_status_filters', 'multigroup_chanfilters'):
                v = n['data'].get(k)
                if isinstance(v, list):
                    for x in v:
                        if isinstance(x, str):
                            filt_types.add(x)
                        elif isinstance(x, dict):
                            filt_types.update(y for y in x.get('channelgroups', []) if isinstance(y, str))
        fcomps = [self.cfg.get(f) for f in sorted(filt_types) if self.cfg.get(f)]
        bytype = collections.OrderedDict()
        for c in fcomps:
            bytype.setdefault(c['type'], []).append(c)
        for t, cs in bytype.items():
            self.covered.add(t)
            d.subsection('Noise filter %s' % ttex(t))
            d.group_table(cs)
        if dbs:
            self.nf_chndb(dbs)

    def nf_chndb(self, dbs):
        '''Summarize OmniChannelNoiseDB components, merging identical ones.'''
        d = self.doc
        big = ('groups', 'bad', 'channel_info', 'top_u_groups', 'top_v_groups', 'top_w_groups')
        d.subsection('Channel noise databases')
        d.para(r'\texttt{OmniChannelNoiseDB} components provide per-channel noise filter '
               r'parameters: coherent noise groups, bad channels and a list of '
               r'\texttt{channel\_info} entries, each applying parameters to a set of channels '
               r'(later entries override earlier ones).')
        d.group_table(dbs, skip=tuple(r'^%s$' % k for k in big))
        rows = []
        for db in dbs:
            dd = db['data']
            groups = dd.get('groups', [])
            sizes = sorted(set(len(g) for g in groups if isinstance(g, list)))
            bad = dd.get('bad', [])
            rows.append((ttex(tn(db)), str(len(groups)), ', '.join(map(str, sizes[:8])),
                         str(len(bad)) if isinstance(bad, list) else '?', str(len(dd.get('channel_info', [])))))
        d.table(['Database', 'Coherent groups', 'Group sizes', 'Bad channels', 'channel_info entries'], rows,
                ['0.34\\linewidth', '0.13\\linewidth', '0.2\\linewidth', '0.12\\linewidth', '0.14\\linewidth'],
                'Channel groups and bad channels per database.')

        def norm_ci(ci):
            out = []
            for ent in ci:
                e = dict(ent)
                ch = e.pop('channels', None)
                e['#channels'] = len(ch) if isinstance(ch, list) else ch
                out.append(e)
            return out
        cigroups = collections.OrderedDict()
        for db in dbs:
            ci = norm_ci(db['data'].get('channel_info', []))
            cigroups.setdefault(hashable(ci), (ci, []))[1].append(tn(db))
        for ci, names in cigroups.values():
            if not ci:
                continue
            rows = []
            for i, ent in enumerate(ci):
                nch = ent.get('#channels')
                items = [(k, v) for k, v in ent.items() if k != '#channels']
                for j, (k, v) in enumerate(items):
                    rows.append((str(i) if j == 0 else '', ('all' if nch is None else str(nch)) if j == 0 else '',
                                 ttex(k), fmt_value(k, v, 'OmniChannelNoiseDB', maxlen=90)))
            d.table(['Entry', 'Channels', 'Parameter', 'Value'], rows,
                    ['0.06\\linewidth', '0.1\\linewidth', '0.2\\linewidth', '0.58\\linewidth'],
                    r'\texttt{channel\_info} entries of %s (channel lists replaced by their size).' % names_tex(names, 4),
                    fontsize=r'\scriptsize')

    # Signal processing ----------------------------------------------

    def sec_sp(self):
        d = self.doc
        sps = self.cfg.of_type('OmnibusSigProc')
        if not sps:
            return
        d.section('Signal processing')
        d.para(r'\texttt{OmnibusSigProc} deconvolves the field and electronics responses from '
               r'the waveforms (2D in time and wire), applying time-domain (Wiener, Gaussian), '
               r'low-frequency (ROI) and wire-domain filters, and finds regions of interest.')
        d.subsection(ttex('OmnibusSigProc'))
        d.group_table(sps)
        filt = self.cfg.of_type('HfFilter', 'LfFilter')
        if filt:
            d.subsection('Filters')
            rows = []
            for c in filt:
                dd = c['data']
                for k, v in dd.items():
                    rows.append((ttex(c.get('name', '')), ttex(k), fmt_value(k, v, c['type'], dd), ttex(c['type'])))
            d.table(['Filter', 'Parameter', 'Value', 'Type'], rows,
                    ['0.22\\linewidth', '0.14\\linewidth', '0.38\\linewidth', '0.14\\linewidth'],
                    'Filter components.  HfFilter: $\\exp(-\\frac{1}{2}(f/\\sigma)^{p})$ (zero at $f=0$ if flag); '
                    'LfFilter: $1-\\exp(-(f/\\tau)^2)$.', fontsize=r'\scriptsize')
            self.plot_filters(filt)
        self.plot_sp_spectra(sps)
        for t in ('L1SPFilter', 'L1SPFilterPD', 'DNNROIFinding'):
            comps = self.cfg.of_type(t)
            if comps:
                d.subsection(ttex(t))
                d.group_table(comps)

    @staticmethod
    def filter_curve(c, f):
        dd = c['data']
        if c['type'] == 'HfFilter':
            return hf_filter(f, dd.get('sigma', 0), dd.get('power', 2), dd.get('flag', True))
        return lf_filter(f, dd.get('tau', 0))

    def plot_filters(self, filt):
        d = self.doc
        timef = [c for c in filt if c['data'].get('max_freq', 1) != 1]
        wiref = [c for c in filt if c['data'].get('max_freq', 1) == 1]
        cats = collections.OrderedDict()
        for c in timef:
            n = c.get('name', '')
            key = 'Low frequency (ROI) filters' if c['type'] == 'LfFilter' else \
                  'Wiener tight' if 'Wiener_tight' in n else 'Wiener wide' if 'Wiener_wide' in n else \
                  'Gaussian' if 'Gaus' in n else 'Other high frequency'
            cats.setdefault(key, []).append(c)
        if cats:
            n = len(cats)
            ncol = min(n, 3)
            nrow = int(math.ceil(n / ncol))
            fig, axes = plt.subplots(nrow, ncol, figsize=(5 * ncol, 3.8 * nrow), squeeze=False)
            for ax, (key, cs) in zip(axes.flat, cats.items()):
                for c in cs:
                    mf = c['data'].get('max_freq', U.MHz)
                    f = np.linspace(0, mf, 1000)
                    if key.startswith('Low'):
                        f = np.logspace(-5, 0, 1000) * mf
                    ls = '-' if c['type'] == 'HfFilter' else '--'
                    ax.plot(f / U.MHz, self.filter_curve(c, f), ls, label=c.get('name', ''), lw=1)
                ax.set_title(key, fontsize=9)
                ax.set_xlabel('frequency [MHz]')
                ax.set_ylabel('filter')
                ax.legend(fontsize=6)
                ax.grid(alpha=0.3)
                if key.startswith('Low'):
                    ax.set_xscale('log')
                    ax.set_xlim(1e-5 * mf / U.MHz, mf / U.MHz)
            for ax in list(axes.flat)[n:]:
                ax.set_visible(False)
            d.figure(fig, 'Time-domain filters vs frequency.  Source: HfFilter/LfFilter components named in legends.')
        if wiref:
            fig, ax = plt.subplots(figsize=(7, 3.5))
            for c in wiref:
                f = np.linspace(0, 1, 500)
                ax.plot(f, self.filter_curve(c, f), label=c.get('name', ''))
            ax.set_xlabel('wire frequency (fraction of Nyquist)')
            ax.set_ylabel('filter')
            ax.legend(fontsize=7)
            ax.grid(alpha=0.3)
            d.figure(fig, 'Wire-domain filters.  Source: HfFilter components named in legend.')

    def plot_sp_spectra(self, sps):
        '''Overlay response spectrum with the filters used, per SP config.'''
        d = self.doc
        groups = collections.OrderedDict()
        for c in sps:
            dd = c['data']
            key = (dd.get('field_response'), dd.get('elecresponse'),
                   tuple(dd.get('Wiener_tight_filters', ['Wiener_tight_U', 'Wiener_tight_V', 'Wiener_tight_W'])),
                   tuple(dd.get('Wiener_wide_filters', ['Wiener_wide_U', 'Wiener_wide_V', 'Wiener_wide_W'])),
                   dd.get('Gaus_wide_filter', 'Gaus_wide'))
            groups.setdefault(key, []).append(c)
        for (frn, ern, wt, ww, gw), cs in groups.items():
            fr = self.field(self.cfg.get(frn)) if frn else None
            erc = self.cfg.get(ern) if ern else None
            if fr is None or erc is None:
                continue
            er = ElecResponse(erc, self.paths)
            if not er.ok:
                continue
            tick = 0.5 * U.us
            for t in self.cfg.of_type('Reframer', 'DepoTransform'):
                tick = t['data'].get('tick', tick)
                break
            fig, axes = plt.subplots(1, len(fr.planes), figsize=(5 * len(fr.planes), 4), squeeze=False)
            for i, ax in enumerate(axes[0]):
                t, v = total_response(fr, i, er)
                if t is None:
                    continue
                step = max(1, int(round(tick / fr.period)))
                vs = v[::step]
                n = 4096
                spec = np.abs(np.fft.rfft(vs, n))
                f = np.fft.rfftfreq(n, tick)
                ax.plot(f / U.MHz, spec / spec.max(), 'k-', lw=1.5, label='|response| (norm.)')
                pid = fr.planes[i]['planeid']
                for names, ls in ((wt, '-'), (ww, '--'), ((gw,), ':')):
                    nm = names[pid] if len(names) > pid else names[-1]
                    fc = self.cfg.get('HfFilter:' + nm)
                    if fc:
                        ax.plot(f / U.MHz, self.filter_curve(fc, f), ls, label=nm)
                ax.set_xlabel('frequency [MHz]')
                ax.set_title('plane %d (%s)' % (pid, PLANE_NAMES.get(pid, '?')))
                ax.legend(fontsize=7)
                ax.grid(alpha=0.3)
            d.figure(fig, r'Magnitude of the response spectrum used in deconvolution (field response %s wire-0 '
                     r'average $\otimes$ %s, sampled at %.2f $\mu$s) with the filters applied by %s.'
                     % (ttex(frn), ttex(ern), tick / U.us, names_tex([tn(c) for c in cs], 3)))

    # Other -------------------------------------------------------------

    COVERED = {'AnodePlane', 'WireSchemaFile', 'FieldResponse', 'ColdElecResponse', 'JsonElecResponse',
               'RCResponse', 'PlaneImpactResponse', 'Drifter', 'DepoSetDrifter', 'DepoTransform', 'Ductor',
               'DepoSplat', 'DepoZipper', 'Reframer', 'AddNoise', 'EmpiricalNoiseModel', 'Digitizer',
               'OmnibusSigProc', 'HfFilter', 'LfFilter', 'OmnibusNoiseFilter', 'OmniChannelNoiseDB',
               'wclsSimDepoSource', 'wclsSimDepoSetSource', 'wclsDepoFluxWriter', 'wclsSimChannelSink',
               'wclsFrameSaver', 'TrackDepos', 'BlipSource', 'DepoBagger', 'DepoSetFilter', 'MultiDuctor',
               'L1SPFilter', 'L1SPFilterPD', 'DNNROIFinding'}

    def sec_other(self):
        d = self.doc
        types = self.cfg.types()
        physics_like = [t for t in sorted(types) if t not in self.covered and not self.is_plumbing(t)]
        if not physics_like:
            return
        d.section('Other components')
        d.para('Components not covered above which may affect results.')
        for t in physics_like:
            d.subsection(ttex(t))
            d.group_table(self.cfg.of_type(t))

    @staticmethod
    def is_plumbing(t):
        return bool(re.match(r'^(Frame|Depo)(Fanin|Fanout|Merger|Splitter|Sync|SetFanout|SetFanin)$|^(Retagger|DumpFrames|'
                             r'DeposOrBust|FrameSync|Random|FftwDFT|TorchDFT|Pgrapher|TbbFlow|MegaAnodePlane|'
                             r'DumpDepos|FrameFileSink|FrameFileSource|DepoFileSink|DepoFileSource|'
                             r'HDF5FrameTap|MagnifySink|ExampleROOTAna|ChannelSelector|ChannelSplitter|'
                             r'FrameFanin|FrameFanout|DepoSetFanout|DepoSetFanin|TensorFileSink)$', t))

    def sec_appendix(self):
        if not self.args.appendix:
            return
        d = self.doc
        d.section('Appendix: all components')
        for t in sorted(self.cfg.types()):
            d.subsection(ttex(t))
            d.group_table(self.cfg.of_type(t))

    def sec_warnings(self):
        if not WARNINGS:
            return
        d = self.doc
        d.section('Warnings')
        d.para(r'\begin{itemize}' + ''.join(r'\item ' + tex(w) for w in WARNINGS) + r'\end{itemize}')


# ----------------------------------------------------------------------
# Main
# ----------------------------------------------------------------------

def parse_args(argv):
    p = argparse.ArgumentParser(description=__doc__.split('\n\n')[0],
                                formatter_class=argparse.RawDescriptionHelpFormatter,
                                epilog=__doc__.split('\n\n', 1)[1])
    p.add_argument('fcl', help='FHiCL file')
    p.add_argument('jsonnet', nargs='?', default=None,
                   help='WCT Jsonnet main file (name as in wcls_main.configs or a path). '
                   'Default: the first wcls_main.configs entry of the selected table')
    p.add_argument('-o', '--output', default='wcls-summary.pdf', help='output PDF file')
    p.add_argument('--block', default=None, help='dotted FHiCL path of the table holding wcls_main')
    p.add_argument('--list', action='store_true', help='list WCLS tables found in the FHiCL and exit')
    p.add_argument('-I', '--fcl-path', action='append', default=[], help='FHiCL include directory (repeatable)')
    p.add_argument('--fcl-tree', action='append', default=[],
                   help='directory tree searched recursively for FHiCL includes by base name (repeatable)')
    p.add_argument('-J', '--jpath', action='append', default=[],
                   help='Jsonnet/WCT search path, searched before $WIRECELL_PATH (repeatable)')
    p.add_argument('-P', '--data-path', action='append', default=[],
                   help='additional directory for WCT data files (repeatable)')
    p.add_argument('--data-tree', action='append', default=[],
                   help='directory tree searched recursively for WCT data files by base name (repeatable)')
    p.add_argument('--map', action='append', default=[],
                   help='make a Jsonnet import prefix point to a directory, eg pgrapher/experiment=/src/DUNEWireCell')
    p.add_argument('--auto', dest='auto', action='store_true', default=True,
                   help='auto-locate FHiCL/Jsonnet/data search paths relative '
                        'to this script + WIRECELL_PATH so the in-tree '
                        'dune10kt-vd example runs without a full DUNE '
                        'environment (default: on)')
    p.add_argument('--no-auto', dest='auto', action='store_false',
                   help='disable --auto path discovery')
    p.add_argument('--dunewirecell', default=None,
                   help='a dunereco DUNEWireCell source directory: maps pgrapher/{experiment,common,dune,ui}, '
                   'searches it for data files and FHiCL includes')
    p.add_argument('--ext-str', action='append', default=[], help='override/add an external string var key=value')
    p.add_argument('--ext-code', action='append', default=[], help='override/add an external code var key=value')
    p.add_argument('--jsonnet-exe', default=os.environ.get('JSONNET', 'jsonnet'), help='jsonnet program')
    p.add_argument('--save-json', default=None, help='also save the compiled WCT configuration to this file')
    p.add_argument('--keep', default=None, help='keep the LaTeX source and figures in this directory')
    p.add_argument('--max-anode-rows', type=int, default=400, help='maximum rows in the per-face table')
    p.add_argument('--no-appendix', dest='appendix', action='store_false',
                   help='do not include the appendix listing all components')
    return p.parse_args(argv)


# ----------------------------------------------------------------------
# Auto-location of search paths
# ----------------------------------------------------------------------

def _ancestors(start, levels):
    '''`start` and up to `levels` parent directories, nearest first.'''
    d = os.path.abspath(start)
    out = [d]
    for _ in range(levels):
        nd = os.path.dirname(d)
        if nd == d:
            break
        d = nd
        out.append(d)
    return out


def _first_glob(ancestors, patterns):
    '''First existing path matching any <ancestor>/<pattern> (bounded globs,
    no recursive **), scanning ancestors nearest-first.'''
    for anc in ancestors:
        for pat in patterns:
            hits = sorted(glob.glob(os.path.join(anc, pat)))
            if hits:
                return hits[0]
    return None


def _wcpath_entries(args_paths):
    '''Explicit paths plus WIRECELL_PATH entries.'''
    return [p for p in (list(args_paths) +
                        os.environ.get('WIRECELL_PATH', '').split(':')) if p]


def autolocate(args):
    '''Best-effort discovery of the Jsonnet / FHiCL / data search paths needed
    to summarize the in-tree dune10kt-vd example WITHOUT a full DUNE shell
    environment.

    Everything discovered here is APPENDED to whatever the user passed
    explicitly (-I/--fcl-tree/-J/-P/--dunewirecell), so explicit flags always
    keep priority and nothing is overridden.  Discovery is relative to this
    script's location (a dunereco checkout) plus WIRECELL_PATH, so it works
    from a bare dunereco + wire-cell-toolkit checkout.  Disable with --no-auto.

    What it tries to find, and why each is needed:
      * --dunewirecell : the DUNEWireCell dir holding this script -- it carries
          dune10kt-vd/, common/, dune/ ... i.e. the pgrapher/experiment tree the
          config Jsonnet imports.
      * --fcl-tree <dunereco> : so the WC FHiCL (#include "wirecell_dune.fcl",
          the per-detector fcls) resolve by basename.
      * -J <wire-cell-toolkit/cfg> : wirecell.jsonnet / pgraph.jsonnet /
          pgrapher/common, imported by the config Jsonnet (unless already on
          WIRECELL_PATH).
      * -P <wire-cell-data> : field-response / wire-geometry .json.bz2 data
          (unless already on WIRECELL_PATH).
      * --fcl-tree <DUNE fcl install> : the dunecore detector-service fcls
          (detectorproperties_dune.fcl, detectorclocks_dune.fcl, ...) that the
          base fcl's structs reference via @local::.  These are NOT in dunereco;
          if a light install tree is nearby we index it so most @local:: resolve.
          A few numbers that live only in the heavy larsoft sim-service chain
          (e.g. dunefd_largeantparameters' DL/DT) cannot be found this way --
          supply them with --ext-code if they end up unresolved.
    '''
    if not args.auto:
        return
    here = os.path.dirname(os.path.abspath(__file__))
    dwc = os.path.dirname(here)                 # scripts/ -> DUNEWireCell/
    ancestors = _ancestors(here, 8)
    found = []

    # 1. --dunewirecell (the DUNEWireCell dir with dune10kt-vd/ etc.)
    if not args.dunewirecell and os.path.isdir(os.path.join(dwc, 'dune10kt-vd')):
        args.dunewirecell = dwc
        found.append('--dunewirecell %s' % dwc)

    # 2. --fcl-tree for the dunereco checkout (WC fcls, resolved by basename).
    reco = next((a for a in ancestors
                 if os.path.basename(a) == 'dunereco'
                 and os.path.isdir(os.path.join(a, 'dunereco'))), None)
    if reco is None:
        reco = os.path.dirname(os.path.dirname(dwc))  # .../DUNEWireCell -> repo
    if reco and os.path.isdir(reco) and reco not in args.fcl_tree:
        args.fcl_tree.append(reco)
        found.append('--fcl-tree %s' % reco)

    # 3. -J for the wire-cell-toolkit cfg library (if not already available).
    have_wct = any(os.path.exists(os.path.join(p, 'wirecell.jsonnet'))
                   for p in _wcpath_entries(args.jpath))
    if not have_wct:
        cfg = _first_glob(ancestors, [
            'wire-cell-toolkit/cfg/wirecell.jsonnet',
            '*/wire-cell-toolkit/cfg/wirecell.jsonnet',
            '*/*/wire-cell-toolkit/cfg/wirecell.jsonnet',
            '*/cfg/wirecell.jsonnet',
        ])
        if cfg:
            d = os.path.dirname(cfg)
            args.jpath.append(d)
            found.append('-J %s' % d)

    # 4. -P for wire-cell-data (if no response data already available).
    have_data = any(glob.glob(os.path.join(p, '*.json.bz2'))
                    for p in _wcpath_entries(args.data_path))
    if not have_data:
        wcd = _first_glob(ancestors, [
            'wire-cell-data', '*/wire-cell-data', '*/*/wire-cell-data',
        ])
        if wcd and os.path.isdir(wcd):
            args.data_path.append(wcd)
            found.append('-P %s' % wcd)

    # 5. --fcl-tree for a nearby DUNE fcl install (dunecore service fcls).
    dp = _first_glob(ancestors, [
        'fcl/detectorproperties_dune.fcl',
        '*/fcl/detectorproperties_dune.fcl',
        '*/*/fcl/detectorproperties_dune.fcl',
        '*/*/*/fcl/detectorproperties_dune.fcl',
        '*/*/*/*/fcl/detectorproperties_dune.fcl',
    ])
    if dp:
        d = os.path.dirname(dp)
        if d not in args.fcl_tree:
            args.fcl_tree.append(d)
            found.append('--fcl-tree %s' % d)

    if found:
        info('auto-located (use --no-auto to disable):\n    ' +
             '\n    '.join(found))
    else:
        info('auto-location found nothing to add (use --no-auto to disable)')


def main(argv=None):
    args = parse_args(argv if argv is not None else sys.argv[1:])
    tmpdir = tempfile.mkdtemp(prefix='wcls-summary-')
    try:
        return run(args, tmpdir)
    finally:
        if not args.keep:
            shutil.rmtree(tmpdir, ignore_errors=True)


def run(args, tmpdir):
    autolocate(args)
    paths = Paths(args, tmpdir)
    incdirs = list(args.fcl_path)
    if args.dunewirecell:
        incdirs.append(os.path.abspath(args.dunewirecell))
    fdoc = FhiclDocument(incdirs, args.fcl_tree)
    fdoc.load(args.fcl)
    blocks = find_wcls_blocks(fdoc)
    if args.list:
        for path, tab in blocks:
            print('%-60s %s' % (path, ' '.join(fhicl_strings(tab['wcls_main'].get('configs')))))
        return 0
    if not blocks:
        raise SystemExit('no table with wcls_main.configs found in %s' % args.fcl)
    alternatives = []
    if args.block:
        sel = [b for b in blocks if b[0] == args.block]
        if not sel:
            raise SystemExit('no WCLS table %s, use --list' % args.block)
    elif args.jsonnet:
        want = args.jsonnet
        sel = [b for b in blocks
               if any(c == want or c.endswith('/' + want) or want.endswith('/' + c) or
                      os.path.basename(c) == os.path.basename(want)
                      for c in fhicl_strings(b[1]['wcls_main'].get('configs')))]
        if not sel:
            info('no WCLS table uses %s, using the first WCLS table' % want)
            sel = blocks
    else:
        sel = blocks
    if len(sel) > 1:
        alternatives = [b[0] for b in sel[1:]]
        info('using %s (also matching: %s)' % (sel[0][0], ', '.join(alternatives)))
    block, btab = sel[0]
    wm = btab['wcls_main']
    extvars, unresolved = make_extvars(wm, args.ext_str, args.ext_code)
    if unresolved:
        msg = '\n'.join('  %s = %s (at %s)' % (k, u.ref, u.src) for k, u in unresolved)
        hint = ('  --ext-code %s=<value>' % k for k, _ in unresolved)
        raise SystemExit(
            'unresolved FHiCL values for external variables:\n%s\n'
            'These @local:: references resolve only against FHiCL that could not\n'
            'be loaded here.  Values from detectorproperties/detectorclocks live\n'
            'in a DUNE fcl install (auto-located if one is near this checkout);\n'
            'a few (e.g. dunefd_largeantparameters DL/DT) live only in the heavy\n'
            'larsoft sim-service chain, which is NOT in dunereco.  Either add the\n'
            'right FHiCL path (-I, --fcl-tree) or supply the values directly, e.g.\n%s'
            % (msg, '\n'.join(hint)))
    fpaths = [os.path.abspath(p) for p in fhicl_strings(wm.get('paths'))]
    configs = fhicl_strings(wm.get('configs'))
    if args.jsonnet:
        jfiles = [args.jsonnet]
    else:
        jfiles = configs
    comps = []
    commands = []
    found_files = []
    for jf in jfiles:
        p = paths.find_jsonnet(jf, fpaths)
        if p is None and args.jsonnet and configs:
            # the named file may be a short form of the configs entry
            for c in configs:
                if os.path.basename(c) == os.path.basename(jf):
                    p = paths.find_jsonnet(c, fpaths)
        if p is None:
            raise SystemExit('Jsonnet file %s not found in search path: %s'
                             % (jf, ':'.join(paths.all_jpaths(fpaths))))
        found_files.append(p)
        c, cmd = compile_jsonnet(args.jsonnet_exe, p, paths.all_jpaths(fpaths), extvars)
        comps += c
        commands.append(cmd)
    if args.save_json:
        with open(args.save_json, 'w') as fp:
            json.dump(comps, fp, indent=1)
    cfg = WCConfig(comps)
    info('%d components, job types: %s' % (len(cfg.comps), ', '.join(t for _, t, _ in detect_jobs(cfg))))

    figdir = args.keep or os.path.join(tmpdir, 'tex')
    os.makedirs(figdir, exist_ok=True)
    doc = Doc(figdir)
    jobinfo = dict(fcl=os.path.abspath(args.fcl), block=block, block_table=btab, alternatives=alternatives,
                   wcls_main=wm, extvars=extvars, jsonnet_files=found_files,
                   command=[_shellquote(x) for x in commands[0]] if commands else [],
                   missing_includes=fdoc.missing)
    summ = Summary(args, paths, cfg, doc, figdir)
    summ.run(jobinfo)
    jobs = ', '.join(t for _, t, _ in detect_jobs(cfg)) or 'unknown'
    title = r'Wire-Cell configuration summary\\[0.5ex]\large %s' % ttex(block)
    subtitle = r'%s\\ %s\\ job types: %s' % (ttex(os.path.basename(args.fcl)),
                                             ', '.join(ttex(os.path.basename(f)) for f in found_files), tex(jobs))
    texsrc = doc.render(title, subtitle)
    texfile = os.path.join(figdir, 'summary.tex')
    with open(texfile, 'w') as fp:
        fp.write(texsrc)
    for i in range(2):
        r = subprocess.run(['pdflatex', '-interaction=nonstopmode', '-halt-on-error', 'summary.tex'],
                           cwd=figdir, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
        if r.returncode != 0:
            sys.stderr.write('\n'.join(r.stdout.splitlines()[-40:]) + '\n')
            if not args.keep:
                keep = tempfile.mkdtemp(prefix='wcls-summary-failed-')
                shutil.copytree(figdir, keep, dirs_exist_ok=True)
                sys.stderr.write('LaTeX sources kept in %s\n' % keep)
            raise SystemExit('pdflatex failed')
    shutil.copy(os.path.join(figdir, 'summary.pdf'), args.output)
    info('wrote %s' % args.output)
    return 0


def _shellquote(s):
    if re.match(r'^[A-Za-z0-9_./=:+-]+$', s):
        return s
    return "'" + s.replace("'", "'\\''") + "'"


if __name__ == '__main__':
    sys.exit(main())
