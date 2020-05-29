"""
I/O functions for working with CmdStan executables.

"""

import os
import re
import threading
import numpy as np
import tempfile


def _rdump_array(key, val):
    c = 'c(' + ', '.join(map(str, val.T.flat)) + ')'
    if (val.size, ) == val.shape:
        return '{key} <- {c}'.format(key=key, c=c)
    else:
        dim = '.Dim = c{0}'.format(val.shape)
        struct = '{key} <- structure({c}, {dim})'.format(key=key, c=c, dim=dim)
        return struct


def rdump(fname, data):
    """Dump a dict of data to a R dump format file.
    """
    with open(fname, 'w') as fd:
        for key, val in data.items():
            if isinstance(val, np.ndarray) and val.size > 1:
                line = _rdump_array(key, val)
            else:
                try:
                    val = val.flat[0]
                except:
                    pass
                line = '%s <- %s' % (key, val)
            fd.write(line)
            fd.write('\n')


def rload(fname):
    """Load a dict of data from an R dump format file.
    """
    with open(fname, 'r') as fd:
        lines = fd.readlines()
    data = {}
    for line in lines:
        lhs, rhs = [_.strip() for _ in line.split('<-')]
        if rhs.startswith('structure'):
            *_, vals, dim = rhs.replace('(', ' ').replace(')', ' ').split('c')
            vals = [float(v) for v in vals.split(',')[:-1]]
            dim = [int(v) for v in dim.split(',')]
            val = np.array(vals).reshape(dim[::-1]).T
        elif rhs.startswith('c'):
            val = np.array([float(_) for _ in rhs[2:-1].split(',')])
        else:
            try:
                val = int(rhs)
            except:
                try:
                    val = float(rhs)
                except:
                    raise ValueError(rhs)
        data[lhs] = val
    return data


def merge_csv_data(*csvs, skip=0):
    """Merge multiple CSV dicts into a single dict.
    """
    data_ = {}
    for csv in csvs:
        for key, val in csv.items():
            # XXX do better
            if key in 'loo loos ks'.split():
                continue
            val = val[skip:]
            if key in data_:
                data_[key] = np.concatenate((data_[key], val), axis=0)
            else:
                data_[key] = val
    return data_


def parse_csv(fname, merge=True):
    """Parse samples from a Stan output CSV file.
    """
    if '*' in fname:
        import glob
        return parse_csv(glob.glob(fname), merge=merge)
    if isinstance(fname, (list, tuple)):
        csv = []
        for _ in fname:
            try:
                csv.append(parse_csv(_))
            except Exception as e:
                print('skipping ', fname, e)
        if merge:
            csv = merge_csv_data(*csv)
        return csv

    lines = []
    with open(fname, 'r') as fd:
        for line in fd.readlines():
            if not line.startswith('#'):
                lines.append(line.strip().split(','))
    names = [field.split('.') for field in lines[0]]
    data = np.array([[float(f) for f in line] for line in lines[1:]])

    namemap = {}
    maxdims = {}
    for i, name in enumerate(names):
        if name[0] not in namemap:
            namemap[name[0]] = []
        namemap[name[0]].append(i)
        if len(name) > 1:
            maxdims[name[0]] = name[1:]

    for name in maxdims.keys():
        dims = []
        for dim in maxdims[name]:
            dims.append(int(dim))
        maxdims[name] = tuple(reversed(dims))

    # data in linear order per Stan, e.g. mat is col maj
    # TODO array is row maj, how to distinguish matrix v array[,]?
    data_ = {}
    for name, idx in namemap.items():
        new_shape = (-1, ) + maxdims.get(name, ())
        data_[name] = data[:, idx].reshape(new_shape)

    return data_


def parse_summary_csv(fname):
    """Parse CSV output of the stansummary program.
    """
    skeys = []
    svals = []
    niter = -1
    with open(fname, 'r') as fd:
        scols = fd.readline().strip().split(',')
        for line in fd.readlines():
            if 'iterations' in line:
                niter_match = re.search(r'(\d+) iterations saved', line)
                if niter_match:
                    niter = int(niter_match.group(1))
                continue
            if line.startswith('#') or '"' not in line:
                continue
            _, k, v = line.split('"')
            skeys.append(k)
            svals.append(np.array([float(_) for _ in v.split(',')[1:]]))
    svals = np.array(svals)

    sdat = {}
    sdims = {}
    for skey, sval in zip(skeys, svals):
        if '[' in skey:
            name, dim = skey.replace('[', ']').split(']')[:-1]
            dim = tuple(int(i) for i in dim.split(','))
            sdims[name] = dim
            if name not in sdat:
                sdat[name] = []
            sdat[name].append(sval)
        else:
            sdat[skey] = sval

    for key in [_ for _ in sdat.keys()]:
        if key in sdims:
            sdat[key] = np.array(sdat[key]).reshape(sdims[key] + (-1, ))

    recs = {}
    dt = [(k, 'f8') for k in scols[1:]]
    for key, val in sdat.items():
        recs[key] = np.rec.array(val, dtype=dt)

    return niter, recs


def compile_model(cs, name):
    if os.path.exists(f'{name}.hpp') and os.path.getmtime(f'{name}.hpp') > os.path.getmtime(f'{name}.stan'):
        return
    cmd = f"set -eux; {cs}/bin/stanc --name=model --include_paths=$PWD --o=$PWD/{name}.hpp $PWD/{name}.stan"
    bash(cmd)
    #assert os.system(cmd)==0
    bash(f"make -C {cs} $PWD/{name}")
    

def diagnose_csvs(cs, *csvs):
    spcsvs = ' '.join(csvs)
    bash(f"{cs}/bin/diagnose {spcsvs}")
    

def run(cs, name, data=None, sampler_args=''):
    with tempfile.TemporaryDirectory() as td:
        rdump(f'{td}/data', data or {})
        bash(f'./{name} sample {sampler_args} '
                         f'data file={td}/data output file={td}/samples')
        diagnose_csvs(cs, f'{td}/samples')
        csv = parse_csv(f'{td}/samples')
    return csv


# this stuff makes bash stdout/stderr show up in jupyter notebook
# https://stackoverflow.com/a/59339154
import signal
import subprocess as sp
import sys


class VerboseCalledProcessError(sp.CalledProcessError):
    def __str__(self):
        if self.returncode and self.returncode < 0:
            try:
                msg = "Command '%s' died with %r." % (
                    self.cmd, signal.Signals(-self.returncode))
            except ValueError:
                msg = "Command '%s' died with unknown signal %d." % (
                    self.cmd, -self.returncode)
        else:
            msg = "Command '%s' returned non-zero exit status %d." % (
                self.cmd, self.returncode)

        return f'{msg}\n' \
               f'Stdout:\n' \
               f'{self.output}\n' \
               f'Stderr:\n' \
               f'{self.stderr}'


def bash(cmd, print_stdout=True, print_stderr=True):
    proc = sp.Popen(cmd, stderr=sp.PIPE, stdout=sp.PIPE, shell=True, universal_newlines=True,
                    executable='/bin/bash')

    all_stdout = []
    all_stderr = []
    while proc.poll() is None:
        for stdout_line in proc.stdout:
            if stdout_line != '':
                if print_stdout:
                    print(stdout_line, end='')
                all_stdout.append(stdout_line)
        for stderr_line in proc.stderr:
            if stderr_line != '':
                if print_stderr:
                    print(stderr_line, end='', file=sys.stderr)
                all_stderr.append(stderr_line)

    stdout_text = ''.join(all_stdout)
    stderr_text = ''.join(all_stderr)
    if proc.wait() != 0:
        raise VerboseCalledProcessError(proc.returncode, cmd, stdout_text, stderr_text)