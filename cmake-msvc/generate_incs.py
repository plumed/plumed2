# generate_incs.py - PLUMED build-time generators, replicated from upstream
# bash/awk scripts (maketools/makecmd, makeModuleMap, xxd, update-version.sh,
# and the Config.inc / Makefile.conf.xxd rules in src/config/Makefile).
#
# Usage (from CMake):
#   python generate_incs.py --src <plumed-src-root> --outdir <gen-dir>
#                           [--completion-only | --config-only]
#
# Outputs into --outdir:
#   PlumedMainEnum.inc / PlumedMainMap.inc
#   GREXEnum.inc       / GREXMap.inc
#   CLToolMainEnum.inc / CLToolMainMap.inc
#   ModuleMap.inc
#   completion.xxd
#   Config.inc / version.h / Makefile.conf.xxd
import re, os, glob, argparse, subprocess

def makecmd(cpp_path, out_dir, base):
    # upstream awk: fields[0]=="case" and fields[1] matches ^cmd_.*:$
    enum, string = [], []
    with open(cpp_path, encoding='utf-8', errors='replace') as fh:
        for line in fh:
            parts = line.split()
            if len(parts) >= 2 and parts[0] == 'case' and re.match(r'^cmd_[A-Za-z0-9_]+:$', parts[1]):
                name = parts[1].rstrip(':')
                enum.append(name)
                string.append(name[4:])
    with open(os.path.join(out_dir, base + 'Enum.inc'), 'w', newline='\n') as fh:
        fh.write(',\n'.join(enum) + ('\n' if enum else ''))
    with open(os.path.join(out_dir, base + 'Map.inc'), 'w', newline='\n') as fh:
        fh.write(',\n'.join('  {"%s",%s}' % (s, e) for s, e in zip(string, enum)) + ('\n' if enum else ''))
    return len(enum)

def make_module_map(src_root, out_path):
    rows = []
    # src_root may be the plumed root or the src/ dir; module.type lives in src/<mod>/
    base = src_root if os.path.isdir(os.path.join(src_root, 'src', 'core')) \
        else os.path.dirname(src_root)
    src_dir = os.path.join(base, 'src')
    pattern = os.path.join(src_dir, '*', 'module.type').replace(os.sep, '/')
    for mtf in glob.glob(pattern):
        moddir = os.path.dirname(mtf)
        for f in glob.glob(os.path.join(moddir, '*.cpp')):
            text = open(f, encoding='utf-8', errors='replace').read()
            for m in re.finditer(r'PLUMED_REGISTER_ACTION\s*\(\s*[0-9a-zA-Z_]+\s*,\s*"([0-9a-zA-Z_]+)"\s*\)', text):
                rows.append((m.group(1), os.path.basename(moddir)))
    with open(out_path, 'w', newline='\n') as fh:
        fh.write(',\n'.join('  {"%s", "%s"}' % (a, b) for a, b in rows) + ('\n' if rows else ''))
    return len(rows)

def to_xxd_bytes(data):
    out, j = [], 0
    for b in data:
        if j != 0 and j % 12 == 0:
            out.append(',\n ')
        elif j != 0:
            out.append(',')
        out.append(' 0x%02x' % b)
        j += 1
    out.append('')
    return ''.join(out)

def make_completion_xxd(src_root, out_dir):
    # upstream: { make-scripts-options.sh ; cat completion.sh ; } | xxd
    # make-scripts-options.sh runs each scripts/*.sh with --options. For the
    # Windows kernel these bash strings are inert data (Completion cltool only
    # prints them), so we run the scripts if bash is available, else emit a
    # minimal placeholder with the standard global options.
    lines = []
    scripts_dir = os.path.join(src_root, 'scripts')
    have_bash = False
    try:
        have_bash = os.system('bash -c "true"') == 0
    except Exception:
        pass
    if have_bash:
        for f in sorted(glob.glob(os.path.join(scripts_dir, '*.sh'))):
            name = os.path.basename(f)[:-3].replace('-', '_')
            r = os.popen('cd "%s" && PLUMED_ROOT="%s" bash scripts/%s --options 2>/dev/null' % (src_root, src_root, name)).read().strip()
            if r:
                lines.append('local cmd_keys_%s="%s"' % (name, r))
    else:
        lines.append('local cmd_keys_plumed="--help --options --version"')
    comp = open(os.path.join(src_root, 'src', 'cltools', 'completion.sh'),
                encoding='utf-8', errors='replace').read()
    data = ('\n'.join(lines) + '\n' + comp).encode('utf-8')
    with open(os.path.join(out_dir, 'completion.xxd'), 'w', newline='\n') as fh:
        fh.write(to_xxd_bytes(data))
    return len(data)

def git_describe(src_root):
    try:
        r = subprocess.run(['git', 'describe', '--long', '--dirty', '--always'],
                           cwd=src_root, capture_output=True, text=True, timeout=30)
        if r.returncode == 0 and r.stdout.strip():
            return r.stdout.strip()
    except Exception:
        pass
    try:
        r = subprocess.run(['git', 'rev-parse', '--short=12', 'HEAD'],
                           cwd=src_root, capture_output=True, text=True, timeout=30)
        if r.returncode == 0 and r.stdout.strip():
            return r.stdout.strip()
    except Exception:
        pass
    return 'Unknown'

def read_version(src_root):
    vf = os.path.join(src_root, 'VERSION.txt')
    if not os.path.isfile(vf):
        return 'Unknown', 'Unknown', 0, 0, 0
    line = ''
    with open(vf, encoding='utf-8', errors='replace') as fh:
        for l in fh:
            l = l.strip()
            if l and not l.startswith('#'):
                line = l
                break
    if not line:
        return 'Unknown', 'Unknown', 0, 0, 0
    m = re.match(r'^([0-9]+)\.([0-9]+)(?:\.([0-9]+))?', line)
    if not m:
        return 'Unknown', line, 0, 0, 0
    major = int(m.group(1))
    minor = int(m.group(2))
    patch = int(m.group(3)) if m.group(3) else 0
    short = '%d.%d' % (major, minor)
    return short, line, major, minor, patch

def make_config(src_root, out_dir):
    # --- Config.inc, following the sed rules in src/config/Makefile for the
    # non-install build (ConfigInstall.inc is not used by the static kernel).
    root = os.path.abspath(src_root).replace('\\', '/')
    text = open(os.path.join(src_root, 'src', 'config', 'Config.inc.in'),
                encoding='utf-8', errors='replace').read()
    text = text.replace('@SOEXT@', 'dll')
    text = text.replace('@ISINSTALLED@', 'false')
    text = text.replace('@htmldir@', 'xxxxNAxxxx')
    text = text.replace('@includedir@', 'xxxxNAxxxx')
    text = text.replace('@program_name@', 'plumed')
    text = text.replace('@PLUMED_ROOT@', root)
    with open(os.path.join(out_dir, 'Config.inc'), 'w', newline='\n') as fh:
        fh.write(text)

    # --- version.h, following src/maketools/update-version.sh
    short, long_, major, minor, patch = read_version(src_root)
    gitv = git_describe(src_root)
    vh = (
        '#ifndef __PLUMED_config_version_h\n'
        '#define __PLUMED_config_version_h\n'
        '#define PLUMED_VERSION_SHORT "%s"\n'
        '#define PLUMED_VERSION_LONG "%s"\n'
        '#define PLUMED_VERSION_GIT "%s"\n'
        '#define PLUMED_VERSION_MAJOR %d\n'
        '#define PLUMED_VERSION_MINOR %d\n'
        '#define PLUMED_VERSION_PATCH %d\n'
        '#endif\n'
    ) % (short, long_, gitv, major, minor, patch)
    with open(os.path.join(out_dir, 'version.h'), 'w', newline='\n') as fh:
        fh.write(vh)

    # --- Makefile.conf.xxd (embedded configuration report, used by getMakefile())
    conf = (
        '# PLUMED build configuration (generated by the CMake build)\n'
        'program_name=plumed\n'
        'version=%s\n'
        'version_git=%s\n'
    ) % (long_, gitv)
    with open(os.path.join(out_dir, 'Makefile.conf.xxd'), 'w', newline='\n') as fh:
        fh.write(to_xxd_bytes(conf.encode('utf-8')))
    return long_, gitv

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--src', required=True)
    ap.add_argument('--outdir', required=True)
    ap.add_argument('--completion-only', action='store_true')
    ap.add_argument('--config-only', action='store_true')
    args = ap.parse_args()
    src_root = os.path.abspath(args.src)
    out_dir = os.path.abspath(args.outdir)
    os.makedirs(out_dir, exist_ok=True)

    if args.config_only:
        ver, gitv = make_config(src_root, out_dir)
        print(f'  Config.inc / version.h / Makefile.conf.xxd (version {ver}, git {gitv})')
        return

    if not args.completion_only:
        core = os.path.join(src_root, 'src', 'core')
        for base, cpp in [('PlumedMain', 'PlumedMain.cpp'), ('GREX', 'GREX.cpp'),
                          ('CLToolMain', 'CLToolMain.cpp')]:
            n = makecmd(os.path.join(core, cpp), out_dir, base)
            print(f'  {base}: {n} commands')
        n = make_module_map(src_root, os.path.join(out_dir, 'ModuleMap.inc'))
        print(f'  ModuleMap: {n} actions')
        ver, gitv = make_config(src_root, out_dir)
        print(f'  Config.inc / version.h / Makefile.conf.xxd (version {ver}, git {gitv})')
    n = make_completion_xxd(src_root, out_dir)
    print(f'  completion.xxd: {n} bytes')

if __name__ == '__main__':
    main()
