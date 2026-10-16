#!/usr/bin/env python3
"""recover_nuclear.py -- bounded, stdlib-only, full-run MEGAHIT + BUSCO(genome) runner.
Prepares real BUSCO output for the existing validated downstream builder. It infers no
phylogeny, admits no markers, and makes no scientific acceptance claim (see NOTES.md)."""
import argparse, hashlib, json, os, subprocess, sys
from datetime import datetime, timezone
from pathlib import Path
import sys
if not __debug__:
    raise RuntimeError('Optimized Python disables required scientific assertions; run without -O.')
_HELPER_DIR = str(Path(__file__).resolve().parent)
if _HELPER_DIR not in sys.path:
    sys.path.insert(0, _HELPER_DIR)

INVARIANTS = [
    "raw reads are never modified; only hashed once and size/mtime_ns-snapshotted before/after",
    "receipt JSON files are hashed and syntax-checked only; schema and claims are not interpreted",
    "MEGAHIT k-list/min-count are program defaults; parameters are reported, not validated as optimal",
    "assembly statistics (N50 etc.) are descriptive only; no quality claim is made",
    "BUSCO completion does not imply marker admission, orthology, or biological truth",
    "the 2661 exclusion is applied downstream, before matrix selection, not by this script",
    "no phylogeny is inferred here; the downstream builder consumes work/busco",
    "no network, no tool installation, no shell=True, no runtime timeout; child nonzero = fail-stop",
]
FASTA_OK = set("ACGTNacgtn")

def utc(): return datetime.now(timezone.utc).isoformat()

def sha256(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""): h.update(chunk)
    return h.hexdigest()

def snap(p):
    s = os.stat(p); return {"path": str(p), "size": s.st_size, "mtime_ns": s.st_mtime_ns}

def atomic_json(path, obj):
    tmp = path.with_name(path.name + ".tmp")
    with open(tmp, "w", encoding="utf-8") as f:
        json.dump(obj, f, indent=2, sort_keys=True); f.write("\n")
    os.replace(tmp, path)

def req_file(p, what):
    p = Path(p)
    if not p.is_file(): raise ValueError(f"{what} is not an existing file: {p}")
    if p.stat().st_size == 0: raise ValueError(f"{what} is empty: {p}")
    return p.resolve()

def req_dir(p, what):
    p = Path(p)
    if not p.is_dir(): raise ValueError(f"{what} is not an existing directory: {p}")
    return p.resolve()

class Runner:
    def __init__(self, work, execute):
        self.work, self.execute, self.path = work, execute, work / "RUN_RECEIPT.json"
        self.receipt = {"tool": "recover_nuclear.py", "execute": execute, "stage": "init", "status": "RUNNING",
                        "utc_start": utc(), "utc_end": None, "commands": [], "hashes": {}, "error": None,
                        "invariants": INVARIANTS, "argv": sys.argv,
                        "runtime_note": "expected runtime unknown; estimate from utc_start/utc_end of commands and logs/"}
        self.flush()
    def flush(self):
        self.receipt["utc_updated"] = utc(); atomic_json(self.path, self.receipt)
    def stage(self, name, status="RUNNING"):
        self.receipt["stage"], self.receipt["status"] = name, status; self.flush()
    def run(self, name, argv, env=None):
        rec = {"stage": name, "argv": argv, "executed": self.execute, "utc_start": utc(), "utc_end": None, "returncode": None}
        self.receipt["commands"].append(rec); self.flush()
        if self.execute:
            logs = self.work / "logs"; logs.mkdir(exist_ok=True)
            with open(logs / f"{name}.stdout.log", "ab") as o, open(logs / f"{name}.stderr.log", "ab") as e:
                for h in (o, e): h.write(f"# {rec['utc_start']} argv={argv}\n".encode()); h.flush()
                rec["returncode"] = subprocess.run(argv, stdout=o, stderr=e, env=env).returncode  # no shell, no timeout
        rec["utc_end"] = utc(); self.flush()
        if self.execute and rec["returncode"] != 0:
            raise RuntimeError(f"{name} exited {rec['returncode']}; see logs/{name}.stderr.log")

def fasta_stats(fa):
    ids, lens, state = set(), [], {"id": None, "len": 0}
    def close():
        if state["id"] is not None:
            if state["len"] == 0: raise ValueError(f"empty sequence for contig {state['id']}")
            lens.append(state["len"])
    with open(fa, "r", encoding="ascii", errors="strict") as f:
        for line in f:
            line = line.rstrip("\r\n")
            if line.startswith(">"):
                close(); tok = line[1:].split()
                if not tok: raise ValueError("empty contig ID")
                if tok[0] in ids: raise ValueError(f"duplicate contig ID {tok[0]}")
                ids.add(tok[0]); state.update(id=tok[0], len=0)
            elif line:
                if state["id"] is None: raise ValueError("sequence before first header")
                if not set(line) <= FASTA_OK: raise ValueError(f"non-ACGTN characters in {state['id']}")
                state["len"] += len(line)
        close()
    if not lens: raise ValueError("no contigs in final.contigs.fa")
    lens.sort(reverse=True); total, acc, n50, n90 = sum(lens), 0, None, None
    for L in lens:
        acc += L
        if n50 is None and 2 * acc >= total: n50 = L
        if n90 is None and 10 * acc >= 9 * total: n90 = L
    ge = lambda t: sum(1 for L in lens if L >= t)
    st = {"total_contigs": len(lens), "total_bases": total, "n50": n50, "n90": n90, "max": lens[0], "min": lens[-1],
          "ge_500bp": ge(500), "ge_1kb": ge(1000), "ge_5kb": ge(5000), "ge_10kb": ge(10000),
          "quality_claim": "none; descriptive statistics only"}
    if st["ge_500bp"] == 0: raise ValueError("no contigs >= 500 bp")
    return st

def execute_recovery(a):
    work = Path(a.work_root)
    if work.exists(): sys.exit(f"ERROR: --work-root must not already exist: {work}")
    work.mkdir(parents=True, exist_ok=False)  # fresh directory is the duplicate-run guard (no fcntl; portable)
    work = work.resolve(); R = Runner(work, a.execute)
    try:
        R.stage("validate")
        if a.threads < 1: raise ValueError("--threads must be >= 1")
        if a.memory_bytes < 1_000_000_000: raise ValueError("--memory-bytes must be >= 1e9 (bytes, not GB)")
        r1, r2 = req_file(a.r1, "--r1"), req_file(a.r2, "--r2")
        if r1 == r2: raise ValueError("--r1 and --r2 resolve to the same file")
        rec_d, rec_q = req_file(a.download_receipt, "--download-receipt"), req_file(a.qc_receipt, "--qc-receipt")
        for p in (rec_d, rec_q): json.loads(p.read_text(encoding="utf-8"))  # syntax check only; schema not interpreted
        megahit, bpy = req_file(a.megahit, "--megahit"), req_file(a.busco_python, "--busco-python")
        bscript, lineage = req_file(a.busco_script, "--busco-script"), req_dir(a.lineage, "--lineage")
        cache = req_dir(a.download_path, "--download-path")
        ppath = req_dir(a.busco_pythonpath, "--busco-pythonpath") if a.busco_pythonpath else None
        R.receipt["params"] = {"source_accession": a.run_accession, "taxon_label": a.taxon_label,
            "threads": a.threads, "memory_bytes": a.memory_bytes, "min_contig_len": 500,
            "megahit_klist_mincount": "program defaults (not passed; not validated as optimal)",
            "r1": str(r1), "r2": str(r2), "lineage": str(lineage), "download_path": str(cache), "megahit": str(megahit),
            "busco_python": str(bpy), "busco_script": str(bscript), "busco_pythonpath": str(ppath) if ppath else None}
        R.stage("provenance")
        R.receipt["raw_before"] = {"r1": snap(r1), "r2": snap(r2)}
        R.receipt["hashes"].update({"download_receipt_sha256": sha256(rec_d), "qc_receipt_sha256": sha256(rec_q),
            "r1_sha256_as_stored": sha256(r1), "r2_sha256_as_stored": sha256(r2)}); R.flush()
        asm = work / "assembly"; contigs = asm / "final.contigs.fa"
        R.stage("megahit")
        R.run("megahit", [str(megahit), "-1", str(r1), "-2", str(r2), "--num-cpu-threads", str(a.threads),
                          "--memory", str(a.memory_bytes), "--min-contig-len", "500", "--out-dir", str(asm)])
        if a.execute:
            if not contigs.is_file() or contigs.stat().st_size == 0: raise RuntimeError(f"missing or empty {contigs}")
            R.stage("assembly_stats")
            R.receipt["assembly_stats"] = fasta_stats(contigs); R.receipt["hashes"]["assembly_sha256"] = sha256(contigs); R.flush()
        bout, run_name = work / "busco", a.taxon_label
        env = dict(os.environ); env["PATH"] = str(bpy.parent) + os.pathsep + env.get("PATH", "")
        if ppath:
            patch_file = ppath / "MEMORY_PATCH.json"
            patch = json.loads(patch_file.read_text())
            patch_code = ppath / "busco/busco_tools/miniprot.py"
            if sha256(patch_code) != patch["patched_sha256"]: raise ValueError("BUSCO memory patch hash mismatch")
            R.receipt["hashes"]["busco_memory_patch_sha256"] = sha256(patch_code)
        if ppath: env["PYTHONPATH"] = str(ppath) + (os.pathsep + env["PYTHONPATH"] if env.get("PYTHONPATH") else "")
        R.receipt["busco_env"] = {"PATH_prepend": str(bpy.parent), "PYTHONPATH_prepend": str(ppath) if ppath else None,
            "note": "patched BUSCO package (representation-only memory patch) is supplied externally; not patched here"}
        R.stage("busco")
        R.run("busco", [str(bpy), str(bscript), "-i", str(contigs), "-m", "genome", "-l", str(lineage), "--offline",
                        "--download_path", str(cache), "-o", run_name, "--out_path", str(bout), "-c", str(a.threads), "--opt-out-run-stats"], env=env)
        if a.execute:
            R.stage("busco_verify"); rd = bout / run_name
            summ = sorted(rd.glob("short_summary*.json"))
            if len(summ) != 1: raise RuntimeError(f"expected exactly one short_summary*.json in {rd}, found {len(summ)}")
            runs = sorted(p for p in rd.glob("run_*") if p.is_dir())
            if [p.name for p in runs] != [f"run_{lineage.name}"]:
                raise RuntimeError(f"expected exactly one run_{lineage.name} dir, found {[p.name for p in runs]}")
            ft, sc = runs[0] / "full_table.tsv", runs[0] / "busco_sequences" / "single_copy_busco_sequences"
            if not ft.is_file() or ft.stat().st_size == 0: raise RuntimeError(f"missing or empty {ft}")
            if not sc.is_dir(): raise RuntimeError(f"missing directory {sc}")
            s = json.loads(summ[0].read_text(encoding="utf-8")); res = s.get("results") if isinstance(s, dict) else None
            if not isinstance(res, dict): raise RuntimeError("short_summary JSON has no 'results' dict")
            n_sc = sum(1 for _ in sc.iterdir())
            R.receipt["busco"] = {"summary_path": str(summ[0]), "full_table_path": str(ft), "single_copy_dir": str(sc),
                "single_copy_files": n_sc, "results": res, "lineage_dataset": s.get("lineage_dataset"),
                "meaning": "execution record only; no marker admitted; 2661 exclusion and selection happen downstream",
                "empty_single_copy_note": "No single-copy sequence files; complete and duplicated counts are recorded in results" if n_sc == 0 else None}
            R.receipt["hashes"].update({"busco_summary_sha256": sha256(summ[0]), "busco_full_table_sha256": sha256(ft)})
        R.stage("finalize")
        R.receipt["raw_after"] = {"r1": snap(r1), "r2": snap(r2)}
        if R.receipt["raw_after"] != R.receipt["raw_before"]: raise RuntimeError("raw read size/mtime_ns changed during run")
        for k, p in (("download_receipt_sha256", rec_d), ("qc_receipt_sha256", rec_q)):
            if sha256(p) != R.receipt["hashes"][k]: raise RuntimeError(f"{k} changed during run")
        R.receipt["utc_end"] = utc()
        R.stage("done", "COMPLETE_EXECUTION_REVIEW_REQUIRED" if a.execute else "DRY_RUN_PLANNED")
        print(f"{R.receipt['status']}: {R.path}")
    except (Exception, KeyboardInterrupt) as e:
        R.receipt["error"] = f"{type(e).__name__}: {e}"; R.receipt["utc_end"] = utc()
        R.stage(R.receipt["stage"], "FAILED")
        print(f"FAILED at stage {R.receipt['stage']}: {R.receipt['error']} (receipt: {R.path})", file=sys.stderr); sys.exit(1)


def parse_arguments(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--config', type=Path, help='JSON config; relative paths resolve beside this file')
    path_keys = ('work_root', 'r1', 'r2', 'download_receipt', 'qc_receipt', 'megahit',
                 'busco_python', 'busco_script', 'lineage', 'download_path', 'busco_pythonpath')
    for key in path_keys:
        ap.add_argument('--' + key.replace('_', '-'))
    for key in ('taxon_label', 'run_accession'):
        ap.add_argument('--' + key.replace('_', '-'))
    ap.add_argument('--threads', type=int)
    ap.add_argument('--memory-bytes', type=int)
    ap.add_argument('--execute', action='store_true', help='Run MEGAHIT and BUSCO; default is read-only preflight')
    a = ap.parse_args(argv)
    config = json.loads(a.config.read_text()) if a.config else {}
    base = a.config.resolve().parent if a.config else Path.cwd()
    for key in (*path_keys, 'taxon_label', 'run_accession', 'threads', 'memory_bytes'):
        value = getattr(a, key)
        if value is None:
            value = config.get(key)
        if key == 'threads' and value is None:
            value = 8
        if key == 'memory_bytes' and value is None:
            value = 24_000_000_000
        if value is None and key != 'busco_pythonpath':
            ap.error('Missing configuration: ' + key)
        if key in path_keys and value is not None:
            path = Path(value).expanduser()
            value = str((path if path.is_absolute() else base / path).resolve())
        setattr(a, key, value)
    return a


def preflight(a):
    from _nuclear_cli import fresh_output, report_plan
    if not isinstance(a.threads, int) or a.threads < 1:
        raise ValueError('threads must be a positive integer')
    if not isinstance(a.memory_bytes, int) or a.memory_bytes < 1_000_000_000:
        raise ValueError('memory_bytes must be an integer >= 1e9')
    files = {key: req_file(getattr(a, key), key) for key in
             ('r1', 'r2', 'download_receipt', 'qc_receipt', 'megahit', 'busco_python', 'busco_script')}
    if files['r1'] == files['r2']:
        raise ValueError('Paired read paths must differ')
    for key in ('download_receipt', 'qc_receipt'):
        json.loads(files[key].read_text())
    for key in ('megahit', 'busco_python'):
        if not os.access(files[key], os.X_OK):
            raise ValueError('Executable permission is required: ' + key)
    dirs = {key: req_dir(getattr(a, key), key) for key in ('lineage', 'download_path')}
    if a.busco_pythonpath:
        patch_root = req_dir(a.busco_pythonpath, 'busco_pythonpath')
        patch = json.loads((patch_root / 'MEMORY_PATCH.json').read_text())
        if sha256(patch_root / 'busco/busco_tools/miniprot.py') != patch['patched_sha256']:
            raise ValueError('BUSCO memory patch hash mismatch')
        dirs['busco_pythonpath'] = patch_root
    if not a.taxon_label or Path(a.taxon_label).name != a.taxon_label or a.taxon_label in ('.', '..'):
        raise ValueError('taxon_label must be a simple nonempty output name')
    work = fresh_output(Path(a.work_root), [*files.values(), *dirs.values()])
    return report_plan('recover_nuclear.py', {**files, **dirs},
                       {'receipt': work / 'RUN_RECEIPT.json', 'assembly': work / 'assembly/final.contigs.fa',
                        'busco_full_table': work / 'busco' / a.taxon_label / ('run_' + Path(a.lineage).name) / 'full_table.tsv'},
                       ['MEGAHIT', 'BUSCO', 'BUSCO lineage dataset'],
                       threads=a.threads, memory_bytes=a.memory_bytes)


def plan(config_path):
    """Return a read-only, JSON-compatible artifact and dependency plan."""
    return preflight(parse_arguments(['--config',str(Path(config_path).resolve())]))


def main(argv=None):
    a = parse_arguments(argv)
    print(json.dumps(preflight(a), indent=2))
    if a.execute:
        execute_recovery(a)
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
