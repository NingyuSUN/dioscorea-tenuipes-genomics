"""Synthetic checks of staging, failure preservation and workflow contracts.

Fake drivers only write JSON/text fixtures; no phylogenetic inference is run.
"""
import importlib.util
import json
import os
import shutil
import subprocess
import sys
import tempfile
import time
import unittest
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
HELPERS = REPO / "workflow/scripts"
sys.path.insert(0, str(HELPERS))
from research_contracts import nuclear_inputs, cp_tree_inputs
from research_run_guard import check_cp_resources, complete
from stage_research_inputs import stage

SMK = shutil.which("snakemake") or str(Path(sys.executable).with_name("snakemake"))


def nuclear_source(root):
    (root / "metadata").mkdir(parents=True)
    (root / "loci").mkdir()
    (root / "metadata/LOCUS_MANIFEST.json").write_text(json.dumps({"loci": [{"gene": "example"}]}))
    (root / "loci/example.faa").write_text('>synthetic\nACDE\n')


def cp_source(root):
    root.mkdir(parents=True)
    executable = root / "fake-iqtree"
    executable.write_text("never executed\n")
    (root / "model.nex").write_text("synthetic model\n")
    (root / "main95.fixed_columns.fa").write_text(">synthetic\nACGT\n")
    config = {"iqtree": str(executable), "threads": 4, "cpu_affinity_max_cpus": 4,
              "model_file": "model.nex", "sequential_panels": ["main95"]}
    (root / "CONFIG.json").write_text(json.dumps(config))
    (root / "PREPARED.json").write_text(json.dumps({"files": {name: {} for name in
        ["CONFIG.json", "model.nex", "main95.fixed_columns.fa"]}}))


class StagingTests(unittest.TestCase):
    def test_inputs_preserved_and_undeclared_files_not_copied(self):
        with tempfile.TemporaryDirectory() as d:
            src, dst = Path(d) / "input", Path(d) / "output"
            nuclear_source(src)
            (src / "loci/unused.faa").write_text("not a declared locus")
            (src / "old_results.txt").write_text("must not be copied")
            stage(src, dst, "nuclear")
            self.assertEqual((src / "loci/example.faa").read_bytes(), (dst / "loci/example.faa").read_bytes())
            self.assertFalse((dst / "old_results.txt").exists())
            self.assertFalse((dst / "loci/unused.faa").exists())
            self.assertIn("loci/example.faa", json.loads((dst / "STAGED.json").read_text())["sha256"])
            with self.assertRaisesRegex(ValueError, "not empty"):
                stage(src, dst, "nuclear")

    def test_missing_input_does_not_create_output(self):
        with tempfile.TemporaryDirectory() as d:
            dst = Path(d) / "out"
            with self.assertRaisesRegex(ValueError, "missing"):
                stage(Path(d) / "absent", dst, "nuclear")
            self.assertFalse(dst.exists())

    def test_nested_paths_and_symlink_alias_refused(self):
        with tempfile.TemporaryDirectory() as d:
            src = Path(d) / "input"
            nuclear_source(src)
            alias = Path(d) / "alias"
            alias.symlink_to(src, target_is_directory=True)
            for dst in [src, src / "out", alias / "out", Path(d)]:
                with self.assertRaisesRegex(ValueError, "non-nested"):
                    stage(src, dst, "nuclear")

    def test_input_contract_rejects_escaping_locus(self):
        with tempfile.TemporaryDirectory() as d:
            src = Path(d) / "input"
            nuclear_source(src)
            (src / "metadata/LOCUS_MANIFEST.json").write_text('{"loci":[{"gene":"../outside"}]}')
            with self.assertRaises(ValueError):
                nuclear_inputs(src)

    def test_cp_contract_includes_nested_files_config_and_binary(self):
        with tempfile.TemporaryDirectory() as d:
            src = Path(d) / "input"
            cp_source(src)
            paths = set(cp_tree_inputs(src))
            self.assertIn(str(src / "main95.fixed_columns.fa"), paths)
            self.assertIn(str(src / "CONFIG.json"), paths)
            self.assertIn(str(src / "fake-iqtree"), paths)
            with self.assertRaisesRegex(ValueError, "allocated"):
                check_cp_resources(src, 1)
            self.assertEqual(check_cp_resources(src, 4)["threads"], 4)

    def test_failed_native_status_never_gets_completion_marker(self):
        with tempfile.TemporaryDirectory() as d:
            p = Path(d)
            native = p / "RUN_STATUS.json"
            native.write_text('{"status":"FAILED_NO_AUTO_RESTART"}')
            with self.assertRaisesRegex(ValueError, "did not complete"):
                complete(native, p / "WORKFLOW_COMPLETE.json", "cp_trees")
            self.assertFalse((p / "WORKFLOW_COMPLETE.json").exists())
            self.assertTrue(native.exists())


@unittest.skipUnless(Path(SMK).is_file(), "Snakemake is needed for synthetic workflow checks")
class WorkflowContractTests(unittest.TestCase):
    def fixture(self, root, locus_failure=False, cp_failure=True):
        for folder in ["workflow/rules", "workflow/scripts", "scripts/research/nuclear", "scripts/research/chloroplast"]:
            (root / folder).mkdir(parents=True, exist_ok=True)
        for name in ["common.smk", "research.smk"]:
            shutil.copy2(REPO / "workflow/rules" / name, root / "workflow/rules" / name)
        for path in HELPERS.glob("*research*.py"):
            shutil.copy2(path, root / "workflow/scripts" / path.name)
        nuclear_source(root / "input/nuclear")
        cp_source(root / "input/cp")
        status = 'PARTIAL_FAILURE_REVIEW_REQUIRED' if locus_failure else 'COMPLETE_REVIEW_PENDING'
        done = [] if locus_failure else [{"gene": "example"}]
        failed = [{"gene": "example"}] if locus_failure else []
        stub = '''import argparse,json
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('--root',type=Path);a,_=p.parse_known_args()
r=a.root;(r/'gene_trees/example').mkdir(parents=True,exist_ok=True)
(r/'gene_trees/example/RECEIPT.json').write_text('partial or finished evidence')
(r/'LOCUS_ALL_STATUS.json').write_text(json.dumps(__RECEIPT__))
'''.replace('__RECEIPT__', repr({"status": status, "requested_loci": 1, "completed": done, "failed": failed}))
        (root / "scripts/research/nuclear/run_loci.py").write_text(stub)
        (root / "scripts/research/chloroplast/run_cp_trees.py").write_text('''import argparse,json
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('--out',type=Path);a,_=p.parse_known_args()
a.out.mkdir(parents=True,exist_ok=True)
(a.out/'CALLED').write_text('fake driver only')
(a.out/'RUN_STATUS.json').write_text(json.dumps(__CP_RECEIPT__))
'''.replace('__CP_RECEIPT__', repr({'status': 'FAILED_NO_AUTO_RESTART'} if cp_failure else {'status': 'FULL_COMPLETE_SCIENTIFIC_REVIEW_PENDING', 'jobs': [{'status': 'COMPLETE_REVIEW_PENDING'}]})))
        config = {"results_dir": str(root / "unused"), "tools": {}, "frozen_dirs": [], "research": {
            "outdir": str(root / "run"), "python": sys.executable, "modules": [],
            "nuclear_input_root": str(root / "input/nuclear"), "reference_tree": "unused.treefile",
            "spatial_plan": "unused.json", "balanced_plan": "unused.json",
            "chloroplast_ils_input_root": "unused", "cp_inputs": str(root / "input/cp")}}
        (root / "config.json").write_text(json.dumps(config))
        (root / "Snakefile").write_text('configfile: "config.json"\ninclude: "workflow/rules/common.smk"\ninclude: "workflow/rules/research.smk"\n')

    def run_smk(self, root, target, cores=1, dry=False):
        command = [SMK, target, "--cores", str(cores)] + (["--dry-run"] if dry else [])
        result = subprocess.run(command, cwd=root, text=True, capture_output=True,
                                env={**os.environ, "XDG_CACHE_HOME": str(root / "cache"), "PYTHONDONTWRITEBYTECODE": "1"})
        return result, result.stdout + result.stderr

    def test_zero_exit_failed_locus_keeps_partial_scientific_files(self):
        with tempfile.TemporaryDirectory(prefix="research failure ") as d:
            root = Path(d)
            self.fixture(root, locus_failure=True)
            result, log = self.run_smk(root, "research_locus_trees")
            self.assertNotEqual(result.returncode, 0, log)
            self.assertIn("Native driver did not complete", log)
            self.assertTrue((root / "run/nuclear/LOCUS_ALL_STATUS.json").exists())
            self.assertTrue((root / "run/nuclear/gene_trees/example/RECEIPT.json").exists())
            self.assertFalse((root / "run/nuclear/WORKFLOW_LOCUS_COMPLETE.json").exists())

    def test_changed_nested_locus_is_visible_to_dag(self):
        with tempfile.TemporaryDirectory(prefix="research change ") as d:
            root = Path(d)
            self.fixture(root)
            result, log = self.run_smk(root, "research_locus_trees")
            self.assertEqual(result.returncode, 0, log)
            locus = root / "input/nuclear/loci/example.faa"
            directory_time = locus.parent.stat().st_mtime_ns
            time.sleep(1.05)
            locus.write_text('>synthetic\nACDEF\n')
            self.assertEqual(directory_time, locus.parent.stat().st_mtime_ns)
            result, log = self.run_smk(root, "research_locus_trees", dry=True)
            self.assertNotIn("Nothing to be done", log)
            self.assertNotEqual(result.returncode, 0, log)
            self.assertIn("ProtectedOutputException", log)

    def test_cp_resource_rejection_happens_before_fake_driver(self):
        with tempfile.TemporaryDirectory(prefix="research cp resources ") as d:
            root = Path(d)
            self.fixture(root)
            result, log = self.run_smk(root, "research_cp_trees", cores=1)
            self.assertNotEqual(result.returncode, 0, log)
            self.assertIn("allocated threads", log)
            self.assertFalse((root / "run/chloroplast_trees/CALLED").exists())

    def test_changed_cp_payload_is_visible_to_dag(self):
        with tempfile.TemporaryDirectory(prefix="research cp change ") as d:
            root = Path(d)
            self.fixture(root, cp_failure=False)
            result, log = self.run_smk(root, "research_cp_trees", cores=4)
            self.assertEqual(result.returncode, 0, log)
            payload = root / "input/cp/main95.fixed_columns.fa"
            directory_time = payload.parent.stat().st_mtime_ns
            time.sleep(1.05)
            payload.write_text('>synthetic\nACGTA\n')
            self.assertEqual(directory_time, payload.parent.stat().st_mtime_ns)
            result, log = self.run_smk(root, "research_cp_trees", cores=4, dry=True)
            self.assertNotIn("Nothing to be done", log)
            self.assertNotEqual(result.returncode, 0, log)
            self.assertIn("ProtectedOutputException", log)

    def test_zero_exit_failed_cp_preserves_native_receipt(self):
        with tempfile.TemporaryDirectory(prefix="research cp failure ") as d:
            root = Path(d)
            self.fixture(root)
            result, log = self.run_smk(root, "research_cp_trees", cores=4)
            self.assertNotEqual(result.returncode, 0, log)
            self.assertIn("Native driver did not complete", log)
            self.assertTrue((root / "run/chloroplast_trees/RUN_STATUS.json").exists())
            self.assertTrue((root / "run/chloroplast_trees/CALLED").exists())
            self.assertFalse((root / "run/chloroplast_trees/WORKFLOW_COMPLETE.json").exists())


if __name__ == '__main__':
    unittest.main()
