import os
from snakemake.exceptions import WorkflowError

REPO = workflow.basedir
R = os.path.abspath(config["results_dir"])
LOGS = f"{R}/logs"


def script(name):
    return os.path.join(REPO, "scripts", name)


def helper(name):
    return os.path.join(REPO, "workflow", "scripts", name)


# Directories of completed runs. The numbered scripts refuse to overwrite
# an existing output directory; this check adds a second barrier so that no
# rule (and no --forcerun) can target them.
for frozen in map(os.path.abspath, config.get("frozen_dirs", [])):
    if R == frozen or R.startswith(frozen + os.sep):
        raise WorkflowError(
            f"results_dir {R} is inside the completed directory {frozen}. "
            "Choose a new, empty results_dir; completed evidence is read-only."
        )

TOOLS = config["tools"]
PY = TOOLS.get("python") or "python"
APY = TOOLS.get("analysis_python") or "python"
RSCRIPT = TOOLS.get("rscript") or "Rscript"
# Directory holding plink / admixture / vcftools; empty means "from PATH"
# (e.g. inside the conda env created by --sdm conda).
POPGEN_BIN = TOOLS.get("popgen_bin") or '$(dirname "$(command -v plink)")'

# Prefix for every population rule: optional PATH extension (e.g. a specific
# bcftools build), strict shell, headless matplotlib.
ENV = (
    (f'export PATH="{TOOLS["extra_path"]}:$PATH"; ' if TOOLS.get("extra_path") else "")
    + "set -euo pipefail; export MPLBACKEND=Agg PYTHONDONTWRITEBYTECODE=1; "
    + 'export MPLCONFIGDIR="${{TMPDIR:-/tmp}}/dt_mpl_$$"; '
)


def fresh(path):
    """Shell snippet: Snakemake pre-creates the parents of declared outputs, but
    scripts 32/33/34/44/53/54/55 call mkdir(exist_ok=False). Remove the tree
    only if it holds nothing but empty directories; otherwise stop, because a
    partial earlier run must be inspected, not overwritten."""
    return (
        f'if [ -d "{path}" ]; then find "{path}" -depth -type d -empty -delete; fi; '
        f'if [ -e "{path}" ]; then echo "Refusing to reuse non-empty {path}; '
        f'inspect or move it first." >&2; exit 1; fi; '
    )
