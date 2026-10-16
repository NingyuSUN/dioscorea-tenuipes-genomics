"""Check every configured historical CLI without loading research inputs."""
import json
import os
import subprocess
import sys
import tempfile
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]


def main():
    catalog = json.loads((REPO / 'workflow/research/PORTABLE_STAGES.json').read_text())
    assert len(catalog['stages']) == 22, 'The reviewed historical scope must remain 22 scripts'
    with tempfile.TemporaryDirectory(prefix='portable CLI ') as directory:
        root = Path(directory)
        env = dict(os.environ, PYTHONDONTWRITEBYTECODE='1', MPLBACKEND='Agg')
        for item in catalog['stages']:
            script = REPO / item['script']
            result = subprocess.run([sys.executable, str(script), '--help'], cwd=root, env=env,
                                    text=True, capture_output=True, timeout=30)
            assert result.returncode == 0, (item['script'], result.stdout, result.stderr)
            assert '--config' in result.stdout and '--execute' in result.stdout, item['script']
            config = REPO / item['example_config']
            assert isinstance(json.loads(config.read_text()), dict), config
        assert not list(root.iterdir()), 'CLI help created files in the working directory'
    print('22/22 portable entry points: safe help, explicit execution flag, valid configuration examples. No analysis executed.')


if __name__ == '__main__':
    main()
