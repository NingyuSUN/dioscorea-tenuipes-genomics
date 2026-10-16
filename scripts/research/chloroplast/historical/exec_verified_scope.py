#!/usr/bin/env python3
"""Verify real cgroup memory controls and bind CPUs before an owned command."""
from pathlib import Path
import argparse
import json
import os


def validate(memory, swap, expected):
    if memory == 'max' or int(memory) != expected or swap != '0':
        raise ValueError('Requested memory/no-swap controls are not effective')


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--memory-bytes', type=int, required=True)
    p.add_argument('--cpus', type=int, required=True)
    p.add_argument('--receipt', type=Path, required=True)
    p.add_argument('command', nargs=argparse.REMAINDER)
    a = p.parse_args()
    command = a.command[1:] if a.command and a.command[0] == '--' else a.command
    if not command or not 1 <= a.cpus <= 4 or a.receipt.exists():
        raise ValueError('Require a command, 1-4 CPUs and a fresh receipt')
    rel = next(x[3:] for x in Path('/proc/self/cgroup').read_text().splitlines() if x.startswith('0::'))
    cg = Path('/sys/fs/cgroup') / rel.lstrip('/')
    memory = (cg / 'memory.max').read_text().strip()
    swap = (cg / 'memory.swap.max').read_text().strip()
    validate(memory, swap, a.memory_bytes)
    chosen = sorted(os.sched_getaffinity(0))[:a.cpus]
    if len(chosen) != a.cpus:
        raise ValueError('Fewer available CPUs than requested')
    os.sched_setaffinity(0, chosen)
    if sorted(os.sched_getaffinity(0)) != chosen:
        raise ValueError('CPU affinity was not applied')
    data = {'status': 'VERIFIED_BEFORE_EXEC', 'cgroup': str(cg),
            'memory_max_bytes': int(memory), 'swap_max_bytes': int(swap),
            'logical_cpu_ids': chosen, 'command': command,
            'cpu_max': (cg / 'cpu.max').read_text().strip() if (cg / 'cpu.max').exists() else None,
            'scope': 'This command and all descendants; coordinator is separately bounded.'}
    with a.receipt.open('x') as f:
        json.dump(data, f, indent=2)
        f.write('\n')
    os.execvpe(command[0], command, os.environ)


if __name__ == '__main__':
    main()
