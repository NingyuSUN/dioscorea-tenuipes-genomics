"""Small, fail-closed I/O and unrooted topology helpers; no inference on import."""
from pathlib import Path
from io import StringIO
from decimal import Decimal
import hashlib, json, math, re, os
from Bio import Phylo

def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()

def save(path, value):
    path = Path(path)
    temporary = path.with_suffix(path.suffix + '.tmp')
    temporary.write_text(json.dumps(value, indent=2, allow_nan=False) + '\n')
    temporary.replace(path)

def fasta(path):
    rows, name = {}, None
    for line in Path(path).read_text().splitlines():
        if line.startswith('>'):
            name = line[1:].split()[0]
            if name in rows: raise ValueError('Duplicate FASTA ID: ' + name)
            rows[name] = ''
        elif line.strip():
            if name is None: raise ValueError('Sequence before FASTA header')
            rows[name] += line.strip().upper()
    if not rows or len({len(s) for s in rows.values()}) != 1:
        raise ValueError('Empty or unequal-length alignment')
    return rows

def write_fasta(path, rows):
    with Path(path).open('w') as out:
        for name, sequence in rows.items():
            out.write('>' + name + '\n')
            for i in range(0, len(sequence), 80): out.write(sequence[i:i+80] + '\n')

def topology(newick):
    tree = Phylo.read(StringIO(newick), 'newick')
    names = [t.name for t in tree.get_terminals()]
    if len(names) != len(set(names)): raise ValueError('Repeated tree tip')
    universe, splits = frozenset(names), set()
    for node in tree.get_nonterminals():
        side = frozenset(t.name for t in node.get_terminals())
        other = universe - side
        if min(len(side), len(other)) < 2: continue
        splits.add(side if min(universe) not in side else other)
    return universe, frozenset(splits)

def has_split(newick, members, expected_tips=None):
    universe, splits = topology(newick)
    if expected_tips is not None and universe != frozenset(expected_tips):
        raise ValueError('Tree label set differs from frozen alignment')
    side = frozenset(members)
    if not side <= universe: raise ValueError('Constraint label absent')
    if min(universe) in side: side = universe - side
    return side in splits

def constraint_newick(tips, members):
    tips, members = set(tips), set(members)
    if not members <= tips or len(members) < 2 or len(tips-members) < 2:
        raise ValueError('Invalid full-label split constraint')
    text = '((' + ','.join(sorted(members)) + '),' + ','.join(sorted(tips-members)) + ');'
    universe, splits = topology(text)
    if universe != tips or len(splits) != 1 or not has_split(text, members, tips):
        raise ValueError('Constraint encodes additional or missing splits')
    return text

def sig_radius(x,p=6):
    value=Decimal(str(x));assert value.is_finite()
    return Decimal(0) if value==0 else Decimal('.5')*Decimal(10)**(value.copy_abs().adjusted()-p+1)

def sitelh(path):
    tokens=Path(path).read_text().split();nt,ns=map(int,tokens[:2]);at=2;rows=[];bounds=[]
    for i in range(1,nt+1):
        assert tokens[at]==f'Tree{i}','Unexpected candidate row order'
        raw=tokens[at+1:at+1+ns];assert len(raw)==ns
        values=[float(x) for x in raw];assert all(math.isfinite(x) and x<=0 for x in values)
        rows.append((tokens[at],values));bounds.append(float(sum((sig_radius(x,6) for x in raw),Decimal(0))));at+=ns+1
    assert at==len(tokens),'Extra or missing site-likelihood tokens'
    return ns,rows,bounds

def user_tree_table(path):
    tail=Path(path).read_text().split('USER TREES',1)[1];result={}
    for line in tail.splitlines():
        fields=[x for x in line.split() if x not in ['+','-']]
        if len(fields)!=8 or not fields[0].isdigit():continue
        try:values=[float(v) for v in fields[1:]]
        except ValueError:continue
        result[int(fields[0])]=dict(zip(['logL','deltaL','bp_RELL','p_KH','p_SH','c_ELW','p_AU'],values))
    if not result:raise ValueError('Missing USER TREES table')
    return result

def topology_id(newick):
    tips,splits=topology(newick)
    canonical=json.dumps({'tips':sorted(tips),'splits':sorted(sorted(s) for s in splits)},separators=(',',':'))
    return hashlib.sha256(canonical.encode()).hexdigest()

def freeze_union(records,tips,constraints,maximum=18):
    grouped={}
    for row in records:
        newick=row['newick'];observed,_=topology(newick)
        assert observed==frozenset(tips),'Candidate label mismatch'
        grouped.setdefault(topology_id(newick),[]).append(row)
    assert len(records)<=maximum and len(grouped)>=2,'Candidate count outside frozen design'
    entries=[]
    for index,(identity,rows) in enumerate(sorted(grouped.items()),1):
        # Branch-length initializations selected deterministically, never by AU.
        newick=min(r['newick'] for r in rows)
        entries.append({'tree_index_1based':index,'topology_sha256':identity,'newick':newick,
                        'origins':[{k:v for k,v in r.items() if k!='newick'} for r in rows],
                        'satisfies_constraints':{name:has_split(newick,group,tips) for name,group in constraints.items()}})
    for name in constraints:assert any(e['satisfies_constraints'][name] for e in entries)
    return entries

def cgroup_check(c):
    rel=next(l[3:] for l in Path('/proc/self/cgroup').read_text().splitlines() if l.startswith('0::'))
    root=Path('/sys/fs/cgroup')/rel.lstrip('/')
    assert all((root/n).is_file() for n in ['memory.max','memory.swap.max','memory.events']),'Required memory cgroup absent'
    mem=(root/'memory.max').read_text().strip();swap=(root/'memory.swap.max').read_text().strip()
    assert mem!='max' and 0<int(mem)<=c['cgroup_memory_max_bytes'],'Use <=4GiB total cgroup'
    assert swap=='0','Swap must be disabled'
    available=sorted(os.sched_getaffinity(0));chosen=available[:c['cpu_affinity_max_cpus']]
    assert chosen;os.sched_setaffinity(0,chosen);assert sorted(os.sched_getaffinity(0))==chosen and len(chosen)<=4
    return {'path':str(root),'memory_max_bytes':int(mem),'swap_max_bytes':int(swap),'logical_cpu_ids':chosen,'cpu_max':(root/'cpu.max').read_text().strip() if (root/'cpu.max').exists() else None,'memory_events_before':(root/'memory.events').read_text()}

