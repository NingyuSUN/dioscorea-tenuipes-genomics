"""Portable-driver contracts using synthetic inputs and fake tools only."""
import contextlib
import hashlib
import importlib.util
import io
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import types
import unittest
from unittest import mock

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
SCRIPTS = ROOT / 'scripts/research'
TARGETS = {
    'full': 'population/historical/run_full_panel.py',
    'plink': 'population/historical/run_plink_arm.py',
    'compare': 'population/historical/compare_full_panel.py',
    'recover': 'nuclear/historical/recover_nuclear.py',
    'main': 'nuclear/historical/run_nuclear_main.py',
    'loci': 'nuclear/historical/prepare_loci.py',
    'translation': 'nuclear/historical/genomic_cds_translation.py',
}


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def load(key):
    path = SCRIPTS / TARGETS[key]
    spec = importlib.util.spec_from_file_location('portable_' + key, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def write_json(path, value):
    path.write_text(json.dumps(value))
    return path


def executable(path, source):
    path.write_text('#!' + sys.executable + '\n' + source)
    path.chmod(0o755)
    return path


class PortableContracts(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)

    def config(self, name, value):
        return write_json(self.root / (name + '.json'), value)

    def test_all_help_and_imports_are_data_free(self):
        for key, rel in TARGETS.items():
            with self.subTest(script=key):
                result = subprocess.run([sys.executable, str(SCRIPTS / rel), '--help'],
                                        cwd=self.root, capture_output=True, text=True)
                self.assertEqual(result.returncode, 0, result.stderr)
                self.assertIn('--execute', result.stdout)
                self.assertIn('--config', result.stdout)
                self.assertTrue(callable(load(key).plan))
        self.assertEqual(list(self.root.iterdir()), [])

    def assert_plan(self, module, config):
        result = module.plan(config)
        stream = io.StringIO()
        with contextlib.redirect_stdout(stream):
            self.assertEqual(module.main(['--config', str(config)]), 0)
        self.assertEqual(json.loads(stream.getvalue()), result)
        self.assertGreaterEqual(result['required_threads'], 1)
        return result

    def test_optimized_python_is_rejected(self):
        for key, relative in TARGETS.items():
            with self.subTest(script=key):
                result = subprocess.run([sys.executable, '-O', str(SCRIPTS / relative), '--help'],
                                        cwd=self.root, capture_output=True, text=True)
                self.assertNotEqual(result.returncode, 0)
                self.assertIn('disables required scientific assertions', result.stderr)

    def make_plink_inputs(self):
        run = self.root / 'fits'
        run.mkdir()
        samples = ['sample' + str(i) for i in range(14)]
        (run / 'paired.bed').write_bytes(b'\x6c\x1b\x01')
        (run / 'paired.bim').write_text('contig marker 0 1 A C\n')
        (run / 'paired.fam').write_text(''.join(f'{s} {s} 0 0 0 -9\n' for s in samples))
        outputs = [{'path': p.name, 'sha256': digest(p)} for p in run.iterdir()]
        write_json(run / 'RUN_RECEIPT.json', {'status': 'COMPLETED_PENDING_ROOT_SCIENTIFIC_AUDIT',
                   'selected_sites': 37566, 'samples': samples, 'output_files': outputs})
        return run, samples

    def test_plink_fake_tool_and_hash_rejection(self):
        run, samples = self.make_plink_inputs()
        tool = executable(self.root / 'fake_plink', '''import sys
from pathlib import Path
a=sys.argv
out=Path(a[a.index('--out')+1]);source=Path(a[a.index('--bfile')+1])
samples=[x.split()[1] for x in Path(str(source)+'.fam').read_text().splitlines()]
Path(str(out)+'.eigenvec').write_text('FID IID PCs\\n'+''.join(s+' '+s+' '+' '.join('1' if i==j else '0' for i in range(14))+'\\n' for j,s in enumerate(samples)))
Path(str(out)+'.eigenval').write_text('\\n'.join(str(14-i) for i in range(14))+'\\n')
Path(str(out)+'.rel').write_text('\\n'.join(' '.join(str(14-i) if i==j else '0' for i in range(14)) for j in range(14))+'\\n')
Path(str(out)+'.log').write_text('37566 variants and 14 people pass filters and QC.\\n')
''')
        config = self.config('plink', dict(run_root='fits', outdir='plink_out', plink='fake_plink', plink_sha256=digest(tool)))
        module = load('plink')
        result = self.assert_plan(module, config)
        self.assertEqual(result['status'], 'PREFLIGHT_PASSED')
        self.assertEqual(result['required_threads'],2)
        self.assertFalse((self.root / 'plink_out').exists())
        with contextlib.redirect_stdout(io.StringIO()):
            self.assertEqual(module.main(['--config', str(config), '--execute']), 0)
        receipt = json.loads((self.root / 'plink_out/RECEIPT.json').read_text())
        self.assertEqual(receipt['status'], 'COMPLETED_PENDING_INDEPENDENT_AUDIT')
        self.assertEqual(receipt['input_hashes_before'], receipt['input_hashes_after'])
        bad = self.config('bad_plink', dict(run_root='fits', outdir='bad_out', plink='fake_plink', plink_sha256='0'*64))
        with self.assertRaisesRegex(ValueError, 'SHA256 mismatch'):
            module.plan(bad)
        self.assertFalse((self.root / 'bad_out').exists())

    def test_full_panel_preflight_requires_exact_modules_and_runner(self):
        module = load('full')
        package = types.ModuleType('pcangsd')
        hashes = {}
        for name in ('shared', 'covariance', 'shared_cy', 'covariance_cy'):
            source = self.root / (name + '.py')
            source.write_text('# Synthetic dependency marker; no mathematical implementation.\n')
            child = types.ModuleType('pcangsd.' + name)
            child.__file__ = str(source)
            setattr(package, name, child)
            hashes[name] = digest(source)
        artifact = self.root / 'bound.dat'
        artifact.write_text('synthetic input identity\n')
        plan = self.config('fit_plan', dict(status='ROOT_APPROVED_FIXED_FULL_PANEL', settings=module.SETTINGS,
                  expected_runner_sha256=digest(SCRIPTS / TARGETS['full']), expected_reference_sha256=module.REFERENCE_SHA,
                  inputs={name: {'path': 'bound.dat', 'sha256': digest(artifact)} for name in ('markers','calls','baseline_mask','common_mask','source_header','extraction_receipt','bed','bim','fam')},
                  samples=['sample'+str(i) for i in range(14)], source_binding=dict(header_sha256='0'*64,VCF_path='synthetic.vcf',historical_body_sha256='0'*64,FAI_sha256='0'*64,source_samples=['sample'+str(i) for i in range(14)])))
        config = self.config('full', dict(plan=plan.name, outdir='fit_out', official_module_sha256=hashes))
        with mock.patch.dict(sys.modules, {'pcangsd': package}):
            result = self.assert_plan(module, config)
            self.assertEqual(result['status'], 'PREFLIGHT_PASSED')
            self.assertFalse((self.root / 'fit_out').exists())
            artifact.write_text('changed\n')
            with self.assertRaisesRegex(ValueError, 'Input identity changed'):
                module.plan(config)
            artifact.write_text('synthetic input identity\n')
            raw_plan=json.loads(plan.read_text());raw_plan['expected_runner_sha256']='0'*64;write_json(plan,raw_plan)
            with self.assertRaisesRegex(ValueError,'SHA256 mismatch'):
                module.plan(config)
            raw_plan['expected_runner_sha256']=digest(SCRIPTS/TARGETS['full']);write_json(plan,raw_plan)
        with mock.patch.dict(sys.modules, {'pcangsd': None}):
            with self.assertRaises(ModuleNotFoundError):
                module.plan(config)
        self.assertFalse((self.root/'fit_out').exists())
        self.assertEqual(module.SETTINGS['selected_sites'], 37566)
        self.assertEqual(module.SETTINGS['seed'], 20260924)
        self.assertEqual(module.parse_pl('0,20,40'), ((0,20,40), 'valid'))
        with self.assertRaises(ValueError):
            module.parse_pl('1,2,3')

    def test_comparison_preflight_and_artificial_matrices(self):
        module = load('compare')
        run = self.root / 'fits'; run.mkdir()
        plink = self.root / 'plink'; plink.mkdir()
        samples = ['group-' + str(i) for i in range(14)]
        values = np.arange(14, 0, -1, dtype=float)
        covariance = np.diag(values)
        arms = {}
        for rep in ('B_HC_e1e6', 'C_PL', 'B_HC_e1e4', 'B_official_HC_e1e6'):
            for rank in (2,1,3):
                name = rep + '_r' + str(rank)
                np.savez(run / (name + '.npz'), covariance=covariance)
                arms[name] = {'all_14_eigenvalues': values.tolist()}
        outputs = [{'path': p.name, 'sha256': digest(p)} for p in run.iterdir()]
        write_json(run / 'RUN_RECEIPT.json', dict(status='COMPLETED_PENDING_ROOT_SCIENTIFIC_AUDIT', selected_sites=37566,
                   arms=arms, samples=samples, output_files=outputs, called_cells=14, missing_cells=0))
        np.savetxt(plink / 'pca.rel', covariance); np.savetxt(plink / 'pca.eigenval', values)
        (plink / 'pca.eigenvec').write_text('FID IID PCs\n'+''.join(s+' '+s+' '+' '.join(str(x) for x in row)+'\n' for s,row in zip(samples,np.eye(14))))
        outputs = [{'path': p.name, 'sha256': digest(p)} for p in plink.iterdir()]
        write_json(plink / 'RECEIPT.json', dict(status='COMPLETED_PENDING_INDEPENDENT_AUDIT', outputs=outputs))
        config = self.config('comparison', dict(run_root='fits', plink_root='plink', outdir='comparison'))
        self.assertEqual(self.assert_plan(module, config)['status'], 'PREFLIGHT_PASSED')
        self.assertFalse((self.root / 'comparison').exists())
        with contextlib.redirect_stdout(io.StringIO()):
            module.main(['--config', str(config), '--execute'])
        result = json.loads((self.root / 'comparison/COMPARISON.json').read_text())
        self.assertEqual(len(result['arms']), 13)
        self.assertTrue(all(x['covariance_max_absolute_difference'] == 0 for x in result['comparisons']))

    def test_locus_preparation_preserves_synthetic_columns_and_missing_rows(self):
        module = load('loci')
        matrix = self.root / 'matrix'; matrix.mkdir(); (matrix / 'trimmed').mkdir()
        taxa = ['Dioscorea_tenuipes_historical', 'Dioscorea_gracillima', 'Oryza_sativa', 'Arabidopsis_thaliana', 'all_missing']
        seqs = {t: [] for t in taxa}; parts=[]; start=1
        for i in range(173):
            length=398 if i<115 else 397
            gene='g'+str(i); records={t: char*length for t,char in zip(taxa,'ACGT-')}
            path=matrix/'trimmed'/(gene+'.faa')
            path.write_text(''.join('>'+t+'\n'+s+'\n' for t,s in records.items()))
            parts.append(dict(gene=gene,start_1based=start,end_1based=start+length-1,length=length,sha256=digest(path)))
            start+=length
            for t,s in records.items():seqs[t].append(s)
        (matrix/'main.faa').write_text(''.join('>'+t+'\n'+''.join(s)+'\n' for t,s in seqs.items()))
        write_json(matrix/'summary.json',dict(genes=173,columns=68796,main_sha256=digest(matrix/'main.faa'),partitions=parts))
        tree=self.root/'reference.nwk';tree.write_text('(A:1,B:1);\n')
        config=self.config('loci',dict(matrix_dir='matrix',reference_tree='reference.nwk',reference_tree_sha256=digest(tree),outdir='loci_out'))
        self.assertEqual(self.assert_plan(module, config)['fixed_loci'],173)
        self.assertFalse((self.root/'loci_out').exists())
        with contextlib.redirect_stdout(io.StringIO()):module.main(['--config',str(config),'--execute'])
        receipt=json.loads((self.root/'loci_out/metadata/LOCUS_MANIFEST.json').read_text())
        self.assertEqual(receipt['columns'],68796)
        self.assertEqual(len(receipt['loci']),173)
        self.assertEqual([x['seed'] for x in receipt['loci']],list(range(20260928,20260928+173)))
        self.assertTrue(all(x['all_missing_taxa_omitted']==['all_missing'] and len(x['retained_taxa'])==4 for x in receipt['loci']))
        with self.assertRaisesRegex(ValueError,'already exists'):module.plan(config)

    def test_translation_preflight_is_read_only_and_hash_bound(self):
        module=load('translation')
        genome=self.root/'genome.fa';genome.write_text('>contig\nACGT\n')
        gff=self.root/'alignments.gff';gff.write_text('##gff-version 3\n')
        query=self.root/'queries.faa';query.write_text('>q0\n'+'A'*(122636-172)+'\n'+''.join(f'>q{i}\nA\n' for i in range(1,173)))
        write_json(self.root/'status.json',dict(exit_code=0,output_sha256=digest(gff)))
        write_json(self.root/'manifest.json',dict(query_sha256=digest(query)))
        config=self.config('translation',dict(genome='genome.fa',genome_sha256=digest(genome),gff='alignments.gff',query_fasta='queries.faa',alignment_status='status.json',query_manifest='manifest.json',outdir='translation_out'))
        result=self.assert_plan(module,config)
        self.assertEqual(result['status'],'PREFLIGHT_PASSED')
        self.assertFalse((self.root/'translation_out').exists())
        self.assertEqual(module.rc(b'ACGTN'),b'NACGT')
        genome.write_text('>contig\nTGCA\n')
        with self.assertRaisesRegex(ValueError,'SHA256 mismatch'):module.plan(config)

    def recovery_config(self):
        for name in ('r1.fastq','r2.fastq'):(self.root/name).write_text('@r\nACGT\n+\nIIII\n')
        for name in ('download.json','qc.json'):write_json(self.root/name,{'synthetic':True})
        (self.root/'embryophyta_odb10').mkdir();(self.root/'cache').mkdir()
        megahit=executable(self.root/'megahit', '''import sys
from pathlib import Path
a=sys.argv;out=Path(a[a.index('--out-dir')+1]);out.mkdir();(out/'final.contigs.fa').write_text('>contig\\n'+'A'*500+'\\n')
''')
        bpy=executable(self.root/'busco_python', '''import sys,json
from pathlib import Path
a=sys.argv;out=Path(a[a.index('--out_path')+1])/a[a.index('-o')+1];out.mkdir(parents=True)
(out/'short_summary.synthetic.json').write_text(json.dumps({'results':{'Complete':1}}))
run=out/'run_embryophyta_odb10';(run/'busco_sequences/single_copy_busco_sequences').mkdir(parents=True);(run/'full_table.tsv').write_text('# synthetic\\n')
''')
        (self.root/'busco_script').write_text('# synthetic script placeholder\n')
        return dict(work_root='recovered',r1='r1.fastq',r2='r2.fastq',download_receipt='download.json',qc_receipt='qc.json',megahit='megahit',busco_python='busco_python',busco_script='busco_script',lineage='embryophyta_odb10',download_path='cache',taxon_label='synthetic',run_accession='SYNTHETIC',threads=1,memory_bytes=1000000000)

    def test_recovery_fake_tools_and_default_no_output(self):
        module=load('recover');config=self.config('recover',self.recovery_config())
        result=self.assert_plan(module,config);self.assertEqual(result['status'],'PREFLIGHT_PASSED')
        self.assertFalse((self.root/'recovered').exists())
        with contextlib.redirect_stdout(io.StringIO()):module.main(['--config',str(config),'--execute'])
        receipt=json.loads((self.root/'recovered/RUN_RECEIPT.json').read_text())
        self.assertEqual(receipt['status'],'COMPLETE_EXECUTION_REVIEW_REQUIRED')
        self.assertEqual(receipt['assembly_stats']['total_bases'],500)
        self.assertEqual(receipt['raw_before'],receipt['raw_after'])

    def test_nuclear_main_preflight_uses_explicit_source_and_tools(self):
        module=load('main');self.recovery_config()
        for tool in ('mafft','iqtree'):executable(self.root/tool,'raise SystemExit(99)\n')
        annotation=self.root/'annotation.faa';annotation.write_text('>protein\nMALW\n')
        write_json(self.root/'annotation.json',dict(size=annotation.stat().st_size,computed_md5=hashlib.md5(annotation.read_bytes()).hexdigest()))
        taxa=[]
        for taxon in ('Dioscorea_tenuipes_historical','Dioscorea_gracillima','Oryza_sativa','Arabidopsis_thaliana'):
            run=self.root/taxon;(run/'run_embryophyta_odb10/busco_sequences/single_copy_busco_sequences').mkdir(parents=True)
            (run/'run_embryophyta_odb10/full_table.tsv').write_text('# Gene Start\tGene End\n')
            (run/'run_embryophyta_odb10/busco_sequences/single_copy_busco_sequences/synthetic.faa').write_text('>synthetic\nMALW\n')
            (run/'logs').mkdir();(run/'logs/miniprot_align_embryophyta_odb10_out.log').write_text('## synthetic raw evidence\n')
            write_json(run/'short_summary.synthetic.json',{'results':{}})
            taxa.append(dict(taxon=taxon,busco_run=taxon))
        taxa.append(dict(taxon='Dioscorea_oppositifolia_WT42',busco_run='generated_by_driver'))
        original_files = {p: digest(p) for row in taxa[:-1] for p in (self.root/row['taxon']).rglob('*') if p.is_file()}
        write_json(self.root/'taxa.json',taxa)
        config=self.config('nuclear',dict(work_root='nuclear_out',taxa_json='taxa.json',annotation_fasta='annotation.faa',annotation_metadata='annotation.json',busco_python='busco_python',busco_script='busco_script',lineage='embryophyta_odb10',download_path='cache',mafft='mafft',iqtree='iqtree',threads=1))
        result=self.assert_plan(module,config)
        self.assertEqual(result['status'],'PREFLIGHT_PASSED')
        self.assertIn(str(self.root/'iqtree'),result['tools'])
        self.assertEqual(result['required_threads'],1)
        self.assertFalse((self.root/'nuclear_out').exists())
        self.assertEqual(module.TAXA[-1]['busco_run'],str(self.root/'nuclear_out/busco/Dioscorea_oppositifolia_WT42'))
        executable(self.root/'busco_python', """import sys,json,hashlib
from pathlib import Path
a=sys.argv
if Path(a[1]).name=='build_busco_matrices.py':
 root=Path(a[a.index('--work-root')+1]);assert '--mafft' in a
 taxa=json.loads(Path(a[a.index('--taxa-json')+1]).read_text());names=[x['taxon'] for x in taxa]
 for row in taxa:
  run=Path(row['busco_run'])
  if run.exists():(run/'miniprot_stop_audit.json').write_text('{}')
 dest=root/'matrices/genome_main_no2661';dest.mkdir(parents=True)
 matrix=dest/'main.faa';matrix.write_text(''.join('>'+name+'\\n'+'M'*60+'\\n' for name in names))
 summary=dict(taxa=names,partitions=[dict(gene='synthetic_gene')],main_sha256=hashlib.sha256(matrix.read_bytes()).hexdigest(),excluded_taxa=[],genes=1,columns=60,required_focal_taxa_per_locus=['Dioscorea_tenuipes_historical','Dioscorea_gracillima'])
 (dest/'summary.json').write_text(json.dumps(summary))
""")
        executable(self.root/'iqtree', """import sys
from pathlib import Path
a=sys.argv;matrix=Path(a[a.index('-s')+1]);out=Path(a[a.index('--prefix')+1]);names=[line[1:] for line in matrix.read_text().splitlines() if line.startswith('>')]
tree='('+','.join(name+':1' for name in names)+');\\n'
Path(str(out)+'.treefile').write_text(tree);Path(str(out)+'.ufboot').write_text(tree*1000);Path(str(out)+'.log').write_text('Analysis results written to synthetic fixture\\n')
""")
        with contextlib.redirect_stdout(io.StringIO()):module.main(['--config',str(config),'--execute'])
        receipt=json.loads((self.root/'nuclear_out/NUCLEAR_RUN.json').read_text())
        self.assertEqual(receipt['status'],'INFERENCE_COMPLETE_SCIENTIFIC_REVIEW_AND_FIGURE_PENDING')
        self.assertEqual([x['stage'] for x in receipt['commands']],['busco_WT42','build_matrix','IQTREE_main'])
        command=receipt['commands'][-1]['argv']
        self.assertEqual(command[command.index('-seed')+1],'20260928')
        self.assertEqual(command[command.index('-B')+1],'1000')
        self.assertEqual(command[command.index('-T')+1],'1')
        self.assertEqual({p: digest(p) for p in original_files}, original_files)
        self.assertTrue(all(not (self.root/row['taxon']/'miniprot_stop_audit.json').exists() for row in taxa[:-1]))
        self.assertTrue(all(str(path) in result['inputs'] for path in original_files))
        staged=json.loads((self.root/'nuclear_out/metadata/INPUT_STAGING.json').read_text())
        self.assertEqual(len(staged['files']),len(original_files))


if __name__ == '__main__':
    unittest.main()
