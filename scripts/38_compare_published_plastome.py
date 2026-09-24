#!/usr/bin/env python3
"""Fetch the published accession and retain versioned evidence for plastome comparison."""
import argparse, hashlib, json, os, subprocess, urllib.request
from pathlib import Path

def main():
    ap=argparse.ArgumentParser();ap.add_argument('root',type=Path);a=ap.parse_args();root=a.root
    out=root/'published_comparison';out.mkdir(exist_ok=False)
    manifest=[]
    for kind,suffix in [('fasta','fasta'),('gbwithparts','gb')]:
        url='https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=PV612491&rettype='+kind+'&retmode=text'
        data=urllib.request.urlopen(url,timeout=90).read();p=out/('PV612491.'+suffix)
        assert data.startswith(b'>') if kind=='fasta' else data.startswith(b'LOCUS')
        p.write_bytes(data);manifest.append(dict(url=url,file=p.name,bytes=len(data),sha256=hashlib.sha256(data).hexdigest()))
    (out/'download_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    mm=os.environ.get('MINIMAP2','minimap2')
    for i in (1,2):
        cmd=[mm,'-t','4','-x','asm5','-c','--cs=long','--secondary=no',str(out/'PV612491.fasta'),str(root/'recovered'/f'Dtenuipes_Izu_reference_candidate.plastome.path{i}.fasta')]
        with (out/f'path{i}.paf').open('w') as f,(out/f'path{i}.log').open('w') as e:subprocess.run(cmd,stdout=f,stderr=e,check=True)
    print('Public accession downloaded and two plastome paths aligned.')
if __name__=='__main__':main()
