#!/usr/bin/env python3
"""Build the three repository pipeline diagrams using Python's standard library.

Run: python docs/figures/build_pipeline_figures.py
SVGs contain editable text and vector geometry; no network or raster assets.
These diagrams describe method/data flow, not a scientific validation result.
"""
from html import escape
from pathlib import Path

OUT = Path(__file__).resolve().parent
NAVY = '#17344F'
TEAL = '#087D86'
AMBER = '#A4660A'
MUTED = '#52687B'
BORDER = '#CCD7E0'
FLOW = '#647D92'

class Figure:
    def __init__(self, name, title, subtitle, height=1180):
        self.name, self.height = name, height
        self.items = [f'''<svg xmlns="http://www.w3.org/2000/svg" width="1680" height="{height}" viewBox="0 0 1680 {height}" role="img" aria-labelledby="title description">
<title id="title">{escape(title)}</title>
<desc id="description">{escape(subtitle)} Solid teal arrows denote the configured workflow. Dashed amber arrows denote separate historical scripts or templates. Dotted blue arrows denote research data flow without an automation claim.</desc>
<defs>
  <marker id="configured" markerWidth="10" markerHeight="10" refX="8" refY="5" orient="auto" markerUnits="userSpaceOnUse"><path d="M0,0 L9,5 L0,10 Z" fill="{TEAL}"/></marker>
  <marker id="historical" markerWidth="10" markerHeight="10" refX="8" refY="5" orient="auto" markerUnits="userSpaceOnUse"><path d="M0,0 L9,5 L0,10 Z" fill="{AMBER}"/></marker>
  <marker id="research" markerWidth="10" markerHeight="10" refX="8" refY="5" orient="auto" markerUnits="userSpaceOnUse"><path d="M0,0 L9,5 L0,10 Z" fill="{FLOW}"/></marker>
</defs>
<style>text {{font-family:Arial,Helvetica,sans-serif;fill:{NAVY};}} .subtitle {{fill:{MUTED};}} .note {{fill:{MUTED};}} </style>
<rect width="1680" height="{height}" fill="white"/>''']
        self.text(64, 66, title, 36, '700')
        self.text(64, 108, subtitle, 21, color=MUTED)
        self.items.append(f'<path d="M64 134 H1616" stroke="{BORDER}" stroke-width="2"/>')
        self.boxes=[]
    def text(self, x,y,s,size=21,weight='400',color=None,anchor='start'):
        color = color or NAVY
        self.items.append(f'<text x="{x}" y="{y}" font-size="{size}" font-weight="{weight}" text-anchor="{anchor}" style="fill:{color}">{escape(s)}</text>')
    def panel(self,x,y,w,h,label,fill='#F4F7FA',accent=NAVY):
        self.items.append(f'<rect x="{x}" y="{y}" width="{w}" height="{h}" rx="18" fill="{fill}"/>')
        self.text(x+24,y+34,label,18,'700',accent)
    def box(self,x,y,w,h,title,lines=(),kind='configured',tag=None):
        c = {'configured':TEAL,'historical':AMBER,'research':NAVY}[kind]
        self.items.append(f'<rect x="{x}" y="{y}" width="{w}" height="{h}" rx="12" fill="white" stroke="{BORDER}" stroke-width="1.6"/>')
        self.items.append(f'<path d="M{x+14} {y+1} H{x+w-14}" stroke="{c}" stroke-width="4" stroke-linecap="round"/>')
        pos=y+35
        self.text(x+22,pos,title,23,'700',c)
        for line in lines:
            pos+=29
            self.text(x+22,pos,line,20,color=MUTED)
        if tag:
            self.text(x+22,y+h-18,tag,16,'700',c)
        self.boxes.append((x,y,w,h,title,lines))
    def arrow(self,pts,kind='configured',head=True):
        c,dash={'configured':(TEAL,''),'historical':(AMBER,'10 7'),'research':(FLOW,'3 7')}[kind]
        d='M'+' L'.join(f'{x} {y}' for x,y in pts)
        self.items.append(f'<path d="{d}" fill="none" stroke="{c}" stroke-width="2.8" stroke-linejoin="round" stroke-linecap="round"'+(f' stroke-dasharray="{dash}"' if dash else '')+(f' marker-end="url(#{kind})"' if head else '')+'/>')
    def legend(self,y,all_styles=True):
        self.items.append(f'<path d="M64 {y-24} H1616" stroke="{BORDER}" stroke-width="1.5"/>')
        if all_styles:
            for x,kind,txt in [(64,'configured','Configured workflow'),(504,'historical','Historical scripts / templates'),(1070,'research','Research data flow only')]:
                self.arrow([(x,y),(x+55,y)],kind)
                self.text(x+72,y+7,txt,18)
        else:
            self.arrow([(64,y),(119,y)],'research')
            self.text(136,y+7,'Research data flow only',18)
            self.text(580,y+7,'Automation status and input contracts: see workflow documentation.',18,color=MUTED)
        self.text(64,y+45,'Method map only. Arrows do not certify completed runs, biological validity or a causal mechanism.',18,color=MUTED)
    def save(self):
        self.items.append('</svg>')
        (OUT/self.name).write_text('\n'.join(self.items)+'\n',encoding='utf-8')


def overview():
    f=Figure('pipeline_overview.svg','Dioscorea tenuipes · analysis map','Assembly, population genomics and organelle evidence connect to later comparative research.',1200)
    f.panel(64,165,1552,200,'01  NUCLEAR ASSEMBLY AND ANNOTATION','#F6F8FB')
    f.box(88,225,340,112,'Sequencing inputs',['Nanopore · Illumina · RNA-seq'],'historical')
    f.box(570,225,490,112,'Assembly + annotation',['Assembly, polishing, QC, gene models'],'historical',tag='SCRIPTS 01–15 · OPTIONAL TEMPLATE')
    f.box(1190,225,400,112,'Nuclear reference',['Reference sequence + annotation'],'historical')
    f.arrow([(428,280),(570,280)],'historical')
    f.arrow([(1060,280),(1190,280)],'historical')

    f.panel(64,395,1552,230,'02  POPULATION GENOMICS','#EFF8F8',TEAL)
    f.box(88,455,340,140,'Population VCF + audit',['Joint calling: script 21','Input checks: scripts 30–31'],'historical')
    f.box(505,455,300,140,'Genotype panels',['Masking + site selection'],'configured',tag='SCRIPT 32')
    f.box(882,455,320,140,'Population structure',['PCA · ADMIXTURE · FST'],'configured',tag='SCRIPTS 33–34')
    f.box(1278,455,312,140,'Related analyses',['Relationships: 44–49','Sample influence: 53–55'],'configured')
    f.arrow([(428,525),(505,525)],'historical')
    f.arrow([(805,525),(882,525)])
    f.arrow([(1202,525),(1278,525)])
    f.arrow([(655,455),(655,440),(1434,440),(1434,455)])
    f.arrow([(1390,337),(1390,380),(50,380),(50,525),(88,525)],'research')
    f.text(716,375,'reference for variant calling',16,color=MUTED)

    f.panel(64,655,1552,200,'03  ORGANELLE EVIDENCE','#FCF7ED',AMBER)
    f.box(88,715,340,112,'Reads + references',['Organelle input evidence'],'historical')
    f.box(570,715,490,112,'Recovery + evidence checks',['Assembly, annotation, read support'],'historical',tag='SCRIPTS 16–20 / 35–43')
    f.box(1190,715,400,112,'Organelle candidates',['Sequence and evidence summaries'],'historical')
    f.arrow([(428,770),(570,770)],'historical')
    f.arrow([(1060,770),(1190,770)],'historical')

    f.text(88,907,'LATER RESEARCH · data relationships; automation status in workflow documentation',18,'700')
    f.box(88,932,348,125,'Nuclear markers',['Alignment + gene trees','Gene concordance'],'research')
    f.box(474,932,348,125,'Chloroplast sensitivity',['Topology · source choice','Alignment / masking choices'],'research')
    f.box(860,932,348,125,'Matched comparison',['D. panthaica comparison','Matched marker definitions'],'research')
    f.box(1246,932,344,125,'Conditional ILS',['Scenario-based interpretation','Assumptions remain explicit'],'research')
    f.legend(1112)
    f.save()


def population():
    f=Figure('pipeline_population.svg','Population genomics · inputs, panels and sensitivity','Original contig coordinates are retained throughout the configured population workflow.',1190)
    f.panel(64,165,1552,235,'INPUT CONTRACT AND GENOTYPE SELECTION','#F6F8FB')
    f.box(88,228,400,142,'Joint population VCF',['Called VCF + sample manifest','Upstream variant calling: script 21'],'historical')
    f.box(638,228,400,142,'Input audit',['Allele / genotype comparison','Sample-quality diagnostics'],'historical',tag='SCRIPTS 30–31')
    f.box(1188,228,402,142,'Masked genotype panels',['Per-sample masks + site filters','Missing calls remain missing'],'configured',tag='SCRIPT 32')
    f.arrow([(488,299),(638,299)],'historical')
    f.arrow([(1038,299),(1188,299)],'historical')

    f.text(88,449,'CONFIGURED ANALYSIS MODULES',18,'700',TEAL)
    f.box(88,507,470,188,'Population structure',['PCA and ADMIXTURE','Pairwise and window FST','QC-aware summaries'],'configured',tag='SCRIPTS 33–34')
    f.box(605,507,470,188,'Individual relationships',['IBS distances · NJ trees · networks','Relatedness and sensitivity panels','Cross-implementation checks'],'configured',tag='SCRIPTS 44–49')
    f.box(1122,507,468,188,'Sample influence',['Leave-one-out and candidate subsets','PCA and FST sensitivity','Comparative summary tables'],'configured',tag='SCRIPTS 53–55')
    f.arrow([(1389,370),(1389,480),(323,480),(323,507)])
    f.arrow([(840,480),(840,507)])
    f.arrow([(1356,480),(1356,507)])

    f.panel(64,751,1552,265,'LATER RESEARCH · separate input contracts and execution status','#F6F8FB')
    f.box(88,820,680,160,'Spatial + balanced-sample FST',['Regional / leave-region-out sensitivity','Balanced subsets and estimator comparisons'],'research',tag='INPUT: RETAINED VARIANT / GENOTYPE PANELS')
    f.box(908,820,682,160,'Genotype-likelihood sensitivity',['Source VCF PL fields + matched SNP panel','Explicit three-genotype likelihood model checks'],'research',tag='PL LIKELIHOODS ARE NOT POSTERIOR PROBABILITIES')
    f.arrow([(323,695),(323,725),(50,725),(50,900),(88,900)],'research')
    f.text(1030,790,'Source VCF PL + matched panel',18,color=MUTED)
    f.arrow([(1249,798),(1249,820)],'research')
    f.text(70,1056,'Interpretation: assess sensitivity to data selection and sample composition; do not infer mechanism from a single plot.',19,color=MUTED)
    f.legend(1100)
    f.save()


def comparative():
    f=Figure('pipeline_comparative.svg','Comparative phylogeny · two evidence streams','Nuclear and chloroplast analyses remain distinct before paired, assumption-aware interpretation.',1190)
    f.panel(64,165,750,485,'NUCLEAR EVIDENCE','#EFF8F8',TEAL)
    f.panel(866,165,750,485,'CHLOROPLAST EVIDENCE','#F6F8FB')
    f.box(88,226,702,105,'Nuclear marker proteins',['Ortholog definitions and source selection'],'research')
    f.box(890,226,702,105,'Chloroplast sequences',['Sequence identity, source provenance and input selection'],'research')
    f.box(88,385,702,105,'Marker alignments + trees',['Alignment / trimming choices; gene and summary trees'],'research')
    f.box(890,385,702,105,'Chloroplast alignment + tree',['Alignment choices and topology estimates'],'research')
    f.arrow([(439,331),(439,385)],'research')
    f.arrow([(1241,331),(1241,385)],'research')
    f.box(88,543,335,82,'Gene concordance',['Across-marker agreement'],'research')
    f.box(455,543,335,82,'Matched comparison',['D. panthaica marker sets'],'research')
    f.box(890,543,335,82,'Topology sensitivity',['Alternative relationships'],'research')
    f.box(1257,543,335,82,'Input sensitivity',['Source choice + masking'],'research')
    f.arrow([(439,490),(439,518),(255,518),(255,543)],'research')
    f.arrow([(439,518),(623,518),(623,543)],'research')
    f.arrow([(1241,490),(1241,518),(1057,518),(1057,543)],'research')
    f.arrow([(1241,518),(1424,518),(1424,543)],'research')

    f.box(405,726,870,123,'Paired interpretation of nuclear and plastid evidence',['Compare supported relationships and sensitivity patterns.','Keep source, marker, alignment and topology uncertainty explicit.'],'research')
    f.arrow([(255,625),(255,681),(840,681),(840,726)],'research')
    f.arrow([(623,625),(623,681)],'research',False)
    f.arrow([(1057,625),(1057,681),(840,681)],'research',False)
    f.arrow([(1424,625),(1424,681),(1057,681)],'research',False)

    f.box(405,907,870,123,'Conditional incomplete lineage sorting (ILS) scenarios',['Simulations and analytic checks under specified assumptions.','Not a calibrated probability that ILS caused the observed pattern.'],'research')
    f.arrow([(840,849),(840,907)],'research')
    f.legend(1100,all_styles=False)
    f.save()

if __name__ == '__main__':
    overview()
    population()
    comparative()
    print('Wrote pipeline_overview.svg, pipeline_population.svg, pipeline_comparative.svg')
