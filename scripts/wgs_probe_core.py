"""Bounded download and FASTQ structural checks; never fabricate reads or scores."""
class LimitedReader:
 def __init__(self,stream,limit):
  import hashlib
  self.stream=stream;self.limit=limit;self.count=0;self.sha=hashlib.sha256()
 def read(self,n=-1):
  remaining=self.limit-self.count
  if remaining<=0:raise RuntimeError('Compressed prefix download budget exceeded')
  data=self.stream.read(min(remaining,n) if n>=0 else remaining);self.count+=len(data);self.sha.update(data);return data

def read_record(f):
 h=f.readline()
 if not h:return None
 s,p,q=f.readline(),f.readline(),f.readline()
 if not s or not p or not q:raise ValueError('Truncated FASTQ record')
 h,s,p,q=(x.rstrip(b'\r\n') for x in [h,s,p,q])
 if not h.startswith(b'@') or not p.startswith(b'+') or not s or len(s)!=len(q):raise ValueError('Invalid FASTQ structure')
 if any(x<33 or x>126 for x in q):raise ValueError('Invalid Phred+33 range')
 return h,s,p,q

def read_id(header):
 x=header.split()[0].removeprefix(b'@')
 return x[:-2] if x.endswith((b'/1',b'/2')) else x
