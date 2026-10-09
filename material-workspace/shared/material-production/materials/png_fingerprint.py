"""Bounded PNG byte/filtered-stream fingerprints; never changes image samples."""
import hashlib,struct,zlib
from .core import ContractError

class PNGSink:
 def __init__(self):
  self.digest=hashlib.sha256();self.bytes=0;self.buffer=bytearray();self.signature=False;self.filtered=hashlib.sha256();self.idat=hashlib.sha256();self.non_idat=hashlib.sha256();self.decompress=zlib.decompressobj();self.chunks=[];self.header=None;self.ended=False
 def write(self,b):
  self.digest.update(b);self.bytes+=len(b);self.buffer.extend(b)
  if not self.signature:
   if len(self.buffer)<8:return len(b)
   if bytes(self.buffer[:8])!=b'\x89PNG\r\n\x1a\n':raise ContractError('Bad PNG signature')
   del self.buffer[:8];self.signature=True
  while len(self.buffer)>=8:
   size=struct.unpack('!I',self.buffer[:4])[0]
   if len(self.buffer)<size+12:break
   raw=bytes(self.buffer[:size+12]);del self.buffer[:size+12];tag=raw[4:8];data=raw[8:-4];crc=struct.unpack('!I',raw[-4:])[0]
   if zlib.crc32(tag+data)&0xffffffff!=crc:raise ContractError('Bad PNG chunk CRC')
   self.chunks.append([tag.decode('ascii'),size])
   if tag==b'IDAT':
    self.idat.update(data)
    while data:
     decoded=self.decompress.decompress(data,1024*1024);self.filtered.update(decoded);data=self.decompress.unconsumed_tail
   else:self.non_idat.update(raw)
   if tag==b'IHDR':self.header=list(struct.unpack('!IIBBBBB',data))
   if tag==b'IEND':self.ended=True
  return len(b)
 def flush(self):pass
 def fingerprint(self):
  if self.buffer or not self.ended or not self.decompress.eof:raise ContractError('Incomplete PNG stream')
  return {'png_sha256':self.digest.hexdigest(),'bytes':self.bytes,'filtered_stream_sha256':self.filtered.hexdigest(),'deflate_stream_sha256':self.idat.hexdigest(),'non_idat_chunks_sha256':self.non_idat.hexdigest(),'IHDR':self.header,'chunks':self.chunks}

def fingerprint_file(path):
 sink=PNGSink()
 with open(path,'rb') as f:
  for b in iter(lambda:f.read(65536),b''):sink.write(b)
 return sink.fingerprint()

def classify(generated,canonical,pixel_equal):
 if generated['png_sha256']==canonical['png_sha256']:return 'exact_png_bytes'
 if not pixel_equal:return 'decoded_pixel_mismatch'
 if generated['filtered_stream_sha256']==canonical['filtered_stream_sha256'] and generated['non_idat_chunks_sha256']==canonical['non_idat_chunks_sha256']:return 'deflate_bitstream_only'
 return 'lossless_container_or_filter_difference'
