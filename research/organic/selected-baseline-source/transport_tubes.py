"""Continuously transported circular tube frames, without tangent threshold flips."""
import math
from mathutils import Vector
from expansion.geometry import Tubes as LegacyTubes
class TransportTubes(LegacyTubes):
 def path(self,points,radius,sides=6,closed=False):
  start=len(self.v);n=len(points);self.components+=1;previous=None
  for i,p in enumerate(points):
   p=Vector(p);before=Vector(points[(i-1)%n] if closed or i else points[i]);after=Vector(points[(i+1)%n] if closed or i<n-1 else points[i]);tangent=(after-before).normalized()
   if previous is None:
    axes=[Vector((1,0,0)),Vector((0,1,0)),Vector((0,0,1))];ref=min(axes,key=lambda a:abs(a.dot(tangent)));a=tangent.cross(ref).normalized()
   else:
    a=previous-tangent*previous.dot(tangent)
    if a.length<1e-9:
     axes=[Vector((1,0,0)),Vector((0,1,0)),Vector((0,0,1))];ref=min(axes,key=lambda v:abs(v.dot(tangent)));a=tangent.cross(ref)
    a.normalize()
    if a.dot(previous)<0:a=-a
   b=tangent.cross(a).normalized();previous=a
   r=radius[i] if isinstance(radius,list) else radius
   for j in range(sides):
    self.v.append(tuple(p+r*(a*math.cos(j*math.tau/sides)+b*math.sin(j*math.tau/sides))));self.uv.append((j/sides,i/max(1,n-1)))
  for i in range(n if closed else n-1):
   for j in range(sides):self.f.append((start+i*sides+j,start+i*sides+(j+1)%sides,start+((i+1)%n)*sides+(j+1)%sides,start+((i+1)%n)*sides+j))
  if not closed:
   self.f.append(tuple(start+j for j in reversed(range(sides))));self.f.append(tuple(start+(n-1)*sides+j for j in range(sides)))
