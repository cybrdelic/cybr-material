"""Loop-domain open-tube UV chart with nondegenerate planar cut caps."""
import numpy as np
def open_tube_chart(points, sides=8):
    p=np.asarray(points,dtype=float)
    if p.ndim!=2 or p.shape[1]!=3 or len(p)<2 or sides<3:raise ValueError("Invalid open tube")
    segment=np.linalg.norm(np.diff(p,axis=0),axis=1)
    if np.any(segment<=0):raise ValueError("Repeated centerline sample")
    length=np.r_[0.,np.cumsum(segment)];v=length/length[-1]
    uv=[]
    for i in range(len(p)-1):
        for j in range(sides):
            uv.extend([(j/sides,v[i]),((j+1)/sides,v[i]),((j+1)/sides,v[i+1]),(j/sides,v[i+1])])
    theta=np.arange(sides)*2*np.pi/sides
    disk=np.stack([.5+.5*np.cos(theta),.5+.5*np.sin(theta)],axis=1)
    uv.extend(disk[::-1]);uv.extend(disk)
    return np.asarray(uv,np.float32),float(length[-1])

