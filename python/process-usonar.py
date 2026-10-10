from matplotlib import pyplot as plt
import numpy as np
import scipy
import netCDF4 as nc
from tqdm import tqdm
import sys

display=False

try:
    dive=sys.argv[1]
except:
    dive="ml18_294a_SNR_echogram_1000.nc"
    
print("Processing %s" % dive)
ds = nc.Dataset(dive)

T=ds.variables["thresholded_data"][:].data.astype(int)
T[:,0:201]=0
T[:,-1]=0
A=ds.variables["Ar"][:].data

if display:
    fig, axes = plt.subplots(2, 1, figsize=(9, 12))
    axes[0].imshow(T[200:500,:].T,aspect='auto');
    axes[1].imshow(np.log(A[200:500,:]).T,aspect='auto');
    plt.show()


W=25
maxecholen=10
mincount=15
vtol=2
vscale=3.9e-3/40e-3
Craw=[]
C=[]
for i in tqdm(range(0,T.shape[0])):
    y=np.where(T[i,200:])[0]
    x=np.ones(y.shape[0])*i
    Craw.append((i,np.vstack([x,y]).T))

    Ci=[]
    y=scipy.ndimage.label(T[i,200:])
    l,n=scipy.ndimage.label(T[i,200:])
    for j in range(n+1):
        blk=np.where(l==j)[0]
        if blk.shape[0]>maxecholen:
            continue
        Ablk=A[i,200+blk[0]:200+blk[-1]+1]
        y=np.sum(blk*Ablk)/np.sum(Ablk)
        Ci.append([i,y])
    C.append((i,np.array(Ci)))

#breakpoint()

candidates=[]
for i in tqdm(range(len(C)-1)):
    candidate_i = []
    reali,C0 = C[i]
    if len(C0)==0:
        # assert(len(Craw[i][1])==0)
        candidates.append((reali,np.array([]).astype(int)))
        continue
    _,C1 = C[i+1]
    if len(C1)==0:
        # assert(len(Craw[i+1][1])==0)
        candidates.append((reali,np.array([]).astype(int)))
        continue
    c0=np.dot(np.array([C0[:,1]]).T,np.ones((1,C1.shape[0])))
    c1=np.dot(np.ones((C0.shape[0],1)),np.array([C1[:,1]]))
    c=(c0-c1).T
    V,U = np.where(np.logical_and(c < 40, c>=0))
    for u,v in zip(list(U),list(V)):
        candidate_i.append(np.hstack([C0[u,:],C1[v,:],np.array([C0[u,1] - C1[v,1]])]))
    candidates.append((reali,np.array(candidate_i).astype(int)))

R=[]
Rraw=[]
for i0 in tqdm(range(len(candidates)-W)):
    reali = (candidates[i0][0]+candidates[i0+W][0])//2
    try:
        candidates_i0=np.vstack([c for _,c in candidates[i0:i0+W] if len(c.shape)>1]).astype(int)
    except:
        continue
    unique, counts = np.unique(candidates_i0[:,4], return_counts=True)
    V=np.argmax(counts)
    Rraw.append(np.array([reali,unique[V]*vscale,counts[V]]))
    if counts[V] > mincount:
        inliers=np.where(abs(candidates_i0[:,4]-unique[V]) < vtol)[0]
        vmean=np.mean(candidates_i0[inliers,4])
        R.append(np.array([reali,vmean*vscale,counts[V]]))
        if (i0 % 1000 == 0):
            Ca=np.vstack([c for _,c in C[i0:i0+W] if len(c)>0])
            plt.scatter(Ca[:,0],Ca[:,1])
            plt.plot([C[i0][0],C[i0+W][0]],[500,500-(C[i0+W][0]-C[i0][0])*vmean],'-r')
            plt.grid()
            plt.title("dive %s\nindex %d:%d Estimated velocity: %.2f" % (dive,i0,reali,vmean))
            plt.savefig(dive[:-3]+("-model-%06d.png"%reali))
            if display:
                plt.show()
            else:
                plt.close()
R=np.array(R)
Rraw=np.array(Rraw)
fig, axes = plt.subplots(2, 1, figsize=(36,36))
axes[0].plot(Rraw[:,0],Rraw[:,1],'-b');
axes[0].plot(R[:,0],R[:,1],'+g');
axes[0].grid()
axes[0].set_title("Velocity [m/s] for dive %s" % dive)
axes[1].plot(Rraw[:,0],Rraw[:,2],'+-b');
axes[1].grid()
axes[1].set_title("Support [count] for dive %s" % dive)
fig.savefig(dive[:-3]+"-vel.png")
if display:
    plt.show()
else:
    plt.close()

