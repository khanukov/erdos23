"""Harden the Grotzsch finding: stronger MaxCut search (upper bound on the
envelope) plus a spectral lower bound (rigorous)."""
import itertools, json, pickle, random, sys, time
import numpy as np
sys.path.insert(0,"/home/user/erdos23/review")
from k7_leg_ceiling import (quotient, canonical, blow_up, compositions, multinomial,
                            invariant, isomorphic, edges_of, grotzsch, triangles)
B="/home/user/erdos23/.work/erdos23_full_replay_bundle"
DEN=181440; N9=1897
with open(f"{B}/public_flagsdp_data/cache_n9.pkl","rb") as fh:
    STATES=[list(m) for _o,m in pickle.load(fh)["states"]]
BUCK={}
for i,a in enumerate(STATES): BUCK.setdefault(invariant(a),[]).append(i)
def index9(a):
    for i in BUCK.get(invariant(a),()):
        if isomorphic(STATES[i],a): return i
    raise KeyError
EDGE=[];MASS=np.zeros((107,N9))
for r in range(107):
    with np.load(f"{B}/k7_compact_v1/root_{r:03d}.npz") as p:
        e=np.asarray(p["edge_raw"],dtype=np.int64); MASS[r]=np.asarray(p["mass_raw"],dtype=np.int64)
    s,a,b=np.nonzero(e); EDGE.append((s,a,b,e[s,a,b]))
MT=MASS.sum(axis=0)
W=[0.10173,0.09717,0.08691,0.10166,0.07628,0.08213,0.07344,0.09352,0.12084,0.10443,0.06189]
W=np.asarray(W); W/=W.sum()
GRO=grotzsch(); assert not triangles(GRO)
KEY={}; P=[];I=[];M=[]
for parts in compositions(9,11):
    k=canonical(*quotient(GRO,parts)); i=KEY.get(k)
    if i is None: i=KEY[k]=index9(blow_up(GRO,parts))
    P.append(parts); I.append(i); M.append(multinomial(parts))
P=np.asarray(P,dtype=float); I=np.asarray(I); M=np.asarray(M,dtype=float)
q=np.zeros(N9); np.add.at(q,I,M*np.exp(P@np.log(W)))
E=edges_of(GRO); pr=[W[u]*W[v] for u,v in E]
de=2*sum(pr); dm=2*min(sum(p for (u,v),p in zip(E,pr) if ((t>>u)&1)==((t>>v)&1)) for t in range(1<<11))
print("Grotzsch blow-up: d_edge=%.6f  d_mono=%.6f  (band upper end 0.3197, 2/25=0.08)"%(de,dm))
print("q9 support: %d states, sum %.12f"%(int((q>0).sum()),q.sum()))

def matrices():
    out=[]
    for r in range(107):
        s,a,b,v=EDGE[r]; w=v*q[s]; keep=w!=0
        if not keep.any(): out.append(None); continue
        aa,bb,ww=a[keep],b[keep],w[keep]
        act=sorted(set(aa.tolist())|set(bb.tolist())); rm={p:i for i,p in enumerate(act)}
        n=len(act); N=np.zeros((n,n))
        np.add.at(N,([rm[x] for x in aa.tolist()],[rm[x] for x in bb.tolist()]),ww)
        out.append(N)
    return out
MATS=matrices()
print("active profile classes per root: max %d"%max(m.shape[0] for m in MATS if m is not None))

def strong_min(N,starts,rng):
    n=N.shape[0]; sym=N+N.T; np.fill_diagonal(sym,0.0); diag=float(np.trace(N)); U=np.triu(N,1)
    best=None
    for t in range(starts):
        c=(np.asarray([rng.random()<0.5 for _ in range(n)],dtype=np.int8) if t else np.zeros(n,dtype=np.int8))
        improved=True
        while improved:                      # steepest-descent single flips
            improved=False
            same=(c[:,None]==c[None,:])
            g=(sym*same).sum(axis=1)-(sym*(~same)).sum(axis=1)
            i=int(np.argmax(g))
            if g[i]>1e-12: c[i]^=1; improved=True
        # Kernighan-Lin style sweep: force each vertex once, keep the best prefix
        for _ in range(3):
            order=list(range(n)); rng.shuffle(order); locked=np.zeros(n,dtype=bool)
            cur=c.copy(); bestc=c.copy()
            bv=diag+float((U*(c[:,None]==c[None,:])).sum())
            for _step in range(n):
                same=(cur[:,None]==cur[None,:])
                g=(sym*same).sum(axis=1)-(sym*(~same)).sum(axis=1)
                g[locked]=-np.inf
                i=int(np.argmax(g))
                if not np.isfinite(g[i]): break
                cur[i]^=1; locked[i]=True
                v=diag+float((U*(cur[:,None]==cur[None,:])).sum())
                if v<bv: bv, bestc = v, cur.copy()
            c=bestc
        val=diag+float((U*(c[:,None]==c[None,:])).sum())
        if best is None or val<best: best=val
    return best

def spectral_min(N):
    """min over rules >= diag + S/2 + (n/4) * lambda_min(B)  with B the symmetrised off-diagonal."""
    n=N.shape[0]; U=np.triu(N,1); Bm=U+U.T; diag=float(np.trace(N)); S=float(U.sum())
    lam=float(np.linalg.eigvalsh(Bm)[0])
    return diag+S/2.0+(n/4.0)*lam

rng=random.Random(2024)
t0=time.time(); env_hi=0.0; env_lo=0.0
for r in range(107):
    N=MATS[r]
    if N is None: continue
    env_hi+=strong_min(N,120,rng); env_lo+=spectral_min(N)
def leg(env): return (10.0/DEN)*(env-(2/25)*float(MT.dot(q)))/10.0
print("\nleg7 upper bound (strong local search over all rules) : %+.6e   [%.0fs]"%(leg(env_hi),time.time()-t0))
print("leg7 lower bound (spectral, rigorous)                : %+.6e"%leg(env_lo))
print("\nstrictly positive under the strong search? ", leg(env_hi)>0)
print("rigorously positive?                       ", leg(env_lo)>0)
