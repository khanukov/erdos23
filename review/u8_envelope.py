#!/usr/bin/env python3
"""Does the 8-root envelope succeed where the 7-root one fails?

U_k = sum over root types of the best profile-rule cut is squeezed between
d_mono (below, by validity) and whatever the certificate needs (above).  Larger
roots express more colourings, so U_k decreases toward d_mono as k grows: an
8-root envelope is TIGHTER than a 7-root one, not looser.  The 7-root envelope
going above 2/25 at the Grotzsch point therefore says nothing about the 8-root
one, and this script measures it.

The shipped 8-root decomposition is the broken one (see INDEPENDENT_REVIEW.md
section 9), so U8 is computed twice: as shipped, and after averaging the
recorded profile pairs over Aut(tau) -- the repair verified in
review/verify_aut_fix.py.  Both are computed by local search over rules, which
UPPER-bounds each per-root minimum, hence upper-bounds U8: the direction that
matters, since U8_upper < 2/25 implies U8 < 2/25.

Two checks per point: `U8 >= d_mono` (validity, required of any envelope) and
`U8 < 2/25` (usefulness, required for the certificate to say anything).

Caveat: the root graphs are not stored in the decomposition and are recovered
from the recorded profiles, so the automorphism groups are recovered too.  The
repair is a patch on the shipped tables, not a regeneration.
"""
from fractions import Fraction
import numpy as np
from scipy.sparse import coo_matrix, csr_matrix
sys.path.insert(0,"/home/user/erdos23/review")
from k7_leg_ceiling import (quotient, canonical, blow_up, compositions, multinomial,
                            invariant, isomorphic, edges_of, grotzsch, petersen, triangles)
from envelope_maxcut_bound import best_cut, exact_cut
B="/home/user/erdos23/.work/erdos23_full_replay_bundle"
N9,N10=1897,12172
with open(f"{B}/public_flagsdp_data/cache_n9.pkl","rb") as fh: CACHE=pickle.load(fh)
STATES=[list(m) for _o,m in CACHE["states"]]
BUCK={}
for i,a in enumerate(STATES): BUCK.setdefault(invariant(a),[]).append(i)
def index9(a):
    for i in BUCK.get(invariant(a),()):
        if isomorphic(STATES[i],a): return i
    raise KeyError
lift=np.load(f"{B}/public_flagsdp_data/c5lift_cache.npz",allow_pickle=True)
cnt=np.rint(np.asarray(lift["Dval"])*10).astype(np.int64)
Dcsc=csr_matrix((cnt,(lift["Drow"],lift["Dcol"])),shape=(N9,N10),dtype=np.int64).tocsc()
COL={}
for j in range(N10):
    s,e=Dcsc.indptr[j],Dcsc.indptr[j+1]
    COL[tuple(sorted(zip(Dcsc.indices[s:e].tolist(),Dcsc.data[s:e].tolist())))]=j

def q10_of(base,W):
    h=len(base); key9={}
    idx9={}
    for parts in compositions(9,h):
        k=canonical(*quotient(base,parts)); i=key9.get(k)
        if i is None: i=key9[k]=index9(blow_up(base,parts))
        idx9[parts]=i
    q=np.zeros(N10); tot=0.0
    for parts in compositions(10,h):
        prof={}
        for i in range(h):
            if parts[i]:
                d=list(parts); d[i]-=1; j=idx9[tuple(d)]
                prof[j]=prof.get(j,0)+parts[i]
        j10=COL[tuple(sorted(prof.items()))]
        w=multinomial(parts)*float(np.prod(np.power(W,parts)))
        q[j10]+=w; tot+=w
    return q/tot

def dens(base,W):
    E=edges_of(base); pr=[W[u]*W[v] for u,v in E]
    return 2*sum(pr), 2*min(sum(p for (u,v),p in zip(E,pr) if ((t>>u)&1)==((t>>v)&1)) for t in range(1<<len(base)))

with open(f"{B}/public_flagsdp_data/u8_decomp.pkl","rb") as fh: ED=pickle.load(fh)
EP={r: tuple(sorted({tuple(p) for p in ED["Rprofiles"][r]},key=lambda p:(len(p),p))) for r in range(410)}
with open(f"{B}/public_flagsdp_data/u8_decomp_all.pkl","rb") as fh: AD=pickle.load(fh)
AP={r: tuple(sorted({tuple(p) for p in AD["Rprofiles"][r]},key=lambda p:(len(p),p))) for r in range(410)}
def recover(profiles):
    tog=set()
    for p in profiles:
        for i,j in itertools.combinations(p,2): tog.add((i,j))
    adj=[0]*8
    for i,j in itertools.combinations(range(8),2):
        if (i,j) not in tog: adj[i]|=1<<j; adj[j]|=1<<i
    return adj
def auts(adj):
    d=[bin(m).count("1") for m in adj]; out=[]; img=[-1]*8
    def go(k,used):
        if k==8: out.append(tuple(img)); return
        for w in range(8):
            if (used>>w)&1 or d[w]!=d[k]: continue
            if all(((adj[k]>>j)&1)==((adj[w]>>img[j])&1) for j in range(k)):
                img[k]=w; go(k+1,used|(1<<w)); img[k]=-1
    go(0,0); return out
ROOTS=[recover(AP[r]) for r in range(410)]   # all-profiles give the better recovery
AUT=[auts(a) for a in ROOTS]
ei=[{p:i for i,p in enumerate(EP[r])} for r in range(410)]
rows=[[] for _ in range(410)];cols=[[] for _ in range(410)]
for s,cs in enumerate(ED["decomp"]):
    for r,pa,pb in cs:
        r=int(r);w=len(EP[r])
        rows[r].append(ei[r][tuple(pa)]*w+ei[r][tuple(pb)]); cols[r].append(s)
EPM=[]
for r in range(410):
    w=len(EP[r])
    M=coo_matrix((np.ones(len(rows[r]),dtype=np.int16),(rows[r],cols[r])),shape=(w*w,N10),dtype=np.int16).tocsr()
    M.sum_duplicates(); EPM.append(M)

# --- precompute, per root: the Aut-closure index and the permutation table ---
CLOSURE=[];PERM=[]
for r in range(410):
    closure=list(EP[r]); pos={p:i for i,p in enumerate(closure)}
    frontier=list(EP[r])
    while frontier:
        nxt=[]
        for g in AUT[r]:
            for p in frontier:
                im=tuple(sorted(g[x] for x in p))
                if im not in pos:
                    pos[im]=len(closure); closure.append(im); nxt.append(im)
        frontier=nxt
    CLOSURE.append(closure)
    PERM.append(np.asarray([[pos[tuple(sorted(g[x] for x in p))] for p in closure]
                            for g in AUT[r]],dtype=np.int64))
print("closure widths: max %d ; |Aut| max %d"%(max(len(c) for c in CLOSURE),
                                               max(len(a) for a in AUT)))

def u8_of(q,symmetrise,starts=40,seed=5):
    rng=random.Random(seed); tot=0.0
    for r in range(410):
        w=len(EP[r]); V=(EPM[r]@q).reshape(w,w)/90.0
        if V.sum()<=0: continue
        if symmetrise:
            n=len(CLOSURE[r]); G=PERM[r]; k=G.shape[0]
            aa,bb=np.nonzero(V); vv=V[aa,bb]
            idx=(G[:,aa]*n+G[:,bb]).ravel()
            wts=np.tile(vv/k,k)
            V=np.bincount(idx,weights=wts,minlength=n*n).reshape(n,n)
        V=(V+V.T)/2.0
        diag=float(np.trace(V)); U=np.triu(V,1); S=float(U.sum())
        Wm=U+U.T
        act=np.nonzero(Wm.sum(axis=1)+np.diag(V))[0]
        if len(act)>1:
            sub=Wm[np.ix_(act,act)]
            cut=exact_cut(sub) if len(act)<=18 else best_cut(sub,starts,rng)
        else:
            cut=0.0
        tot+=diag+2.0*(S-cut)
    return tot

LO,HI=0.2486,0.3197
def maximise_dmono(base,seed,restarts=5,steps=300):
    h=len(base); rng=random.Random(seed); best=(-1,None)
    for t in range(restarts):
        W=np.full(h,1.0/h) if t==0 else np.asarray([rng.random()+1e-3 for _ in range(h)]); W/=W.sum()
        de,dm=dens(base,W); cur=dm if LO<=de<=HI else -1
        sc=0.3
        for st in range(steps):
            i,j=rng.randrange(h),rng.randrange(h)
            if i==j: continue
            T=W.copy(); d=sc*rng.random()*T[i]; T[i]-=d; T[j]+=d
            de,dm=dens(base,T)
            if not LO<=de<=HI: continue
            if dm>cur: W,cur=T,dm
            if st%50==49: sc*=0.65
        if cur>best[0]: best=(cur,W.copy())
    return best
GRO=grotzsch(); PET=petersen()
for g in (GRO,PET): assert not triangles(g)
cands=[("Grotzsch(leg7-worst)",GRO,np.asarray([0.10173,0.09717,0.08691,0.10166,0.07628,0.08213,0.07344,0.09352,0.12084,0.10443,0.06189]))]
for nm,g in [("Grotzsch",GRO),("Petersen",PET)]+[("tf9#%d"%i,STATES[i]) for i in (1866,1329,906,1607,1791,171)]:
    v,W=maximise_dmono(g,hash(nm)&0xffff if False else sum(map(ord,nm)))
    if v>0: cands.append((nm+"(dmono-max)",g,W))
print("\n%-24s %8s %9s %11s %11s %10s"%("graphon","d_edge","d_mono","U8 broken","U8 repaired","leg8"))
for nm,g,W in cands:
    W=np.asarray(W,dtype=float); W=W/W.sum()
    de,dm=dens(g,W)
    q=q10_of(g,W)
    raw=u8_of(q,False,starts=12); fix=u8_of(q,True,starts=12)
    print("%-24s %8.4f %9.6f %11.6f %11.6f %10.5f   valid:%s  <2/25:%s"%(
        nm,de,dm,raw,fix,fix-0.08,fix>=dm-1e-9,fix<0.08),flush=True)
