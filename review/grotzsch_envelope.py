"""At the Grotzsch point the certificate's K7 pool gives leg7 > 0.
Is that an artefact of the small rule pool, or does the true envelope
(minimum over ALL profile rules) also exceed 2/25 there?"""
import itertools, json, pickle, random, sys, time
from math import factorial
import numpy as np
sys.path.insert(0,"/home/user/erdos23/review")
from k7_leg_ceiling import (quotient, canonical, blow_up, compositions, multinomial,
                            invariant, isomorphic, edges_of, grotzsch, petersen, cycle,
                            triangles, degrees)
B="/home/user/erdos23/.work/erdos23_full_replay_bundle"
LO,HI=0.2486,0.3197; DEN=181440; N9=1897
with open(f"{B}/public_flagsdp_data/cache_n9.pkl","rb") as fh:
    STATES=[list(m) for _o,m in pickle.load(fh)["states"]]
BUCK={}
for i,a in enumerate(STATES): BUCK.setdefault(invariant(a),[]).append(i)
def index9(a):
    for i in BUCK.get(invariant(a),()):
        if isomorphic(STATES[i],a): return i
    raise KeyError
cert=json.loads(open(f"{B}/erdos23_global_exact_dual.json").read())
EDGE=[];MASS=np.zeros((107,N9));POOL={}
for r in range(107):
    with np.load(f"{B}/k7_compact_v1/root_{r:03d}.npz") as p:
        e=np.asarray(p["edge_raw"],dtype=np.int64); MASS[r]=np.asarray(p["mass_raw"],dtype=np.int64)
    s,a,b=np.nonzero(e); EDGE.append((s,a,b,e[s,a,b],int(e.shape[1])))
    rows=[]
    for d in cert["descriptors"]:
        if d["kind"]=="k7" and int(d["root"])==r:
            rl=np.asarray(d["rule"],dtype=np.uint8)
            sm=np.triu(np.equal.outer(rl,rl)).astype(np.int64)
            rows.append(np.tensordot(e,sm,axes=([1,2],[0,1])))
    POOL[r]=np.asarray(rows,dtype=np.float64)
MT=MASS.sum(axis=0)
KEY={}
def table(base):
    P=[];I=[];M=[]
    for parts in compositions(9,len(base)):
        k=canonical(*quotient(base,parts)); i=KEY.get(k)
        if i is None: i=KEY[k]=index9(blow_up(base,parts))
        P.append(parts); I.append(i); M.append(multinomial(parts))
    return np.asarray(P,dtype=float),np.asarray(I),np.asarray(M,dtype=float)
def q9_of(tab,w):
    P,I,M=tab; v=M*np.exp(P@np.log(np.maximum(w,1e-300))); q=np.zeros(N9); np.add.at(q,I,v); return q
def dens(base,w):
    E=edges_of(base); pr=[w[u]*w[v] for u,v in E]
    return 2*sum(pr), 2*min(sum(p for (u,v),p in zip(E,pr) if ((t>>u)&1)==((t>>v)&1)) for t in range(1<<len(base)))
def leg7_pool(q):
    return (10.0/DEN)*(sum(float(POOL[r].dot(q).min()) for r in range(107))-(2/25)*float(MT.dot(q)))/10.0

def min_rule_localsearch(N, starts, rng):
    """Upper bound on min over ALL rules of sum_{a<=b} N_ab [c_a=c_b] (N upper-triangular)."""
    k=N.shape[0]
    sym=N+N.T
    np.fill_diagonal(sym, 0.0)
    diag=float(np.trace(N))
    best=None
    for t in range(starts):
        c=np.zeros(k,dtype=np.int8) if t==0 else (np.asarray([rng.random()<0.5 for _ in range(k)],dtype=np.int8))
        while True:
            same=(c[:,None]==c[None,:])
            # gain of flipping vertex i = (weight to same side) - (weight to other side)
            g=(sym*same).sum(axis=1)-(sym*(~same)).sum(axis=1)
            i=int(np.argmax(g))
            if g[i]<=1e-12: break
            c[i]^=1
        val=diag+float((np.triu(N,1)*(c[:,None]==c[None,:])).sum())
        if best is None or val<best: best=val
    return best

def leg7_heur(q, starts=40, seed=1):
    rng=random.Random(seed); env=0.0
    for r in range(107):
        s,a,b,v,k=EDGE[r]
        w=v*q[s]
        keep=w!=0
        if not keep.any(): continue
        aa,bb,ww=a[keep],b[keep],w[keep]
        act=sorted(set(aa.tolist())|set(bb.tolist()))
        rm={p:i for i,p in enumerate(act)}
        n=len(act); N=np.zeros((n,n))
        np.add.at(N,( [rm[x] for x in aa.tolist()],[rm[x] for x in bb.tolist()]),ww)
        env+=min_rule_localsearch(N,starts,rng)
    return (10.0/DEN)*(env-(2/25)*float(MT.dot(q)))/10.0

GRO=grotzsch(); assert not triangles(GRO)
tab=table(GRO); h=11
rng=random.Random(11); best=(-9.9,None)
for att in range(8):
    w=np.full(h,1.0/h) if att==0 else np.asarray([rng.random()+1e-3 for _ in range(h)]); w/=w.sum()
    de,_=dens(GRO,w); cur=leg7_pool(q9_of(tab,w)) if LO<=de<=HI else -9.9
    sc=0.3
    for s in range(400):
        i,j=rng.randrange(h),rng.randrange(h)
        if i==j: continue
        t=w.copy(); d=sc*rng.random()*t[i]; t[i]-=d; t[j]+=d
        de,_=dens(GRO,t)
        if not LO<=de<=HI: continue
        v=leg7_pool(q9_of(tab,t))
        if v>cur: w,cur=t,v
        if s%50==49: sc*=0.7
    if cur>best[0]: best=(cur,w.copy())
val,w=best; q=q9_of(tab,w); de,dm=dens(GRO,w)
print("Grotzsch blow-up, weights:",np.round(w,5).tolist())
print("d_edge=%.6f  d_mono=%.6f   (band [%.4f,%.4f], 2/25=0.08)"%(de,dm,LO,HI))
print("leg7 with the certificate's K7 pool : %+.6e"%val)
t0=time.time()
hv=leg7_heur(q)
print("leg7 with local-search over ALL rules: %+.6e   [%.0fs]"%(hv,time.time()-t0))
print("\n=> the true 7-root envelope at this point exceeds 2/25:" , hv>0)
PET=petersen(); tp=table(PET); wp=np.full(10,0.1); qp=q9_of(tp,wp)
print("\ncontrol, Petersen: pool %+.6e   all-rules %+.6e"%(leg7_pool(qp),leg7_heur(qp)))
