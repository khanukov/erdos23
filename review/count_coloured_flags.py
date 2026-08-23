"""Size of the two-coloured triangle-free flag algebra, order by order.

A max-cut proof carries the cut as part of the structure, so the LP/SDP
variables are triangle-free graphs with their vertices 2-coloured, up to
isomorphism *and* colour swap.  Burnside over Aut(G) x Z2 counts them without
enumerating colourings.
"""
import itertools, pickle, sys, time

sys.setrecursionlimit(10000)


def degrees(a):
    return [bin(m).count("1") for m in a]


def automorphisms(adj):
    n = len(adj); d = degrees(adj); out = []; img = [-1] * n

    def go(k, used):
        if k == n:
            out.append(tuple(img)); return
        for w in range(n):
            if (used >> w) & 1 or d[w] != d[k]:
                continue
            if all(((adj[k] >> j) & 1) == ((adj[w] >> img[j]) & 1) for j in range(k)):
                img[k] = w; go(k + 1, used | (1 << w)); img[k] = -1
    go(0, 0)
    return out


def cycle_lengths(perm):
    n = len(perm); seen = [False] * n; out = []
    for v in range(n):
        if seen[v]:
            continue
        L = 0; u = v
        while not seen[u]:
            seen[u] = True; u = perm[u]; L += 1
        out.append(L)
    return out


def coloured_count(adj):
    aut = automorphisms(adj)
    plain = 0; swapped = 0
    for p in aut:
        cyc = cycle_lengths(p)
        plain += 2 ** len(cyc)
        if all(L % 2 == 0 for L in cyc):
            swapped += 2 ** len(cyc)
    return (plain + swapped) // (2 * len(aut)) if (plain + swapped) % (2 * len(aut)) == 0 \
        else (plain + swapped) / (2 * len(aut))


def canonical_key(adj):
    n = len(adj)
    d = degrees(adj)
    col = [(d[v], tuple(sorted(d[w] for w in range(n) if (adj[v] >> w) & 1))) for v in range(n)]
    order = sorted(range(n), key=lambda v: col[v])
    cells, s = [], 0
    for i in range(1, n + 1):
        if i == n or col[order[i]] != col[order[s]]:
            cells.append(order[s:i]); s = i
    best = None
    for ch in itertools.product(*[itertools.permutations(c) for c in cells]):
        perm = [v for c in ch for v in c]
        key = tuple(1 if (adj[perm[x]] >> perm[y]) & 1 else 0
                    for x in range(n) for y in range(x + 1, n))
        if best is None or key < best:
            best = key
    return best


def extend(graphs):
    """All triangle-free graphs on n+1 vertices, from those on n."""
    out = {}
    for adj in graphs:
        n = len(adj)
        for mask in range(1 << n):
            ok = True
            for i in range(n):
                if not (mask >> i) & 1:
                    continue
                for j in range(i + 1, n):
                    if (mask >> j) & 1 and (adj[i] >> j) & 1:
                        ok = False; break
                if not ok:
                    break
            if not ok:
                continue
            new = [adj[i] | (((mask >> i) & 1) << n) for i in range(n)] + [mask]
            out.setdefault(canonical_key(new), new)
    return list(out.values())


graphs = [[0]]
print(f"{'n':>3} {'triangle-free':>14} {'2-coloured':>12}   {'growth':>7}")
prev = None
for n in range(1, 10):
    if n > 1:
        graphs = extend(graphs)
    started = time.time()
    coloured = sum(coloured_count(a) for a in graphs)
    growth = "" if prev is None else f"{coloured / prev:6.1f}x"
    print(f"{n:>3} {len(graphs):>14} {coloured:>12}   {growth:>7}   [{time.time() - started:.0f}s]",
          flush=True)
    prev = coloured
