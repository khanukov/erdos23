/*
 * blocks10.c -- order-10 flag moment blocks of the triangle-free flag algebra.
 *
 * For a type sigma of size s (a labelled triangle-free graph on s vertices, one
 * representative per isomorphism class) and f = (10 - s) / 2 free vertices,
 * the moment block of an order-10 state G is the integer matrix
 *
 *     M_sigma(G)[F1, F2] = #{ (t, A) : t an ordered s-tuple of V(G) inducing
 *                             exactly the labelled graph sigma, A a subset of
 *                             V(G) \ t of size f, (t, A) forms the sigma-flag F1
 *                             and (t, V(G) \ t \ A) forms the sigma-flag F2 }.
 *
 * For a graphon W with state vector q (q_G = probability that 10 random points
 * induce G), sum_G q_G M_sigma(G) / (10!/(10-s)! C(10-s, f)) is the matrix
 * E_x[ 1[x induces sigma] phi(x) phi(x)^T ] with phi_F(x) the probability that
 * x plus f random points form F.  It is positive semidefinite, so for every
 * vector v the row  r_G = v^T M_sigma(G) v  satisfies  <r, q> >= 0.
 *
 * Modes:
 *   blocks10 index  states.bin tuples.bin           (once) cache the kept tuples
 *   blocks10 moment states.bin tuples.bin q.bin out.bin
 *       out: for each level (s = 0,2,4,6 as requested) and type: int32 nF, then
 *            nF*nF float64 entries of sum_G q_G M_sigma(G) / normaliser
 *   blocks10 rows   states.bin tuples.bin vecs.bin out.bin
 *       vecs: int32 count; per vector: int32 s, int32 type, int32 nF, int64[nF]
 *       out : per vector int64[n_states] exact numerators v^T M_sigma(G) v
 *   blocks10 info   states.bin                      print the levels and sizes
 *
 * states.bin: int32 n_states, then n_states*10 uint16 adjacency bitmasks.
 * The optional final argument "levels" (default "0246") selects the type sizes.
 */
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define NV 10
#define MAXF 5

static int n_states;
static uint16_t *adj;

/* ---------- permutations ---------- */
static int n_perms[8];
static int *perms[8];

static void gen_perms(int k) {
    int n = 1;
    for (int i = 2; i <= k; i++) n *= i;
    perms[k] = malloc(sizeof(int) * (size_t)n * (k ? k : 1));
    n_perms[k] = n;
    int a[8], c[8] = {0};
    for (int i = 0; i < k; i++) a[i] = i;
    int count = 0;
    memcpy(perms[k] + (size_t)count * k, a, sizeof(int) * k);
    count++;
    int i = 0;
    while (i < k) {
        if (c[i] < i) {
            if (i % 2 == 0) { int t = a[0]; a[0] = a[i]; a[i] = t; }
            else { int t = a[c[i]]; a[c[i]] = a[i]; a[i] = t; }
            memcpy(perms[k] + (size_t)count * k, a, sizeof(int) * k);
            count++;
            c[i]++;
            i = 0;
        } else {
            c[i] = 0;
            i++;
        }
    }
    if (count != n) { fprintf(stderr, "permutation count mismatch\n"); exit(1); }
}

static inline int pidx(int i, int j) {
    if (i > j) { int t = i; i = j; j = t; }
    return j * (j - 1) / 2 + i;
}

/* pattern on k vertices: bit pidx(i,j) = edge ij */
static int relabel(int p, int k, const int *pi) {
    int out = 0;
    for (int i = 0; i < k; i++)
        for (int j = i + 1; j < k; j++)
            if (p >> pidx(i, j) & 1) out |= 1 << pidx(pi[i], pi[j]);
    return out;
}

static int has_triangle(int p, int k) {
    for (int i = 0; i < k; i++)
        for (int j = i + 1; j < k; j++)
            for (int l = j + 1; l < k; l++)
                if ((p >> pidx(i, j) & 1) && (p >> pidx(i, l) & 1) && (p >> pidx(j, l) & 1))
                    return 1;
    return 0;
}

/* ---------- levels ---------- */
typedef struct {
    int s, f;
    int n_types;
    int *type_of;      /* 2^C(s,2): representative pattern -> type, else -1 */
    int *rep;          /* type -> pattern */
    int *aut;          /* type -> |Aut| */
    int flag_bits;     /* f*s + C(f,2) */
    int *n_flags;      /* per type */
    int **flag_of;     /* per type: pattern -> flag index of its canonical form, -1 if invalid */
    int n_splits;      /* C(10-s, f) */
    int *split_a;      /* n_splits * f  indices into rem[] */
    int *split_b;      /* n_splits * f */
    double normaliser; /* 10!/(10-s)! * C(10-s,f) */
} Level;

static int n_levels;
static Level levels[4];

static int binom(int n, int k) {
    if (k < 0 || k > n) return 0;
    long r = 1;
    for (int i = 1; i <= k; i++) r = r * (n - k + i) / i;
    return (int)r;
}

/* flag pattern: nb_i at bits [i*s, (i+1)*s), free-free edge (i,j) at f*s + pidx(i,j) */
static int flag_relabel(int p, int s, int f, const int *pi) {
    int out = 0;
    for (int i = 0; i < f; i++) {
        int nb = (p >> (i * s)) & ((1 << s) - 1);
        out |= nb << (pi[i] * s);
    }
    for (int i = 0; i < f; i++)
        for (int j = i + 1; j < f; j++)
            if (p >> (f * s + pidx(i, j)) & 1) out |= 1 << (f * s + pidx(pi[i], pi[j]));
    return out;
}

static int flag_valid(int p, int s, int f, int sigma) {
    int nb[MAXF];
    for (int i = 0; i < f; i++) nb[i] = (p >> (i * s)) & ((1 << s) - 1);
    /* neighbourhoods independent in sigma */
    for (int i = 0; i < f; i++)
        for (int a = 0; a < s; a++)
            for (int b = a + 1; b < s; b++)
                if ((nb[i] >> a & 1) && (nb[i] >> b & 1) && (sigma >> pidx(a, b) & 1)) return 0;
    /* adjacent free vertices share no root neighbour; free graph triangle-free */
    int ff = p >> (f * s);
    if (has_triangle(ff, f)) return 0;
    for (int i = 0; i < f; i++)
        for (int j = i + 1; j < f; j++)
            if ((ff >> pidx(i, j) & 1) && (nb[i] & nb[j])) return 0;
    return 1;
}

static void build_level(Level *L, int s) {
    L->s = s;
    L->f = (NV - s) / 2;
    int f = L->f;
    int tb = s * (s - 1) / 2;
    int ntp = 1 << tb;
    L->type_of = malloc(sizeof(int) * ntp);
    L->n_types = 0;
    int *reps = malloc(sizeof(int) * ntp), *auts = malloc(sizeof(int) * ntp);
    for (int p = 0; p < ntp; p++) {
        L->type_of[p] = -1;
        if (has_triangle(p, s)) continue;
        int canon = p, aut = 0;
        for (int k = 0; k < n_perms[s]; k++) {
            int r = relabel(p, s, perms[s] + (size_t)k * s);
            if (r < canon) canon = r;
            if (r == p) aut++;
        }
        if (canon == p) {
            reps[L->n_types] = p;
            auts[L->n_types] = aut;
            L->type_of[p] = L->n_types++;
        }
    }
    L->rep = reps;
    L->aut = auts;
    L->flag_bits = f * s + f * (f - 1) / 2;
    int nfp = 1 << L->flag_bits;
    L->n_flags = malloc(sizeof(int) * L->n_types);
    L->flag_of = malloc(sizeof(int *) * L->n_types);
    for (int t = 0; t < L->n_types; t++) {
        int sigma = reps[t];
        int *tab = malloc(sizeof(int) * nfp);
        int *rep_index = malloc(sizeof(int) * nfp);
        int n = 0;
        for (int p = 0; p < nfp; p++) {
            rep_index[p] = -1;
            if (!flag_valid(p, s, f, sigma)) continue;
            int canon = p;
            for (int k = 0; k < n_perms[f]; k++) {
                int r = flag_relabel(p, s, f, perms[f] + (size_t)k * f);
                if (r < canon) canon = r;
            }
            if (canon == p) rep_index[p] = n++;
        }
        for (int p = 0; p < nfp; p++) {
            tab[p] = -1;
            if (!flag_valid(p, s, f, sigma)) continue;
            int canon = p;
            for (int k = 0; k < n_perms[f]; k++) {
                int r = flag_relabel(p, s, f, perms[f] + (size_t)k * f);
                if (r < canon) canon = r;
            }
            tab[p] = rep_index[canon];
            if (tab[p] < 0) { fprintf(stderr, "canonical flag missing\n"); exit(1); }
        }
        free(rep_index);
        L->flag_of[t] = tab;
        L->n_flags[t] = n;
    }
    /* splits of the 10-s remaining vertices into A (size f) and B (the rest) */
    int m = NV - s;
    L->n_splits = binom(m, f);
    L->split_a = malloc(sizeof(int) * L->n_splits * f);
    L->split_b = malloc(sizeof(int) * L->n_splits * f);
    int c = 0;
    for (int mask = 0; mask < (1 << m); mask++) {
        if (__builtin_popcount(mask) != f) continue;
        int ia = 0, ib = 0;
        for (int i = 0; i < m; i++) {
            if (mask >> i & 1) L->split_a[c * f + ia++] = i;
            else L->split_b[c * f + ib++] = i;
        }
        c++;
    }
    if (c != L->n_splits) { fprintf(stderr, "split count mismatch\n"); exit(1); }
    double norm = 1.0;
    for (int i = 0; i < s; i++) norm *= (NV - i);
    L->normaliser = norm * L->n_splits;
}

/* ---------- kept tuples ---------- */
/* entry: type << 24 | roots packed 4 bits each (position k at bits 4k) */
static uint32_t **kept;        /* per level: entries */
static uint64_t **kept_off;    /* per level: n_states + 1 offsets */

static void enum_roots(const uint16_t *a, Level *L, int k, int used, int pattern,
                       int *roots, uint32_t **buf, size_t *n, size_t *cap) {
    if (k == L->s) {
        int t = L->type_of[pattern];
        if (t < 0) return;
        uint32_t packed = (uint32_t)t << 24;
        for (int i = 0; i < L->s; i++) packed |= (uint32_t)roots[i] << (4 * i);
        if (*n == *cap) { *cap = *cap ? *cap * 2 : 1024; *buf = realloc(*buf, sizeof(uint32_t) * *cap); }
        (*buf)[(*n)++] = packed;
        return;
    }
    for (int v = 0; v < NV; v++) {
        if (used >> v & 1) continue;
        int p = pattern;
        for (int i = 0; i < k; i++)
            if (a[roots[i]] >> v & 1) p |= 1 << pidx(i, k);
        roots[k] = v;
        enum_roots(a, L, k + 1, used | (1 << v), p, roots, buf, n, cap);
    }
}

static void build_index(const char *path) {
    kept = malloc(sizeof(uint32_t *) * n_levels);
    kept_off = malloc(sizeof(uint64_t *) * n_levels);
    for (int l = 0; l < n_levels; l++) {
        uint32_t *buf = NULL;
        size_t n = 0, cap = 0;
        kept_off[l] = malloc(sizeof(uint64_t) * (n_states + 1));
        for (int st = 0; st < n_states; st++) {
            kept_off[l][st] = n;
            int roots[8];
            enum_roots(adj + (size_t)st * NV, &levels[l], 0, 0, 0, roots, &buf, &n, &cap);
        }
        kept_off[l][n_states] = n;
        kept[l] = buf;
        fprintf(stderr, "level s=%d: %zu kept tuples (%.1f per state)\n",
                levels[l].s, n, (double)n / n_states);
    }
    FILE *fp = fopen(path, "wb");
    if (!fp) { perror(path); exit(1); }
    fwrite(&n_states, sizeof(int), 1, fp);
    fwrite(&n_levels, sizeof(int), 1, fp);
    for (int l = 0; l < n_levels; l++) {
        fwrite(&levels[l].s, sizeof(int), 1, fp);
        fwrite(kept_off[l], sizeof(uint64_t), n_states + 1, fp);
        fwrite(kept[l], sizeof(uint32_t), kept_off[l][n_states], fp);
    }
    fclose(fp);
}

static void load_index(const char *path) {
    FILE *fp = fopen(path, "rb");
    if (!fp) { perror(path); exit(1); }
    int ns, nl;
    if (fread(&ns, sizeof(int), 1, fp) != 1 || fread(&nl, sizeof(int), 1, fp) != 1 ||
        ns != n_states || nl != n_levels) {
        fprintf(stderr, "index does not match states/levels\n"); exit(1);
    }
    kept = malloc(sizeof(uint32_t *) * n_levels);
    kept_off = malloc(sizeof(uint64_t *) * n_levels);
    for (int l = 0; l < n_levels; l++) {
        int s;
        if (fread(&s, sizeof(int), 1, fp) != 1 || s != levels[l].s) {
            fprintf(stderr, "index level mismatch\n"); exit(1);
        }
        kept_off[l] = malloc(sizeof(uint64_t) * (n_states + 1));
        if (fread(kept_off[l], sizeof(uint64_t), n_states + 1, fp) != (size_t)n_states + 1) {
            fprintf(stderr, "short index\n"); exit(1);
        }
        size_t n = kept_off[l][n_states];
        kept[l] = malloc(sizeof(uint32_t) * (n ? n : 1));
        if (fread(kept[l], sizeof(uint32_t), n, fp) != n) { fprintf(stderr, "short index\n"); exit(1); }
    }
    fclose(fp);
}

/* ---------- passes ---------- */
typedef void (*visit_fn)(int level, int type, int f1, int f2, void *ctx);

static void walk_state(int st, visit_fn visit, void *ctx) {
    const uint16_t *a = adj + (size_t)st * NV;
    for (int l = 0; l < n_levels; l++) {
        Level *L = &levels[l];
        int s = L->s, f = L->f, m = NV - s;
        for (uint64_t e = kept_off[l][st]; e < kept_off[l][st + 1]; e++) {
            uint32_t packed = kept[l][e];
            int t = packed >> 24;
            int roots[8], inroot = 0;
            for (int i = 0; i < s; i++) { roots[i] = (packed >> (4 * i)) & 15; inroot |= 1 << roots[i]; }
            int rem[NV], nb[NV], r = 0;
            for (int v = 0; v < NV; v++) {
                if (inroot >> v & 1) continue;
                int b = 0;
                for (int i = 0; i < s; i++) if (a[v] >> roots[i] & 1) b |= 1 << i;
                nb[r] = b;
                rem[r++] = v;
            }
            (void)m;
            const int *tab = L->flag_of[t];
            for (int c = 0; c < L->n_splits; c++) {
                const int *A = L->split_a + c * f, *B = L->split_b + c * f;
                int pa = 0, pb = 0;
                for (int i = 0; i < f; i++) { pa |= nb[A[i]] << (i * s); pb |= nb[B[i]] << (i * s); }
                for (int i = 0; i < f; i++)
                    for (int j = i + 1; j < f; j++) {
                        if (a[rem[A[i]]] >> rem[A[j]] & 1) pa |= 1 << (f * s + pidx(i, j));
                        if (a[rem[B[i]]] >> rem[B[j]] & 1) pb |= 1 << (f * s + pidx(i, j));
                    }
                int f1 = tab[pa], f2 = tab[pb];
                if (f1 < 0 || f2 < 0) { fprintf(stderr, "invalid flag in a triangle-free state\n"); exit(1); }
                visit(l, t, f1, f2, ctx);
            }
        }
    }
}

/* moment */
typedef struct { double w; double ***M; } MomentCtx;
static void moment_visit(int l, int t, int f1, int f2, void *ctx) {
    MomentCtx *c = ctx;
    c->M[l][t][(size_t)f1 * levels[l].n_flags[t] + f2] += c->w;
}

/* rows */
typedef struct {
    int n_vecs;
    int *vec_level, *vec_type;
    int64_t **vec;
    int64_t *acc;      /* n_vecs * n_states */
    int st;
    int **by_type_start; /* per level: per type, index into order[] */
    int *order;          /* vectors grouped by (level,type) */
} RowsCtx;
static void rows_visit(int l, int t, int f1, int f2, void *ctx) {
    RowsCtx *c = ctx;
    int start = c->by_type_start[l][t], end = c->by_type_start[l][t + 1];
    for (int k = start; k < end; k++) {
        int j = c->order[k];
        c->acc[(size_t)j * n_states + c->st] += c->vec[j][f1] * c->vec[j][f2];
    }
}

static void read_states(const char *path) {
    FILE *fp = fopen(path, "rb");
    if (!fp) { perror(path); exit(1); }
    if (fread(&n_states, sizeof(int), 1, fp) != 1) { fprintf(stderr, "bad states\n"); exit(1); }
    adj = malloc(sizeof(uint16_t) * (size_t)n_states * NV);
    if (fread(adj, sizeof(uint16_t), (size_t)n_states * NV, fp) != (size_t)n_states * NV) {
        fprintf(stderr, "short states\n"); exit(1);
    }
    fclose(fp);
    for (int st = 0; st < n_states; st++) {
        const uint16_t *a = adj + (size_t)st * NV;
        for (int u = 0; u < NV; u++)
            for (int v = u + 1; v < NV; v++) {
                if (((a[u] >> v) & 1) != ((a[v] >> u) & 1)) { fprintf(stderr, "asymmetric state %d\n", st); exit(1); }
                if ((a[u] >> v & 1) && (a[u] & a[v])) { fprintf(stderr, "state %d has a triangle\n", st); exit(1); }
            }
    }
}

int main(int argc, char **argv) {
    if (argc < 3) { fprintf(stderr, "usage: see header\n"); return 1; }
    const char *mode = argv[1];
    read_states(argv[2]);
    const char *level_spec = "0246";
    if (!strcmp(mode, "info") && argc > 3) level_spec = argv[3];
    if (!strcmp(mode, "index") && argc > 4) level_spec = argv[4];
    if ((!strcmp(mode, "moment") || !strcmp(mode, "rows")) && argc > 6) level_spec = argv[6];
    for (int k = 0; k <= 6; k++) gen_perms(k);
    n_levels = 0;
    for (const char *c = level_spec; *c; c++) build_level(&levels[n_levels++], *c - '0');
    if (!strcmp(mode, "info")) {
        for (int l = 0; l < n_levels; l++) {
            Level *L = &levels[l];
            printf("level s=%d f=%d: %d types, splits %d, normaliser %.0f\n", L->s, L->f, L->n_types,
                   L->n_splits, L->normaliser);
            for (int t = 0; t < L->n_types; t++)
                printf("  type %d pattern %d aut %d flags %d\n", t, L->rep[t], L->aut[t], L->n_flags[t]);
        }
        return 0;
    }
    if (!strcmp(mode, "index")) { build_index(argv[3]); return 0; }
    load_index(argv[3]);
    if (!strcmp(mode, "moment")) {
        FILE *fq = fopen(argv[4], "rb");
        if (!fq) { perror(argv[4]); return 1; }
        double *q = malloc(sizeof(double) * n_states);
        if (fread(q, sizeof(double), n_states, fq) != (size_t)n_states) { fprintf(stderr, "short q\n"); return 1; }
        fclose(fq);
        MomentCtx ctx;
        ctx.M = malloc(sizeof(double **) * n_levels);
        for (int l = 0; l < n_levels; l++) {
            ctx.M[l] = malloc(sizeof(double *) * levels[l].n_types);
            for (int t = 0; t < levels[l].n_types; t++)
                ctx.M[l][t] = calloc((size_t)levels[l].n_flags[t] * levels[l].n_flags[t], sizeof(double));
        }
        for (int st = 0; st < n_states; st++) {
            if (q[st] == 0.0) continue;
            ctx.w = q[st];
            walk_state(st, moment_visit, &ctx);
        }
        FILE *fo = fopen(argv[5], "wb");
        if (!fo) { perror(argv[5]); return 1; }
        for (int l = 0; l < n_levels; l++)
            for (int t = 0; t < levels[l].n_types; t++) {
                int nF = levels[l].n_flags[t];
                size_t nn = (size_t)nF * nF;
                for (size_t i = 0; i < nn; i++) ctx.M[l][t][i] /= levels[l].normaliser;
                fwrite(&nF, sizeof(int), 1, fo);
                fwrite(ctx.M[l][t], sizeof(double), nn, fo);
            }
        fclose(fo);
        return 0;
    }
    if (!strcmp(mode, "rows")) {
        FILE *fv = fopen(argv[4], "rb");
        if (!fv) { perror(argv[4]); return 1; }
        RowsCtx ctx;
        if (fread(&ctx.n_vecs, sizeof(int), 1, fv) != 1) { fprintf(stderr, "bad vecs\n"); return 1; }
        ctx.vec_level = malloc(sizeof(int) * ctx.n_vecs);
        ctx.vec_type = malloc(sizeof(int) * ctx.n_vecs);
        ctx.vec = malloc(sizeof(int64_t *) * ctx.n_vecs);
        for (int j = 0; j < ctx.n_vecs; j++) {
            int s, t, nF;
            if (fread(&s, sizeof(int), 1, fv) != 1 || fread(&t, sizeof(int), 1, fv) != 1 ||
                fread(&nF, sizeof(int), 1, fv) != 1) { fprintf(stderr, "bad vecs\n"); return 1; }
            int l = -1;
            for (int k = 0; k < n_levels; k++) if (levels[k].s == s) l = k;
            if (l < 0 || t < 0 || t >= levels[l].n_types || nF != levels[l].n_flags[t]) {
                fprintf(stderr, "vector %d: bad level/type/size\n", j); return 1;
            }
            ctx.vec_level[j] = l;
            ctx.vec_type[j] = t;
            ctx.vec[j] = malloc(sizeof(int64_t) * nF);
            if (fread(ctx.vec[j], sizeof(int64_t), nF, fv) != (size_t)nF) { fprintf(stderr, "short vecs\n"); return 1; }
            for (int i = 0; i < nF; i++)
                if (ctx.vec[j][i] > 1000000 || ctx.vec[j][i] < -1000000) {
                    fprintf(stderr, "vector entries must be at most 1e6 in absolute value\n"); return 1;
                }
        }
        fclose(fv);
        ctx.by_type_start = malloc(sizeof(int *) * n_levels);
        ctx.order = malloc(sizeof(int) * (ctx.n_vecs ? ctx.n_vecs : 1));
        int pos = 0;
        for (int l = 0; l < n_levels; l++) {
            ctx.by_type_start[l] = malloc(sizeof(int) * (levels[l].n_types + 1));
            for (int t = 0; t < levels[l].n_types; t++) {
                ctx.by_type_start[l][t] = pos;
                for (int j = 0; j < ctx.n_vecs; j++)
                    if (ctx.vec_level[j] == l && ctx.vec_type[j] == t) ctx.order[pos++] = j;
            }
            ctx.by_type_start[l][levels[l].n_types] = pos;
        }
        ctx.acc = calloc((size_t)ctx.n_vecs * n_states, sizeof(int64_t));
        for (int st = 0; st < n_states; st++) {
            ctx.st = st;
            walk_state(st, rows_visit, &ctx);
        }
        FILE *fo = fopen(argv[5], "wb");
        if (!fo) { perror(argv[5]); return 1; }
        fwrite(ctx.acc, sizeof(int64_t), (size_t)ctx.n_vecs * n_states, fo);
        fclose(fo);
        return 0;
    }
    fprintf(stderr, "unknown mode %s\n", mode);
    return 1;
}
