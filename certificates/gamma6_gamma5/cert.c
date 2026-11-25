/* cert.c -- certificate checker (exact integer arithmetic) for
 *
 *     gamma_n >= THETA :  every configuration of n distinct points in the plane spans
 *                         C(n,2) lines whose directions (mod pi) have some cyclic gap >= THETA.
 *
 * THETA is given in units of 1/5120 degree (an integer), so all thresholds are integers.
 *
 * Normalisation (similarity invariance; reflections allowed):
 *   p0=(0,0), p1=(1,0) is a longest pair  =>  every |pi - pj| <= 1;
 *   the other points p2..p_{n-1} are labelled so that x2 <= x3 <= ... <= x_{n-1};
 *   reflection in the x-axis gives y2 >= 0; reflection x -> 1-x (swap p0,p1, relabel) gives x2 + x_{n-1} <= 1.
 *   Hence (pk) lies in [0,1] x [-1,1].
 * Coordinates are integers in units of 2^-K (exact dyadic boxes).  Boxes are closed.
 *
 * Suppose, for contradiction, that all cyclic gaps are < THETA  (units: N = 921600 = pi).
 * With L = C(n,2) lines and E = L*THETA - N > 0 the following hold (proved in the report):
 *   (H) for two lines p,q the ccw arc l(p->q) from d_p to d_q satisfies
 *          m*THETA - E < l < m*THETA,  where m = 1 + #{directions strictly inside the arc}  (1 <= m <= L-1).
 * Pruning tests (each a consequence of (H) for every configuration in the box):
 *   D : box violates the normalisation (exact integer tests).
 *   P : two lines whose direction intervals lie in a common arc of length <= THETA - E.
 *   J : for disjoint intervals p -> q, no admissible m in [1+jmin, 1+jmax] is compatible with [lmin,lmax].
 *   M : relative to a line p, the rank sets R_q = {m : window m meets [lmin,lmax]} (all ranks if unknown)
 *       admit no perfect matching lines <-> ranks 1..L-1 (Glover's greedy, exact for interval sets).
 * Direction intervals: the directions of the Minkowski difference box B_j - B_i (if it avoids 0) form the arc
 * spanned by its 4 corners; each corner angle is bracketed between consecutive entries of a table of integer
 * vectors v_m ~ 2^62 (cos(m pi/N), sin(m pi/N)) (angle error < 2^-61 rad), by exact __int128 cross products.
 * Every bracket is widened by one table unit on each side.  So all quantities are integers and all tests exact.
 *
 * usage: cert n THETA_units K nproc rank [d0]
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdint.h>
#include <time.h>

typedef __int128 i128;
typedef long long ll;
#define MAXN 7
#define MAXL 21
#define NT 921600LL               /* table units per pi; 1 unit = 1/5120 degree */

static ll *TC, *TS;              /* table */
static int n, L, D;
static ll TH, E;                 /* THETA, E in units */
static int K; static ll ONE;     /* coordinate unit 2^-K; ONE = 2^K */
static int PA[MAXL], PB[MAXL], LIDX[MAXN][MAXN];
static ll cnt_nodes, cnt_D, cnt_P, cnt_J, cnt_M, cnt_unres, maxdepth;

/* ---- angle bracket of integer vector (x,y), y >= 0, (x,y) != 0: returns m with
        angle(v_m) <= angle(x,y) <= angle(v_{m+1}),  0 <= m <= NT-1 ---- */
static inline i128 crs(ll m, ll x, ll y) { return (i128)TC[m] * y - (i128)TS[m] * x; }
static ll bracket(ll x, ll y) {
    double a = atan2((double)y, (double)x); /* guess only */
    ll m = (ll)floor(a / M_PI * NT);
    if (m < 0) m = 0; if (m > NT - 1) m = NT - 1;
    while (m > 0 && crs(m, x, y) < 0) m--;
    while (m < NT - 1 && crs(m + 1, x, y) > 0) m++;
    if (crs(m, x, y) < 0 || crs(m + 1, x, y) > 0) { fprintf(stderr, "bracket failure\n"); exit(2); }
    return m;
}

/* direction interval of difference box [ax0,ax1]x[ay0,ay1]: returns 0 if unknown (box contains 0
   or width >= NT); else lo in [0,NT), w >= 0 such that every direction (mod pi, units) lies in [lo, lo+w]. */
static int dirint(ll ax0, ll ax1, ll ay0, ll ay1, ll *lo, ll *w) {
    if (ax0 <= 0 && ax1 >= 0 && ay0 <= 0 && ay1 >= 0) return 0;
    ll cx[4] = {ax0, ax1, ax0, ax1}, cy[4] = {ay0, ay0, ay1, ay1};
    ll mn = 1LL << 60, mx = -(1LL << 60);
    if (ay0 > 0 || ay1 < 0) {                  /* open upper / lower half plane */
        int s = (ay0 > 0) ? 1 : -1;
        for (int c = 0; c < 4; c++) { ll m = bracket(s * cx[c], s * cy[c]);
            if (m < mn) mn = m; if (m + 1 > mx) mx = m + 1; }
    } else {                                    /* open right / left half plane, y range contains 0 */
        int s = (ax0 > 0) ? 1 : -1;
        for (int c = 0; c < 4; c++) { ll x = s * cx[c], y = s * cy[c];
            if (y >= 0) { ll m = bracket(x, y); if (m < mn) mn = m; if (m + 1 > mx) mx = m + 1; }
            else { ll m = bracket(x, -y); if (-(m + 1) < mn) mn = -(m + 1); if (-m > mx) mx = -m; } }
    }
    mn -= 1; mx += 1;                           /* absorb table error (< 2^-61 rad << 1 unit) */
    if (mx - mn >= NT) return 0;
    ll l = mn % NT; if (l < 0) l += NT;
    *lo = l; *w = mx - mn; return 1;
}

/* box: X[4k..4k+3] = x_lo, x_hi, y_lo, y_hi of point k+2 */
static void getpts(const ll *X, ll *x0, ll *x1, ll *y0, ll *y1) {
    x0[0] = x1[0] = y0[0] = y1[0] = 0; x0[1] = x1[1] = ONE; y0[1] = y1[1] = 0;
    for (int k = 2; k < n; k++) { x0[k] = X[4*(k-2)]; x1[k] = X[4*(k-2)+1]; y0[k] = X[4*(k-2)+2]; y1[k] = X[4*(k-2)+3]; }
}

static int test_domain(const ll *X) {
    if (X[3] < 0) return 1;                                         /* y2 >= 0 */
    for (int k = 2; k < n - 1; k++) if (X[4*(k-2)] > X[4*(k-1)+1]) return 1;   /* x sorted */
    if (X[0] + X[4*(n-3)] > ONE) return 1;                          /* x2 + x_{n-1} <= 1 */
    ll x0[MAXN], x1[MAXN], y0[MAXN], y1[MAXN]; getpts(X, x0, x1, y0, y1);
    i128 one2 = (i128)ONE * ONE;
    for (int i = 0; i < n; i++) for (int j = i + 1; j < n; j++) {
        ll a0 = x0[j]-x1[i], a1 = x1[j]-x0[i], b0 = y0[j]-y1[i], b1 = y1[j]-y0[i];
        ll mx = (a0 > 0) ? a0 : (a1 < 0 ? -a1 : 0), my = (b0 > 0) ? b0 : (b1 < 0 ? -b1 : 0);
        if ((i128)mx*mx + (i128)my*my > one2) return 1;             /* |pi-pj| <= 1 */
    }
    return 0;
}

static ll gLo[MAXL], gW[MAXL]; static int gKn[MAXL];

/* returns 0 = alive, 2 = P, 3 = J, 4 = M */
static int test_dirs(const ll *X) {
    ll x0[MAXN], x1[MAXN], y0[MAXN], y1[MAXN]; getpts(X, x0, x1, y0, y1);
    int nunk = 0;
    for (int l = 0; l < L; l++) { int i = PA[l], j = PB[l];
        gKn[l] = dirint(x0[j]-x1[i], x1[j]-x0[i], y0[j]-y1[i], y1[j]-y0[i], &gLo[l], &gW[l]);
        if (!gKn[l]) nunk++; }
    ll *lo = gLo, *w = gW; int *kn = gKn;
    ll rel[MAXL];
    for (int p = 0; p < L; p++) if (kn[p]) {
        for (int r = 0; r < L; r++) if (kn[r]) { ll t = lo[r] - lo[p]; if (t < 0) t += NT; rel[r] = t; }
        int ra[MAXL], rb[MAXL];
        for (int q = 0; q < L; q++) { ra[q] = 1; rb[q] = L - 1; }
        for (int q = 0; q < L; q++) if (q != p && kn[q]) {
            ll s = rel[q];
            if (!(w[p] < s && s + w[q] < NT)) {        /* not certainly ccw-after p */
                if (p < q) {
                    ll h1 = (w[p] > s + w[q]) ? w[p] : s + w[q];               /* arc starting at lo_p */
                    ll h2 = (w[q] > NT - s + w[p]) ? w[q] : NT - s + w[p];     /* arc starting at lo_q */
                    ll h = h1 < h2 ? h1 : h2;
                    if (h <= TH - E) return 2;
                }
                continue;
            }
            ll lmin = s - w[p], lmax = s + w[q];
            int jmin = 0, jmax = nunk;
            for (int r = 0; r < L; r++) if (r != p && r != q && kn[r]) {
                ll t = rel[r];
                if (t > w[p] && t + w[r] < s) jmin++;
                if (t < s + w[q] || t + w[r] > NT) jmax++;
            }
            int ok = 0;
            for (int jj = jmin; jj <= jmax; jj++) if (lmin < (jj+1)*TH && lmax > (jj+1)*TH - E) { ok = 1; break; }
            if (!ok) return 3;
            ll a = lmin / TH + 1, b = (lmax + E - 1) / TH;
            if (a < 1) a = 1; if (b > L - 1) b = L - 1;
            if (a > b) return 3;
            ra[q] = (int)a; rb[q] = (int)b;
        }
        /* Glover greedy: ranks 1..L-1, lines q != p */
        int used[MAXL]; for (int q = 0; q < L; q++) used[q] = (q == p);
        for (int m = 1; m < L; m++) { int bq = -1;
            for (int q = 0; q < L; q++) if (!used[q] && ra[q] <= m && (bq < 0 || rb[q] < rb[bq])) bq = q;
            if (bq < 0 || rb[bq] < m) return 4;
            used[bq] = 1; }
    }
    return 0;
}

static ll bsz(const ll *X, int i) { if (i < 2) return -1; const ll *b = X + 4*(i-2);
    ll a = b[1]-b[0], c = b[3]-b[2]; return a > c ? a : c; }

static int choose_split(const ll *X) {   /* widest direction interval (unknown = NT); split its larger endpoint box */
    /* Lines with both endpoints fixed (only p0p1) have no box to split and are skipped (referee fix:
       previously line 0 could be chosen, giving bi = 1 and an out-of-bounds index). Skipping it does not
       change any choice in which some other line is strictly wider, which was the case in every certified run. */
    ll ww = -1; int bi = -1;
    for (int l = 0; l < L; l++) { int i = PA[l], j = PB[l];
        if (i < 2 && j < 2) continue;
        ll wl = gKn[l] ? gW[l] : NT;
        if (wl > ww) { ww = wl; bi = (bsz(X, i) > bsz(X, j)) ? i : j; } }
    if (bi < 2 || bi >= n) { fprintf(stderr, "choose_split: invalid point index %d\n", bi); exit(3); }
    int idx = 4*(bi-2); if (X[idx+3]-X[idx+2] > X[idx+1]-X[idx]) idx += 2;
    if (idx < 0 || idx + 1 >= 2*D) { fprintf(stderr, "choose_split: index out of range %d\n", idx); exit(3); }
    return idx;
}

static FILE *unresf, *samplef; static ll sample_every = 0, sample_ctr = 0;
static void sample_leaf(const ll *X, int code) {   /* optional: dump every k-th pruned leaf for independent re-checking */
    if (!samplef || ++sample_ctr % sample_every) return;
    fprintf(samplef, "%d", code); for (int t = 0; t < 2*D; t++) fprintf(samplef, " %lld", X[t]); fprintf(samplef, "\n"); }
static void rec(ll *X, int depth) {
    cnt_nodes++; if (depth > maxdepth) maxdepth = depth;
    if (test_domain(X)) { cnt_D++; sample_leaf(X, 1); return; }
    int r = test_dirs(X);
    if (r) sample_leaf(X, r);
    if (r == 2) { cnt_P++; return; } if (r == 3) { cnt_J++; return; } if (r == 4) { cnt_M++; return; }
    int c = choose_split(X);
    if (X[c+1] - X[c] < 2 || depth > 200) { cnt_unres++;
        if (unresf) { for (int t = 0; t < 2*D; t++) fprintf(unresf, "%lld ", X[t]); fprintf(unresf, "\n"); fflush(unresf); }
        return; }
    ll Y[4*MAXN]; memcpy(Y, X, sizeof(ll)*2*D);
    ll m = X[c] + (X[c+1] - X[c]) / 2;
    Y[c+1] = m; rec(Y, depth + 1);
    Y[c+1] = X[c+1]; Y[c] = m; rec(Y, depth + 1);
}

/* soundness self-test: follow the branch containing a given grid point; report any pruning */
static int trace(ll *P, int maxd) {
    ll X[4*MAXN];
    for (int k = 0; k < n - 2; k++) { X[4*k] = 0; X[4*k+1] = ONE; X[4*k+2] = -ONE; X[4*k+3] = ONE; }
    for (int d = 0; d < maxd; d++) {
        for (int c = 0; c < D; c++) if (P[c] < X[2*c] || P[c] > X[2*c+1]) { printf("point left box?!\n"); return -1; }
        if (test_domain(X)) return 100 + d;
        int r = test_dirs(X); if (r) return 1000 * r + d;
        int c = choose_split(X); if (X[c+1] - X[c] < 2) return 0;
        ll m = X[c] + (X[c+1] - X[c]) / 2;
        if (P[c/2 * 1] <= m && 0) {}
        int coord = c / 2;  /* c = 2*coord */
        if (P[coord] <= m) X[c+1] = m; else X[c] = m;
    }
    return 0;
}

int main(int argc, char **argv) {
    if (argc >= 2 && !strcmp(argv[1], "trace")) {  /* cert trace n THETA K maxd < points */
        n = atoi(argv[2]); TH = atoll(argv[3]); K = atoi(argv[4]); int maxd = atoi(argv[5]);
        ONE = 1LL << K; L = n*(n-1)/2; D = 2*(n-2); E = L*TH - NT;
        { int l = 0; for (int i = 0; i < n; i++) for (int j = i+1; j < n; j++) { PA[l] = i; PB[l] = j; l++; } }
        char fn[256]; snprintf(fn, sizeof fn, "table_%lld.bin", NT); FILE *f = fopen(fn, "rb");
        TC = malloc(sizeof(ll)*(NT+1)); TS = malloc(sizeof(ll)*(NT+1));
        for (ll m = 0; m <= NT; m++) { ll v[2]; if (fread(v, sizeof(ll), 2, f) != 2) return 1; TC[m] = v[0]; TS[m] = v[1]; }
        ll P[2*MAXN]; int cnt = 0, bad = 0;
        while (1) { for (int c = 0; c < D; c++) if (scanf("%lld", &P[c]) != 1) goto done;
            int r = trace(P, maxd); cnt++; if (r) { bad++; printf("PRUNED config %d code %d\n", cnt, r); } }
        done: printf("traced %d configs, pruned (unsound) %d\n", cnt, bad); return bad != 0;
    }
    if (argc < 6) { fprintf(stderr, "usage: cert n THETA_units K nproc rank [d0]\n"); return 1; }
    n = atoi(argv[1]); TH = atoll(argv[2]); K = atoi(argv[3]); int nproc = atoi(argv[4]), rank = atoi(argv[5]);
    int d0 = argc > 6 ? atoi(argv[6]) : 14;
    if (n < 4 || n > MAXN || K < 4 || K > 44) { fprintf(stderr, "bad n or K\n"); return 1; }
    ONE = 1LL << K; L = n*(n-1)/2; D = 2*(n-2); E = L*TH - NT;
    if (E <= 0 || TH - E <= 0) { fprintf(stderr, "need 0 < E < THETA\n"); return 1; }
    { int l = 0; for (int i = 0; i < n; i++) for (int j = i+1; j < n; j++) { PA[l] = i; PB[l] = j; LIDX[i][j] = LIDX[j][i] = l; l++; } }
    /* load table */
    char fn[256]; snprintf(fn, sizeof fn, "table_%lld.bin", NT);
    FILE *f = fopen(fn, "rb"); if (!f) { fprintf(stderr, "no table\n"); return 1; }
    TC = malloc(sizeof(ll)*(NT+1)); TS = malloc(sizeof(ll)*(NT+1));
    for (ll m = 0; m <= NT; m++) { ll v[2]; if (fread(v, sizeof(ll), 2, f) != 2) { fprintf(stderr, "short table\n"); return 1; } TC[m] = v[0]; TS[m] = v[1]; }
    fclose(f);
    /* sanity: table strictly ccw-increasing, endpoints exact */
    if (TC[0] != (1LL<<62) || TS[0] != 0 || TC[NT] != -(1LL<<62) || TS[NT] != 0) { fprintf(stderr, "table endpoints\n"); return 1; }
    for (ll m = 0; m < NT; m++) if ((i128)TC[m]*TS[m+1] - (i128)TS[m]*TC[m+1] <= 0) { fprintf(stderr, "table order\n"); return 1; }
    char un[256]; snprintf(un, sizeof un, "unres_n%d_t%lld_r%d.txt", n, TH, rank); unresf = fopen(un, "w");
    if (getenv("SAMPLE")) { sample_every = atoll(getenv("SAMPLE")); char sf[256];
        snprintf(sf, sizeof sf, "leaves_n%d_t%lld_r%d.txt", n, TH, rank); samplef = fopen(sf, "w"); }
    clock_t t0 = clock();
    /* breadth-first expansion to depth d0 (deterministic), then DFS on own share */
    int cap = 1 << 18; ll *Q = malloc(sizeof(ll) * 2 * D * cap), *Q2 = malloc(sizeof(ll) * 2 * D * cap);
    int nq = 1;
    for (int k = 0; k < n - 2; k++) { Q[4*k] = 0; Q[4*k+1] = ONE; Q[4*k+2] = -ONE; Q[4*k+3] = ONE; }
    ll pre_nodes = 0, preD = 0, preP = 0, preJ = 0, preM = 0;
    for (int d = 0; d < d0; d++) {
        int n2 = 0;
        for (int b = 0; b < nq; b++) {
            ll *X = Q + (ll)b * 2 * D; pre_nodes++;
            if (test_domain(X)) { preD++; continue; }
            int r = test_dirs(X);
            if (r == 2) { preP++; continue; } if (r == 3) { preJ++; continue; } if (r == 4) { preM++; continue; }
            int c = choose_split(X); ll m = X[c] + (X[c+1]-X[c]) / 2;
            if (n2 + 2 > cap) { fprintf(stderr, "queue cap\n"); return 1; }
            ll *A = Q2 + (ll)n2 * 2 * D; memcpy(A, X, sizeof(ll)*2*D); A[c+1] = m; n2++;
            ll *B = Q2 + (ll)n2 * 2 * D; memcpy(B, X, sizeof(ll)*2*D); B[c] = m; n2++;
        }
        ll *t = Q; Q = Q2; Q2 = t; nq = n2;
    }
    printf("n=%d THETA=%lld units (=%.10f deg) E=%lld K=%d nproc=%d rank=%d d0=%d frontier=%d\n",
           n, TH, TH/5120.0, E, K, nproc, rank, d0, nq);
    if (rank == 0) printf("pre-phase (counted by rank 0 only): nodes=%lld D=%lld P=%lld J=%lld M=%lld\n", pre_nodes, preD, preP, preJ, preM);
    for (int b = rank; b < nq; b += nproc) rec(Q + (ll)b * 2 * D, d0);
    double secs = (double)(clock() - t0) / CLOCKS_PER_SEC;
    printf("DFS: boxes=%d..%d step %d nodes=%lld D=%lld P=%lld J=%lld M=%lld UNRESOLVED=%lld maxdepth=%lld cpu=%.1fs\n",
           rank, nq, nproc, cnt_nodes, cnt_D, cnt_P, cnt_J, cnt_M, cnt_unres, maxdepth, secs);
    return cnt_unres ? 3 : 0;
}
