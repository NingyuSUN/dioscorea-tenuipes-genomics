/* Read-only exact-allele VCF audit, HTSlib 1.13. Phasing ignored in GT
 * comparison. Diploid and haploid GT supported; unsupported ploidy,
 * malformed/sortedness errors and duplicate exact variants fail explicitly. No
 * input indices are created. Build: cc -O3 -I<htslib-1.13> this.c -lhts -o
 * compare_population_genotypes Run: compare_population_genotypes left.vcf.gz
 * right.vcf.gz output_prefix
 */
#include <htslib/vcf.h>
#include <inttypes.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#define U uint64_t
#define FAIL(msg)                                                              \
  do {                                                                         \
    fprintf(stderr, "ERROR: %s\n", msg);                                       \
    exit(2);                                                                   \
  } while (0)
typedef struct {
  U n, missing, partial, non_diploid, homref, het, homalt, dp_n, dp_sum, dp0ref,
      gq_n, gq_sum, good, changed, filled;
} Stats;
typedef struct {
  htsFile *f;
  bcf_hdr_t *h;
  bcf1_t *pending;
  int eof, last_rid;
  hts_pos_t last_pos;
  int ns;
  U n;
  Stats *s;
  int32_t *gt, *dp, *gq;
  int ngt, ndp, ngq;
} Reader;
typedef struct {
  bcf1_t **r;
  int n;
  int rid;
  hts_pos_t pos;
} Group;
static char *alias(const char *s) {
  char *p = strdup(s);
  size_t n = strlen(p);
  if (n > 9 && !strncmp(p, "Dte_", 4) && !strcmp(p + n - 5, "_L001")) {
    memmove(p, p + 4, n - 3);
    p[n - 9] = 0;
    for (char *q = p; *q; q++)
      if (*q == '_')
        *q = '-';
  }
  return p;
}
static int genotype(int32_t *g, int stride, int *a, int *miss, int *hr) {
  int n = 0;
  *miss = 0;
  *hr = 1;
  for (int j = 0; j < stride; j++) {
    if (g[j] == bcf_int32_vector_end)
      break;
    int v = bcf_gt_is_missing(g[j]) ? -1 : bcf_gt_allele(g[j]);
    if (n >= 2)
      FAIL("Ploidy above two is unsupported");
    a[n++] = v;
    if (v < 0)
      *miss = 1;
    if (v != 0)
      *hr = 0;
  }
  if (!n) {
    *miss = 1;
    *hr = 0;
  }
  if (n == 2 && a[0] > a[1]) {
    int t = a[0];
    a[0] = a[1];
    a[1] = t;
  }
  return n;
}
static Reader open_reader(const char *path) {
  Reader a = {0};
  a.f = hts_open(path, "r");
  if (!a.f)
    FAIL("Cannot open VCF");
  hts_set_threads(a.f, 2);
  a.h = bcf_hdr_read(a.f);
  if (!a.h)
    FAIL("Cannot read header");
  a.ns = bcf_hdr_nsamples(a.h);
  if (a.ns < 1)
    FAIL("No samples");
  a.s = calloc(a.ns, sizeof(Stats));
  a.last_rid = -1;
  a.last_pos = -1;
  return a;
}
static int fetch(Reader *a, bcf1_t *r) {
  int ret = bcf_read(a->f, a->h, r);
  if (ret == -1) {
    a->eof = 1;
    return 0;
  }
  if (ret < 0 || r->errcode)
    FAIL("Malformed/truncated record");
  if (r->rid < a->last_rid || (r->rid == a->last_rid && r->pos < a->last_pos))
    FAIL("Unsorted VCF");
  a->last_rid = r->rid;
  a->last_pos = r->pos;
  return 1;
}
static void consume(Reader *a, bcf1_t *r) {
  a->n++;
  int nt = bcf_get_genotypes(a->h, r, &a->gt, &a->ngt);
  if (nt <= 0 || nt % a->ns)
    FAIL("GT absent or invalid");
  int stride = nt / a->ns;
  int nd = bcf_get_format_int32(a->h, r, "DP", &a->dp, &a->ndp),
      nq = bcf_get_format_int32(a->h, r, "GQ", &a->gq, &a->ngq);
  if (nd > 0 && nd != a->ns)
    FAIL("DP is not scalar");
  if (nq > 0 && nq != a->ns)
    FAIL("GQ is not scalar");
  for (int i = 0; i < a->ns; i++) {
    Stats *s = &a->s[i];
    int g[2] = {-1, -1}, miss, hr;
    int n = genotype(a->gt + i * stride, stride, g, &miss, &hr);
    int dp = nd > 0 ? a->dp[i] : -1, gq = nq > 0 ? a->gq[i] : -1;
    s->n++;
    s->missing += miss;
    s->partial += miss && n == 2 && g[1] >= 0;
    s->non_diploid += n != 2;
    if (dp >= 0) {
      s->dp_n++;
      s->dp_sum += dp;
    }
    if (gq >= 0) {
      s->gq_n++;
      s->gq_sum += gq;
    }
    if (!miss) {
      s->homref += hr;
      s->het += n == 2 && g[0] != g[1];
      s->homalt += !hr && (n != 2 || g[0] == g[1]);
      s->dp0ref += hr && dp <= 0;
      s->good += dp >= 5 && gq >= 20;
    }
  }
}
static void free_group(Group *g) {
  for (int i = 0; i < g->n; i++)
    bcf_destroy(g->r[i]);
  free(g->r);
  memset(g, 0, sizeof(*g));
}
static Group next_group(Reader *a) {
  Group g = {0};
  bcf1_t *r = a->pending;
  a->pending = NULL;
  if (!r) {
    if (a->eof)
      return g;
    r = bcf_init();
    if (!fetch(a, r)) {
      bcf_destroy(r);
      return g;
    }
  }
  g.rid = r->rid;
  g.pos = r->pos;
  while (1) {
    consume(a, r);
    bcf_unpack(r, BCF_UN_STR);
    g.r = realloc(g.r, (g.n + 1) * sizeof(*g.r));
    g.r[g.n++] = r;
    r = bcf_init();
    if (!fetch(a, r)) {
      bcf_destroy(r);
      break;
    }
    if (r->rid != g.rid || r->pos != g.pos) {
      a->pending = r;
      break;
    }
  }
  return g;
}
static int alleles_equal(bcf1_t *x, bcf1_t *y) {
  if (x->n_allele != y->n_allele)
    return 0;
  for (int j = 0; j < x->n_allele; j++)
    if (strcmp(x->d.allele[j], y->d.allele[j]))
      return 0;
  return 1;
}
static void duplicates(Group *g) {
  for (int i = 0; i < g->n; i++)
    for (int j = 0; j < i; j++)
      if (alleles_equal(g->r[i], g->r[j]))
        FAIL("Duplicate exact variant");
}
static void save_stats(Reader *a, const char *pre, const char *label) {
  char path[4096];
  snprintf(path, sizeof path, "%s.%s.samples.tsv", pre, label);
  FILE *f = fopen(path, "wx");
  if (!f)
    FAIL("Output exists or cannot be created");
  fprintf(f, "sample\trecords\tmissing\tpartial_missing\tnon_"
             "diploid\thomref\thet\thomalt\tdp_observed\tdp_sum\thomref_dp_"
             "zero_or_missing\tgq_observed\tgq_sum\tcalled_dp5_gq20\tchanged_"
             "gt\tmissing_to_homref\n");
  for (int i = 0; i < a->ns; i++) {
    Stats s = a->s[i];
    char *al = alias(a->h->samples[i]);
    fprintf(f, "%s", al);
    free(al);
    U v[] = {s.n,    s.missing, s.partial, s.non_diploid, s.homref,
             s.het,  s.homalt,  s.dp_n,    s.dp_sum,      s.dp0ref,
             s.gq_n, s.gq_sum,  s.good,    s.changed,     s.filled};
    for (int j = 0; j < 15; j++)
      fprintf(f, "\t%" PRIu64, v[j]);
    fputc('\n', f);
  }
  if (fclose(f))
    FAIL("Output write failed");
}
int main(int argc, char **argv) {
  if (argc != 4)
    FAIL("Usage: compare left.vcf.gz right.vcf.gz output_prefix");
  Reader a = open_reader(argv[1]), b = open_reader(argv[2]);
  if (a.ns != b.ns || a.h->n[BCF_DT_CTG] != b.h->n[BCF_DT_CTG])
    FAIL("Different sample/contig counts");
  for (int i = 0; i < a.ns; i++) {
    char *x = alias(a.h->samples[i]), *y = alias(b.h->samples[i]);
    if (strcmp(x, y))
      FAIL("Sample order/alias mismatch");
    free(x);
    free(y);
  }
  for (int i = 0; i < a.h->n[BCF_DT_CTG]; i++)
    if (strcmp(bcf_hdr_id2name(a.h, i), bcf_hdr_id2name(b.h, i)))
      FAIL("Contig order mismatch");
  Group x = next_group(&a), y = next_group(&b);
  U common = 0, left = 0, right = 0, changed = 0, filled = 0, compared = 0;
  int32_t *ga = NULL, *gb = NULL;
  int na = 0, nb = 0;
  U it = 0;
  while (x.n || y.n) {
    if (x.n)
      duplicates(&x);
    if (y.n)
      duplicates(&y);
    int cmp = !x.n             ? 1
              : !y.n           ? -1
              : x.rid != y.rid ? (x.rid < y.rid ? -1 : 1)
              : x.pos < y.pos  ? -1
              : x.pos > y.pos  ? 1
                               : 0;
    if (cmp < 0) {
      left += x.n;
      free_group(&x);
      x = next_group(&a);
    } else if (cmp > 0) {
      right += y.n;
      free_group(&y);
      y = next_group(&b);
    } else {
      int *used = calloc(y.n, sizeof(int));
      for (int i = 0; i < x.n; i++) {
        int j = 0;
        for (; j < y.n; j++)
          if (alleles_equal(x.r[i], y.r[j]))
            break;
        if (j == y.n) {
          left++;
          continue;
        }
        used[j] = 1;
        common++;
        int n1 = bcf_get_genotypes(a.h, x.r[i], &ga, &na),
            n2 = bcf_get_genotypes(b.h, y.r[j], &gb, &nb);
        if (n1 <= 0 || n2 <= 0 || n1 % a.ns || n2 % b.ns)
          FAIL("Invalid GT");
        for (int k = 0; k < a.ns; k++) {
          int u[2] = {-1, -1}, v[2] = {-1, -1}, m1, m2, h1, h2;
          int p = genotype(ga + k * (n1 / a.ns), n1 / a.ns, u, &m1, &h1),
              q = genotype(gb + k * (n2 / b.ns), n2 / b.ns, v, &m2, &h2);
          compared++;
          if (p != q || u[0] != v[0] || u[1] != v[1]) {
            changed++;
            b.s[k].changed++;
            if (m1 && h2) {
              filled++;
              b.s[k].filled++;
            }
          }
        }
      }
      for (int j = 0; j < y.n; j++)
        if (!used[j])
          right++;
      free(used);
      free_group(&x);
      free_group(&y);
      x = next_group(&a);
      y = next_group(&b);
    }
    if (++it % 1000000 == 0)
      fprintf(stderr,
              "groups=%" PRIu64 " matched=%" PRIu64 " changed=%" PRIu64 "\n",
              it, common, changed);
  }
  save_stats(&a, argv[3], "left");
  save_stats(&b, argv[3], "right");
  printf("{\"complete\":true,\"htslib\":\"%s\",\"left_records\":%" PRIu64
         ",\"right_records\":%" PRIu64 ",\"matched_records\":%" PRIu64
         ",\"left_only\":%" PRIu64 ",\"right_only\":%" PRIu64
         ",\"compared_genotypes\":%" PRIu64 ",\"changed_gt\":%" PRIu64
         ",\"missing_to_homref\":%" PRIu64 "}\n",
         hts_version(), a.n, b.n, common, left, right, compared, changed,
         filled);
  if (hts_close(a.f) || hts_close(b.f))
    FAIL("VCF close failed");
  return 0;
}
