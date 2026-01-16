#include <iostream>
#include <vector>
#include <set>
#include <vector>
#include <algorithm>
#include <cmath>
#include <unordered_map>
#include <unordered_set>
#include <map>
#include <queue>
#include <stack>
#include <deque>
#include <random>
#include <bitset>
#include <cassert>
#include <fstream>
#include <chrono>
#include <sstream>
#include <omp.h>
#include <iomanip>

using namespace std;
using ll = long long;
using pii = pair<int, int>;
using vi = vector<int>;
using vs = vector<string>;

string g;
vi sa;
string bwt;
int C[5];
vector<vi> occ;

int ctoi(char c) {
    if (c == '$') return 0;
    if (c == 'A') return 1;
    if (c == 'C') return 2;
    if (c == 'G') return 3;
    if (c == 'T') return 4;
    assert(0);
}

vi build_sa(string s) {
    int n = (int)s.size();
    vi c1(n, -1), c2(n, -1), srt1(n, -1), srt2(n, -1), tmp(n, -1), inv(n, -1), cnt(n + 1, 0);

    vector<pair<char, int>> pre;
    pre.reserve(n);
    for (int i = 0; i < n; i++) pre.emplace_back(s[i], i);
    sort(pre.begin(), pre.end());

    int cur = 0, mx = -1;
    for (auto [ch, ind] : pre) {
        srt1[ind] = cur++;
        c1[ind] = (cur == 1) ? 0 : c1[pre[cur - 2].second] + (ch != pre[cur - 2].first);
        mx = max(mx, c1[ind]);
    }

    for (int l = 1; l <= n; l <<= 1) {
        fill(cnt.begin(), cnt.begin() + mx + 1, 0);
        for (int i = 0; i < n; i++) if (c1[i] < mx) cnt[c1[i] + 1]++;
        for (int j = 1; j <= mx; j++) cnt[j] += cnt[j - 1];

        for (int i = 0; i < n; i++) {
            int j = (i + l >= n) ? i + l - n : i + l;
            tmp[i] = cnt[c1[j]]++;
        }
        for (int i = 0; i < n; i++) inv[tmp[i]] = i;

        for (int j = mx; j > 0; j--) cnt[j] = cnt[j - 1];
        cnt[0] = 0;

        for (int i = 0; i < n; i++) {
            int ind = inv[i];
            srt2[ind] = cnt[c1[ind]]++;
        }
        mx = 0;
        for (int i = 0; i < n; i++) inv[srt2[i]] = i;
        for (int i = 0; i < n; i++) {
            int ind = inv[i];
            int j = (ind + l >= n) ? ind + l - n : ind + l;
            if (i == 0) c2[ind] = 0;
            else {
                int j2 = (inv[i - 1] + l >= n) ? inv[i - 1] + l - n : inv[i - 1] + l;
                c2[ind] = c2[inv[i - 1]] + (c1[ind] != c1[inv[i - 1]] || c1[j] != c1[j2]);
            }
            mx = max(mx, c2[ind]);
        }
        swap(c1, c2);
        swap(srt1, srt2);
        if (mx == n - 1) break;
    }

    vi res(n);
    for (int i = 0; i < n; i++) res[srt1[i]] = i;
    return res;
}

string build_bwt(const string& s) {
    int n = (int)s.size();
    string res(n, ' ');
    for (int i = 0; i < n; i++) res[i] = (sa[i] == 0) ? '$' : s[sa[i] - 1];
    return res;
}

pii bsearch(const string& p) {
    int n = (int)g.size(), lo = 0, hi = n;
    for (int i = (int)p.size() - 1; i >= 0; i--) {
        int idx = ctoi(p[i]);
        lo = C[idx] + occ[idx][lo];
        hi = C[idx] + occ[idx][hi];
        if (lo >= hi) return {-1, -1};
    }
    return {lo, hi};
}

string rc(string s) {
    reverse(s.begin(), s.end());
    for (char& c : s) {
        if (c == 'A') c = 'T';
        else if (c == 'T') c = 'A';
        else if (c == 'C') c = 'G';
        else if (c == 'G') c = 'C';
    }
    return s;
}

pii longest_clean(const string& s) {
    int best_st = 0, best_len = 0, cur_st = 0;
    for (int i = 0; i <= (int)s.size(); i++) {
        if (i == (int)s.size() || s[i] == 'N') {
            int len = i - cur_st;
            if (len > best_len) { best_len = len; best_st = cur_st; }
            cur_st = i + 1;
        }
    }
    return {best_st, best_len};
}

bool verify(const string& rd, int st, int len) {
    if (st < 0 || st + len > (int)g.size() - 1) return false;
    for (int i = 0; i < len; i++) if (rd[i] != 'N' && rd[i] != g[st + i]) return false;
    return true;
}

const int MAX_RDS = 23000000, MAX_OCC = 100000, RLEN = 100;
static char io_buf[1 << 20];

int main() {
    ios_base::sync_with_stdio(false);
    cin.tie(nullptr);

    auto t0 = chrono::high_resolution_clock::now();

    ifstream fa("../GCF_000005845.2_ASM584v2_genomic.fna");
    string line;
    while (getline(fa, line)) if (!line.empty() && line[0] != '>') g += line;
    g += '$';
    fa.close();
    int glen = (int)g.size() - 1;

    sa = build_sa(g);
    bwt = build_bwt(g);
    occ.assign(5, vi(glen + 2, 0));
    for (int i = 0; i <= glen; i++) {
        for (int j = 0; j < 5; j++) occ[j][i + 1] = occ[j][i];
        int idx = ctoi(bwt[i]);
        if (idx >= 0) occ[idx][i + 1]++;
    }
    C[0] = 0;
    for (int i = 1; i < 5; i++) C[i] = C[i - 1] + occ[i - 1][glen + 1];

    vs rds;
    rds.reserve(MAX_RDS);
    ifstream fq("../ERR022075_1.fastq");
    fq.rdbuf()->pubsetbuf(io_buf, sizeof(io_buf));
    string cur;
    while (getline(fq, line)) {
        getline(fq, cur);
        getline(fq, line);
        getline(fq, line);
        rds.push_back(std::move(cur));
    }
    fq.close();

    int n = (int)rds.size();
    vi cnt1(n, -1), cnt2(n, -1), diff(glen + 1, 0);
    vs buf(n);

    #pragma omp parallel for schedule(dynamic, 10000)
    for (int i = 0; i < n; i++) {
        const string& seq = rds[i];

        auto [sub_st, sub_len] = longest_clean(seq);
        string sub = seq.substr(sub_st, sub_len);
        auto [lo1, hi1] = bsearch(sub);
        int occ1 = (lo1 == -1) ? 0 : hi1 - lo1;

        string rseq = rc(seq);
        auto [rc_st, rc_len] = longest_clean(rseq);
        string rc_sub = rseq.substr(rc_st, rc_len);
        auto [lo2, hi2] = bsearch(rc_sub);
        int occ2 = (lo2 == -1) ? 0 : hi2 - lo2;
        if (occ1 + occ2 > MAX_OCC) continue;

        vi pos1, pos2;

        if (lo1 != -1) {
            for (int j = lo1; j < hi1; j++) {
                int p = sa[j] - sub_st;
                if (verify(seq, p, RLEN)) pos1.push_back(p);
            }
        }

        if (lo2 != -1) {
            for (int j = lo2; j < hi2; j++) {
                int p = sa[j] - rc_st;
                if (verify(rseq, p, RLEN)) pos2.push_back(p);
            }
        }

        int c1 = (int)pos1.size(), c2 = (int)pos2.size();
        cnt1[i] = c1;
        cnt2[i] = c2;

        ostringstream oss;
        oss << i << " | " << c1 << " ";
        bool first = true;
        for (int p : pos1) { if (!first) oss << ","; oss << p; first = false; }
        oss << " | " << c2 << " ";
        first = true;
        for (int p : pos2) { if (!first) oss << ","; oss << p; first = false; }
        oss << "\n";
        buf[i] = oss.str();

        if (c1 + c2 == 1) {
            int pos = (c1 == 1) ? pos1[0] : pos2[0];
            #pragma omp atomic
            diff[pos]++;
            #pragma omp atomic
            diff[pos + RLEN]--;
        }
    }

    ofstream res("../results.txt");
    res << "ID | ORIG | COMP\n";
    for (int i = 0; i < n; i++) if (cnt1[i] >= 0) res << buf[i];
    res.close();

    ll skip = 0, unmap = 0, uniq = 0, m2 = 0, m3 = 0;
    for (int i = 0; i < n; i++) {
        int c1 = cnt1[i], c2 = cnt2[i];
        int tot = c1 + c2;
        if (tot < 0) { skip++; continue; }
        if (tot == 0) unmap++;
        else if (tot == 1) uniq++;
        else if (tot == 2) m2++;
        else m3++;
    }

    int covered = 0;
    int pref_sum = 0;
    for (int i = 0; i < glen; i++) { pref_sum += diff[i]; if (pref_sum > 0) covered++; }
    double cov = 100.0 * covered / glen;

    auto t1 = chrono::high_resolution_clock::now();
    double sec = chrono::duration_cast<chrono::milliseconds>(t1 - t0).count() / 1000.0;

    ll mapped = uniq + m2 + m3;
    ll proc = n - skip;

    cout << "Genome size: " << glen << " chars" << endl;
    cout << "Total reads: " << n << endl;
    cout << "Unmapped: " << unmap << " (" << 100.0 * unmap / proc << "%)" << endl;
    cout << "Mapped: " << mapped << " (" << 100.0 * mapped / proc << "%)" << endl;
    cout << "  Unique (1 occurrence): " << uniq << " (" << 100.0 * uniq / proc << "%)" << endl;
    cout << "  Double (2 occurrences): " << m2 << " (" << 100.0 * m2 / proc << "%)" << endl;
    cout << "  Multi (3+ occurrences): " << m3 << " (" << 100.0 * m3 / proc << "%)" << endl;
    cout << "Alignment quality: 100%, exact matching only, \'N\' is a wildcard" << endl;
    cout << "Genome coverage: " << cov << "%" << endl;
    cout << "Total time: " << fixed << setprecision(1) << sec << " s" << endl;
}