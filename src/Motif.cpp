#include "Motif.hpp"
#include "Utils.hpp"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdint>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <utility>

namespace {

constexpr double MIN_P = 1e-12;

std::string trim(std::string s) {
    const auto b = s.find_first_not_of(" \t\r\n");
    if (b == std::string::npos) return {};
    const auto e = s.find_last_not_of(" \t\r\n");
    return s.substr(b, e-b+1);
}

std::string unspace(std::string s) {
    s.erase(
        std::remove_if(s.begin(), s.end(),
            [](unsigned char c) { return std::isspace(c); }),
        s.end()
    );
    return s;
}

int base_id(char b) {
    switch (std::toupper((unsigned char) b)) {
        case 'A': return 0;
        case 'C': return 1;
        case 'G': return 2;
        case 'T': return 3;
        default: return -1;
    }
}

char comp(char b) {
    switch (std::toupper((unsigned char) b)) {
        case 'A': return 'T';
        case 'C': return 'G';
        case 'G': return 'C';
        case 'T': return 'A';
        default: return 'N';
    }
}

std::vector<double> nums(const std::string &s) {
    std::istringstream in(s);
    std::vector<double> vs;
    double v = 0.0;
    while (in >> v) vs.push_back(v);
    if (in.fail() && !in.eof())
        logger("Error: invalid MEME probability row");
    return vs;
}

int get_width(const std::string &s) {
    const auto pos = s.find("w=");
    if (pos == std::string::npos)
        logger("Error: MEME matrix has no width (w=)");
    std::istringstream in(s.substr(pos+2));
    int w = 0;
    in >> w;
    if (!in || w <= 0)
        logger("Error: invalid MEME matrix width");
    return w;
}

}  // namespace

class motif_loader {
    public:
    static motif_matrix from_values(
        std::string n, int w, bool fp, bool fm,
        std::array<double, 4> bg,
        std::vector<std::array<double, 4>> rows
    ) {
        return motif_matrix(
            std::move(n), w, fp, fm, bg, std::move(rows)
        );
    }
};

motif_matrix::motif_matrix(
    std::string n, int w, bool fp, bool fm,
    std::array<double, 4> bg,
    std::vector<std::array<double, 4>> rows
): n_(std::move(n)), w_(w), fp_(fp), fm_(fm),
   bg_(bg), rows_(std::move(rows)) {
    lod_.resize(rows_.size());
    for (size_t i = 0; i < rows_.size(); i ++) {
        for (size_t j = 0; j < 4; j ++) {
            const auto p = std::max(rows_[i][j], MIN_P);
            lod_[i][j] = std::log2(p/bg_[j]);
        }
    }
}

const std::string &motif_matrix::name() const { return n_; }
int motif_matrix::width() const { return w_; }
bool motif_matrix::supports(char st) const {
    return st == '+' ? fp_ : st == '-' ? fm_ : false;
}

int motif_matrix::base_index(char b) {
    return base_id(b);
}

double motif_matrix::score(const std::string &s) const {
    if ((int) s.size() != w_) return -std::numeric_limits<double>::infinity();

    double sc = 0.0;
    for (int i = 0; i < w_; i ++) {
        sc += position_score(s[(size_t) i], i);
    }
    return sc;
}

double motif_matrix::position_score(char b, int p) const {
    const auto id = base_index(b);
    return id < 0 ? 0.0 : lod_[(size_t) p][(size_t) id];
}

double motif_matrix::position_score(int b, int p) const {
    return b < 0 ? 0.0 : lod_[(size_t) p][(size_t) b];
}

std::vector<motif_matrix> motif_matrix::load(std::istream &in) {
    std::vector<motif_matrix> ps;
    std::string s;
    std::string abc = "ACGT";
    std::array<double, 4> bg = {0.25, 0.25, 0.25, 0.25};
    bool fp = true, fm = true;
    bool in_bg = false;
    bool has_ver = false;

    std::string name;
    int w = 0;
    std::vector<std::array<double, 4>> rows;

    auto flush = [&]() {
        if (name.empty()) return;
        if (w <= 0 || (int) rows.size() != w)
            logger("Error: incomplete MEME motif '" + name + "'");

        double sum = 0.0;
        for (double v : bg) {
            if (v <= 0.0)
                logger("Error: MEME background probabilities must be positive");
            sum += v;
        }
        for (auto &v : bg) v /= sum;

        ps.push_back(motif_loader::from_values(
            name, w, fp, fm, bg, rows
        ));
        name.clear();
        w = 0;
        rows.clear();
    };

    while (std::getline(in, s)) {
        s = trim(s);
        if (s.empty() || s[0] == '#') continue;

        if (startswith(lower(s), "meme version")) {
            std::istringstream ss(s);
            std::string key, ver;
            ss >> key >> ver >> ver;
            if (ver != "4")
                logger("Error: only MEME version 4 motif files are supported");
            has_ver = true;
            continue;
        }

        if (w > 0 && !name.empty()) {
            const auto vs = nums(s);
            if ((int) vs.size() != (int) abc.size())
                logger("Error: MEME matrix row length does not match ALPHABET");

            std::array<double, 4> r = {0.0, 0.0, 0.0, 0.0};
            for (size_t i = 0; i < vs.size(); i ++) {
                const auto id = base_id(abc[i]);
                if (vs[i] < 0.0)
                    logger("Error: negative probability in MEME matrix");
                r[(size_t) id] = vs[i];
            }
            rows.push_back(r);
            if ((int) rows.size() == w) flush();
            continue;
        }

        if (startswith(s, "ALPHABET=")) {
            abc = unspace(s.substr(9));
            if (abc.size() != 4 ||
                    base_id(abc[0]) < 0 || base_id(abc[1]) < 0 ||
                    base_id(abc[2]) < 0 || base_id(abc[3]) < 0)
                logger("Error: only A/C/G/T MEME alphabets are supported");
            continue;
        }

        if (startswith(lower(s), "strands:")) {
            fp = fm = false;
            std::istringstream ss(s.substr(8));
            std::string st;
            while (ss >> st) {
                if (st == "+") fp = true;
                else if (st == "-") fm = true;
            }
            if (!fp && !fm)
                logger("Error: MEME motif specifies no valid strand");
            continue;
        }

        if (startswith(lower(s), "background letter frequencies")) {
            in_bg = true;
            continue;
        }

        if (in_bg) {
            std::istringstream ss(s);
            std::string b;
            double v = 0.0;
            bool ok = false;
            while (ss >> b >> v) {
                const auto id = base_id(b[0]);
                if (id < 0)
                    logger("Error: invalid base in MEME background");
                bg[(size_t) id] = v;
                ok = true;
            }
            in_bg = !ok;
            if (ok) continue;
        }

        if (startswith(s, "MOTIF ")) {
            flush();
            std::istringstream ss(s.substr(6));
            ss >> name;
            if (name.empty())
                logger("Error: MEME motif has an empty name");
            continue;
        }

        if (startswith(lower(s), "letter-probability matrix:")) {
            if (name.empty())
                logger("Error: MEME matrix appears before a MOTIF declaration");
            w = get_width(s);
            rows.clear();
            rows.reserve((size_t) w);
            continue;
        }
    }

    flush();
    if (!has_ver) logger("Error: motif file is not in MEME 4 format");
    if (ps.empty()) logger("Error: no motifs found in MEME file");
    return ps;
}

motif_hit::motif_hit(
    size_t idx, std::string h,
    int is, int ie, std::string m, char st,
    int s, int e, double sc, std::string seq
): idx_(idx), h_(std::move(h)), m_(std::move(m)), seq_(std::move(seq)),
   is_(is), ie_(ie), s_(s), e_(e), st_(st), sc_(sc) {}

size_t motif_hit::region_index() const { return idx_; }
const std::string &motif_hit::host() const { return h_; }
int motif_hit::igr_start() const { return is_; }
int motif_hit::igr_end() const { return ie_; }
const std::string &motif_hit::motif() const { return m_; }
char motif_hit::strand() const { return st_; }
int motif_hit::start() const { return s_; }
int motif_hit::end() const { return e_; }
double motif_hit::score() const { return sc_; }
const std::string &motif_hit::sequence() const { return seq_; }

std::string motif_scanner::reverse_complement(const std::string &s) {
    std::string r(s.rbegin(), s.rend());
    for (auto &b : r) b = comp(b);
    return r;
}

std::vector<motif_hit> motif_scanner::scan(
    const region &r, size_t idx,
    const std::vector<motif_matrix> &ms, double cut
) {
    std::vector<motif_hit> hits;
    const auto &seq = r.sequence();
    const auto n = r.source_len();
    if (n <= 0) return hits;
    const auto loop = r.circ() && (int) seq.size() == n;

    std::vector<std::int8_t> ids(seq.size());
    for (size_t i = 0; i < seq.size(); i ++)
        ids[i] = (std::int8_t) motif_matrix::base_index(seq[i]);

    for (const auto &m : ms) {
        const auto w = m.width();
        if (w <= 0 || w > (int) seq.size()) continue;
        const auto fp = m.supports('+');
        const auto fm = m.supports('-');
        const auto limit = loop ? n : (int) seq.size()-w+1;

        for (int o = 0; o < limit; o ++) {
            double ps = 0.0;
            double ns = 0.0;
            if (fp && fm) {
                for (int p = 0; p < w; p ++) {
                    auto at = o+p;
                    if (loop && at >= n) at -= n;
                    const auto id = ids[(size_t) at];
                    if (id >= 0) {
                        ps += m.position_score(id, p);
                        ns += m.position_score(3-id, w-p-1);
                    }
                }
            } else if (fp) {
                for (int p = 0; p < w; p ++) {
                    auto at = o+p;
                    if (loop && at >= n) at -= n;
                    const auto id = ids[(size_t) at];
                    if (id >= 0) ps += m.position_score(id, p);
                }
            } else {
                for (int p = 0; p < w; p ++) {
                    auto at = o+p;
                    if (loop && at >= n) at -= n;
                    const auto id = ids[(size_t) at];
                    if (id >= 0)
                        ns += m.position_score(3-id, w-p-1);
                }
            }

            const auto hs = (r.start()+o) % n + 1;
            const auto he = (r.start()+o+w-1) % n + 1;
            auto window = [&]() {
                if (!loop || o+w <= n)
                    return seq.substr((size_t) o, (size_t) w);
                std::string s = seq.substr((size_t) o);
                s.append(seq, 0, (size_t) (o+w-n));
                return s;
            };

            if (fp && ps > cut) {
                    hits.emplace_back(
                        idx, r.host(), r.start()+1, r.end(), m.name(), '+',
                        hs, he, ps, window()
                    );
            }
            if (fm && ns > cut) {
                    hits.emplace_back(
                        idx, r.host(), r.start()+1, r.end(), m.name(), '-',
                        hs, he, ns,
                        reverse_complement(window())
                    );
            }
        }
    }
    return hits;
}
