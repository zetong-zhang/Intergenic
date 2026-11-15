#include "Classes.hpp"
#include "Utils.hpp"

#include <algorithm>
#include <cctype>
#include <sstream>
#include <stdexcept>
#include <tuple>
#include <utility>

namespace {

constexpr int ADJACENT_DISTANCE = 30;

std::vector<std::string> cols(const std::string &s) {
    std::vector<std::string> out;
    std::istringstream in(s);
    std::string v;
    while (std::getline(in, v, '\t')) out.push_back(v);
    return out;
}

float num(const std::string &s) {
    return s.empty() || s == "." ? 1.0F : std::stof(s);
}

bool is_circ(const std::string &h) {
    return lower(h).find("circular") != std::string::npos;
}

bool is_lin(const std::string &h) {
    return lower(h).find("linear") != std::string::npos;
}

std::string trim(const std::string &s) {
    const auto b = s.find_first_not_of(" \t");
    if (b == std::string::npos) return {};
    const auto e = s.find_last_not_of(" \t");
    return s.substr(b, e-b+1);
}

bool location_bounds(
    const std::string &loc, int &s, int &e
) {
    // GenBank/EMBL locations may use complement(), join(), or order(). For
    // masking, the first coordinate through the last coordinate represents
    // the complete target span; this also preserves a circular wrapped span.
    std::vector<int> ns;
    for (size_t i = 0; i < loc.size();) {
        if (!std::isdigit((unsigned char) loc[i])) {
            i ++;
            continue;
        }
        size_t j = i;
        while (j < loc.size() &&
                std::isdigit((unsigned char) loc[j])) j ++;
        try {
            ns.push_back(std::stoi(loc.substr(i, j-i)));
        } catch (const std::exception &) {
            return false;
        }
        i = j;
    }
    if (ns.empty()) return false;
    s = ns.front()-1;
    e = ns.back();
    return s >= 0 && e > 0;
}

}  // namespace

feature::feature(
    std::string h, std::string t, int s, int e, float sc, char st
): h_(std::move(h)), t_(std::move(t)), s_(s), e_(e), sc_(sc), st_(st) {}

const std::string &feature::host() const { return h_; }
const std::string &feature::type() const { return t_; }
int feature::start() const { return s_; }
int feature::end() const { return e_; }
float feature::score() const { return sc_; }
char feature::strand() const { return st_; }
int feature::length() const { return e_-s_; }
bool feature::wrapped() const { return s_ > e_; }

feature_loader::feature_loader(
    const std::unordered_set<std::string> &ts, float cut
): ts_(ts), cut_(cut) {}

feature feature_loader::parse(
    const std::string &s, const std::string &fmt
) {
    const auto v = cols(s);
    const auto kind = lower(fmt);
    try {
        if (kind == "gff") {
            if (v.size() != 9)
                logger("Error: incorrect GFF file format (n_cols != 9).");
            if (v[3].empty() || v[4].empty() || v[6].empty())
                logger("Error: missing required GFF field (" + s + ")");
            return feature(v[0], v[2], std::stoi(v[3])-1,
                           std::stoi(v[4]), num(v[5]), v[6][0]);
        }
        if (kind == "bed") {
            if (v.size() < 3)
                logger("Error: incorrect BED file format (n_cols < 3).");
            const auto sc = v.size() >= 5 ? num(v[4]) : 1.0F;
            const auto st = v.size() >= 6 && !v[5].empty() ?
                            v[5][0] : '.';
            return feature(v[0], "CDS", std::stoi(v[1]),
                           std::stoi(v[2]), sc, st);
        }
    } catch (const std::exception &) {
        logger("Error: invalid " + fmt + " line (" + s + ")");
    }
    logger("Error: unsupported annotation file format");
}

void feature_loader::load(
    std::istream &in, const std::string &fmt
) {
    lst_.clear();
    seen_.clear();
    std::string s;
    while (std::getline(in, s)) {
        if (!s.empty() && s.back() == '\r') s.pop_back();
        if (s.empty() || s[0] == '#') continue;
        const auto t = lower(s);
        if (startswith(t, "track") || startswith(t, "browser"))
            continue;

        auto f = parse(s, fmt);
        seen_.insert(f.host());
        if (accepts(f))
            lst_[f.host()].push_back(std::move(f));
    }
}

bool feature_loader::provides(const std::string &h) const {
    return seen_.count(h) != 0;
}

bool feature_loader::accepts(const feature &f) const {
    return ts_.count(f.type()) && f.score() >= cut_;
}

const std::vector<feature> &feature_loader::features(
    const std::string &h
) const {
    static const std::vector<feature> none;
    const auto it = lst_.find(h);
    return it == lst_.end() ? none : it->second;
}

std::vector<feature> feature_loader::resolve(
    const std::string &h, const std::vector<feature> &emb
) const {
    if (provides(h)) return features(h);

    std::vector<feature> ms;
    for (const auto &f : emb) {
        if (accepts(f)) ms.push_back(f);
    }
    return ms;
}

region::region(
    std::string h, int s, int e,
    std::string seq, std::string near, int len, bool circ
): h_(std::move(h)), seq_(std::move(seq)),
   near_(std::move(near)), s_(s), e_(e), n_(len), circ_(circ) {}

const std::string &region::host() const { return h_; }
int region::start() const { return s_; }
int region::end() const { return e_; }
const std::string &region::sequence() const { return seq_; }
const std::string &region::adjacent() const { return near_; }
int region::source_len() const { return n_; }
bool region::circ() const { return circ_; }
int region::length() const { return (int) seq_.size(); }
bool region::wrapped() const { return s_ > e_; }

bool scaffold::next_line(std::istream &in, std::string &s) {
    if (!pend_.empty()) {
        s.swap(pend_);
        return true;
    }
    if (!std::getline(in, s)) return false;
    if (!s.empty() && s.back() == '\r') s.pop_back();
    return true;
}

void scaffold::normalize(std::string &seq) {
    size_t n = 0;
    for (char raw : seq) {
        const auto c = (unsigned char) raw;
        if (std::isalpha(c) || c == '-' || c == '.')
            seq[n++] = (char) std::toupper(c);
    }
    seq.resize(n);
}

std::string scaffold::record_name(const std::string &h) {
    size_t b = startswith(h, "LOCUS") ? 5u : 2u;
    while (b < h.size() && std::isspace((unsigned char) h[b])) b ++;
    const auto e = h.find_first_of(" \t;", b);
    return h.substr(b, e == std::string::npos ?
                    std::string::npos : e-b);
}

void scaffold::read_header(
    const std::string &h, bool circ
) {
    size_t b = h[0] == '>' || h[0] == '@' ? 1u : 0u;
    while (b < h.size() && std::isspace((unsigned char) h[b])) b ++;
    const auto e = h.find_first_of(" \t", b);
    name_ = h.substr(b, e == std::string::npos ?
                     std::string::npos : e-b);
    if (name_.empty()) logger("Error: invalid sequence header line");

    if (is_circ(h)) circ_ = true;
    else if (is_lin(h)) circ_ = false;
    else circ_ = circ;
}

void scaffold::parse_feature_line(const std::string &s) {
    if (s.size() <= 21 || s.compare(0, 5, "     ") != 0 ||
            std::isspace((unsigned char) s[5])) return;

    const auto t = trim(s.substr(5, 15));
    const auto loc = trim(s.substr(21));
    int b = 0, e = 0;
    if (!t.empty() && !loc.empty() && location_bounds(loc, b, e)) {
        emb_.emplace_back(name_, t, b, e, 1.0F, '.');
    }
}

bool scaffold::is_adjacent(
    int s, int e, const feature &f
) const {
    // "Adjacent" means an endpoint is at most 30 bp from either IGR boundary.
    // The circular comparison is performed with wrapped coordinate copies so
    // both IGRs and features may cross the origin.
    const auto n = (int) seq_.size();
    auto fs = f.start();
    auto fe = f.end();
    if (fe < fs) fe += n;

    if (!circ_) {
        const auto lo = std::max(0, s-ADJACENT_DISTANCE);
        const auto hi = std::min(n, e+ADJACENT_DISTANCE);
        return fs <= hi && fe >= lo;
    }

    auto re = e;
    if (re < s) re += n;
    const auto lo = s-ADJACENT_DISTANCE;
    const auto hi = re+ADJACENT_DISTANCE;
    for (auto d = -n; d <= n; d += n) {
        if (fs+d <= hi && fe+d >= lo) return true;
    }
    return false;
}

std::string scaffold::adjacent_features(
    int s, int e, const std::vector<feature> &ms
) const {
    // Emit complete, user-readable target coordinates in stable genomic order.
    // Output coordinates are 1-based and inclusive.
    using coordinate = std::tuple<std::string, int, int, char>;
    std::vector<coordinate> cs;
    for (const auto &f : ms) {
        if (is_adjacent(s, e, f))
            cs.emplace_back(
                f.host(), f.start(), f.end(), f.strand()
            );
    }
    std::sort(cs.begin(), cs.end());
    cs.erase(
        std::unique(cs.begin(), cs.end()),
        cs.end()
    );
    if (cs.empty()) return ".";

    std::ostringstream out;
    for (size_t i = 0; i < cs.size(); i ++) {
        if (i) out << ',';
        const auto &[h, b, f, st] = cs[i];
        out << h << ':' << b+1 << ".." << f << '(' << st << ')';
    }
    return out.str();
}

void scaffold::mask(
    std::vector<char> &cov, const feature &f
) const {
    // Clamp annotations to the sequence, and mask both pieces of a wrapped
    // feature on a circular scaffold.
    if (f.length() == 0) return;
    const auto n = (int) seq_.size();
    if (circ_ && f.wrapped()) {
        const auto b = std::clamp(f.start(), 0, n);
        const auto e = std::clamp(f.end(), 0, n);
        std::fill(cov.begin()+b, cov.end(), true);
        std::fill(cov.begin(), cov.begin()+e, true);
        return;
    }

    const auto b = std::clamp(f.start(), 0, n);
    const auto e = std::clamp(f.end(), 0, n);
    if (b < e) std::fill(cov.begin()+b, cov.begin()+e, true);
}

std::string scaffold::subseq(int s, int e) const {
    if (s <= e) return seq_.substr((size_t) s, (size_t) (e-s));
    return seq_.substr((size_t) s) + seq_.substr(0, (size_t) e);
}

std::vector<region> scaffold::extract(
    const std::vector<feature> &ms, int min
) const {
    // Build the unmasked runs first. On a circular sequence, runs touching the
    // two sequence ends are one region crossing the origin.
    const auto n = (int) seq_.size();
    if (n <= 0) return {};

    std::vector<char> cov((size_t) n, false);
    for (const auto &f : ms) mask(cov, f);

    std::vector<std::pair<int, int>> runs;
    for (int i = 0; i < n;) {
        while (i < n && cov[(size_t) i]) i ++;
        if (i == n) break;
        const auto s = i;
        while (i < n && !cov[(size_t) i]) i ++;
        runs.emplace_back(s, i);
    }

    if (circ_ && runs.size() >= 2 && runs.front().first == 0 &&
            runs.back().second == n) {
        runs.front().first = runs.back().first;
        runs.pop_back();
    }

    std::vector<region> results;
    for (const auto &[s, e] : runs) {
        // Apply the minimum length after merging circular end runs so the
        // reported length matches the complete biological IGR.
        const auto len = s <= e ? e-s : e+n-s;
        if (len >= min)
            results.emplace_back(
                name_, s, e, subseq(s, e),
                adjacent_features(s, e, ms), n, circ_
            );
    }
    return results;
}

bool scaffold::load(std::istream &in, bool circ) {
    seq_.clear();
    name_.clear();
    emb_.clear();
    circ_ = circ;

    std::string s;
    if (!next_line(in, s)) return false;
    while (s.empty()) {
        if (!next_line(in, s)) return false;
    }

    if (s[0] == '>') {
        // FASTA: sequence continues until the next header.
        read_header(s, circ);
        while (next_line(in, s)) {
            if (!s.empty() && s[0] == '>') {
                pend_ = s;
                break;
            }
            seq_ += s;
        }
        normalize(seq_);
    } else if (s[0] == '@') {
        // FASTQ: ignore quality values after verifying their length.
        read_header(s, circ);
        bool sep = false;
        while (next_line(in, s)) {
            if (!s.empty() && s[0] == '+') {
                sep = true;
                break;
            }
            seq_ += s;
        }
        if (!sep) logger("Error: incomplete FASTQ record");
        normalize(seq_);

        std::string q;
        while (q.size() < seq_.size() && next_line(in, s)) q += s;
        if (q.size() != seq_.size())
            logger("Error: FASTQ sequence and quality lengths differ");
    } else if (startswith(s, "LOCUS") || startswith(s, "ID")) {
        // GenBank/EMBL: parse embedded feature locations until ORIGIN/SQ,
        // then collect sequence lines until the record terminator.
        name_ = record_name(s);
        if (name_.empty()) logger("Error: invalid LOCUS/ID line");
        if (is_circ(s)) circ_ = true;
        else if (is_lin(s)) circ_ = false;
        else circ_ = circ;

        bool seq_on = false;
        bool feat_on = false;
        while (next_line(in, s)) {
            if (!seq_on && startswith(s, "FT")) {
                parse_feature_line("     " + s.substr(5));
                continue;
            }
            if (!seq_on && startswith(s, "FEATURES")) {
                feat_on = true;
                continue;
            }
            if (!seq_on && (startswith(s, "ORIGIN") || startswith(s, "SQ"))) {
                seq_on = true;
                feat_on = false;
            } else if (s == "//") {
                break;
            } else if (!seq_on && feat_on) {
                parse_feature_line(s);
            } else if (seq_on) {
                seq_ += s;
            }
        }
        if (!seq_on)
            logger("Error: no sequence section in GenBank/EMBL record");
        normalize(seq_);
    } else {
        logger("Error: unsupported input sequence format.");
    }
    return true;
}

const std::string &scaffold::name() const { return name_; }

const std::vector<feature> &scaffold::embedded_features() const {
    return emb_;
}
