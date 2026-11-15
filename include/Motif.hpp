#pragma once

#include "Classes.hpp"

#include <array>
#include <istream>
#include <string>
#include <vector>

/*
 * A MEME 4 letter-probability matrix. Scores are the sum of log2(P(base |
 * motif) / P(base | background)) over the motif width.
 */
class motif_matrix {
    public:
    static std::vector<motif_matrix> load(std::istream &in);

    const std::string &name() const;
    int width() const;
    bool supports(char st) const;
    double score(const std::string &s) const;
    double position_score(char b, int p) const;
    double position_score(int b, int p) const;

    private:
    motif_matrix(
        std::string n, int w, bool fp, bool fm,
        std::array<double, 4> bg,
        std::vector<std::array<double, 4>> rows
    );

    static int base_index(char b);

    std::string n_;
    int w_;
    bool fp_, fm_;
    std::array<double, 4> bg_;
    std::vector<std::array<double, 4>> rows_;
    std::vector<std::array<double, 4>> lod_;

    friend class motif_loader;
    friend class motif_scanner;
};

class motif_hit {
    public:
    motif_hit(
        size_t idx, std::string h,
        int is, int ie, std::string m, char st,
        int s, int e, double sc, std::string seq
    );

    size_t region_index() const;
    const std::string &host() const;
    int igr_start() const;
    int igr_end() const;
    const std::string &motif() const;
    char strand() const;
    int start() const;
    int end() const;
    double score() const;
    const std::string &sequence() const;

    private:
    size_t idx_;
    std::string h_, m_, seq_;
    int is_, ie_, s_, e_;
    char st_;
    double sc_;
};

class motif_scanner {
    public:
    static std::vector<motif_hit> scan(
        const region &r, size_t idx,
        const std::vector<motif_matrix> &ms, double cut
    );

    private:
    static std::string reverse_complement(const std::string &s);
};
