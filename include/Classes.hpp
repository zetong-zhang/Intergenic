#pragma once

#include <istream>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

/*
 * Annotation intervals are stored with 0-based, end-exclusive coordinates.
 * A wrapped interval has start > end and represents bases from start to the
 * end of its sequence followed by bases from the beginning of the sequence.
 */
class feature {
    public:
    feature(std::string h, std::string t, int s, int e, float sc, char st);

    const std::string &host() const;
    const std::string &type() const;
    int start() const;
    int end() const;
    float score() const;
    char strand() const;
    int length() const;
    bool wrapped() const;

    private:
    std::string h_, t_;
    int s_, e_;
    float sc_;
    char st_;
};

/*
 * Loads target features from GFF/BED and records which sequence names appear
 * in the external file. A recorded name takes precedence over embedded
 * GenBank/EMBL features for that sequence, even if all external features are
 * filtered out by type or score.
 */
class feature_loader {
    public:
    feature_loader() = default;
    feature_loader(
        const std::unordered_set<std::string> &ts, float cut
    );

    void load(std::istream &in, const std::string &fmt);
    std::vector<feature> resolve(
        const std::string &h, const std::vector<feature> &emb
    ) const;

    private:
    static feature parse(const std::string &s, const std::string &fmt);
    bool accepts(const feature &f) const;
    bool provides(const std::string &h) const;
    const std::vector<feature> &features(const std::string &h) const;

    std::unordered_map<std::string, std::vector<feature>> lst_;
    std::unordered_set<std::string> seen_;
    std::unordered_set<std::string> ts_;
    float cut_ = 0.0F;
};

/*
 * An extracted intergenic region. Coordinates and adjacency lists follow the
 * same conventions as feature intervals. The sequence is stored in output
 * orientation, including the joined sequence for an origin-spanning region.
 */
class region {
    public:
    region(
        std::string h, int s, int e,
        std::string seq, std::string near, int len, bool circ
    );

    const std::string &host() const;
    int start() const;
    int end() const;
    const std::string &sequence() const;
    const std::string &adjacent() const;
    int source_len() const;
    bool circ() const;
    int length() const;
    bool wrapped() const;

    private:
    std::string h_, seq_, near_;
    int s_, e_, n_;
    bool circ_;
};

/*
 * Reads one FASTA/FASTQ/GenBank/EMBL record at a time. The next record remains
 * pending in FASTA streams, allowing arbitrarily large multi-record inputs to
 * be handled incrementally.
 */
class scaffold {
    public:
    bool load(std::istream &in, bool circ = false);
    const std::string &name() const;
    const std::vector<feature> &embedded_features() const;
    std::vector<region> extract(
        const std::vector<feature> &ms, int min
    ) const;

    private:
    bool next_line(std::istream &in, std::string &s);
    void read_header(const std::string &h, bool circ);
    void parse_feature_line(const std::string &s);
    bool is_adjacent(int s, int e, const feature &f) const;
    std::string adjacent_features(
        int s, int e, const std::vector<feature> &ms
    ) const;
    void mask(std::vector<char> &cov, const feature &f) const;
    std::string subseq(int s, int e) const;

    static void normalize(std::string &seq);
    static std::string record_name(const std::string &h);

    std::string pend_, name_, seq_;
    std::vector<feature> emb_;
    bool circ_ = false;
};
