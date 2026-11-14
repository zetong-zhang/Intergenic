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
    feature(std::string host, std::string type, int start, int end,
            float score, char strand);

    const std::string &host() const;
    const std::string &type() const;
    int start() const;
    int end() const;
    float score() const;
    char strand() const;
    int length() const;
    bool wrapped() const;

    private:
    std::string host_, type_;
    int start_, end_;
    float score_;
    char strand_;
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
        const std::unordered_set<std::string> &types, float threshold
    );

    void load(std::istream &handle, const std::string &format);
    std::vector<feature> resolve(
        const std::string &host, const std::vector<feature> &embedded
    ) const;

    private:
    static feature parse(const std::string &line, const std::string &format);
    bool accepts(const feature &item) const;
    bool provides(const std::string &host) const;
    const std::vector<feature> &features(const std::string &host) const;

    std::unordered_map<std::string, std::vector<feature>> features_;
    std::unordered_set<std::string> provided_;
    std::unordered_set<std::string> types_;
    float threshold_ = 0.0F;
};

/*
 * An extracted intergenic region. Coordinates and adjacency lists follow the
 * same conventions as feature intervals. The sequence is stored in output
 * orientation, including the joined sequence for an origin-spanning region.
 */
class region {
    public:
    region(
        std::string host, int start, int end,
        std::string sequence, std::string adjacent
    );

    const std::string &host() const;
    int start() const;
    int end() const;
    const std::string &sequence() const;
    const std::string &adjacent() const;
    int length() const;
    bool wrapped() const;

    private:
    std::string host_, sequence_, adjacent_;
    int start_, end_;
};

/*
 * Reads one FASTA/FASTQ/GenBank/EMBL record at a time. The next record remains
 * pending in FASTA streams, allowing arbitrarily large multi-record inputs to
 * be handled incrementally.
 */
class scaffold {
    public:
    bool load(std::istream &handle, bool default_circular = false);
    const std::string &name() const;
    const std::vector<feature> &embedded_features() const;
    std::vector<region> extract(
        const std::vector<feature> &masks, int minlen
    ) const;

    private:
    bool next_line(std::istream &handle, std::string &line);
    void read_header(const std::string &header, bool default_circular);
    void parse_feature_line(const std::string &line);
    bool is_adjacent(int start, int end, const feature &item) const;
    std::string adjacent_features(
        int start, int end, const std::vector<feature> &masks
    ) const;
    void mask(std::vector<char> &covered, const feature &item) const;
    std::string subseq(int start, int end) const;

    static void normalize(std::string &sequence);
    static std::string record_name(const std::string &header);

    std::string pending_, name_, sequence_;
    std::vector<feature> embedded_features_;
    bool circular_ = false;
};
