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

std::vector<std::string> fields(const std::string &line) {
    std::vector<std::string> items;
    std::istringstream stream(line);
    std::string item;
    while (std::getline(stream, item, '\t')) items.push_back(item);
    return items;
}

float score_of(const std::string &value) {
    return value.empty() || value == "." ? 1.0F : std::stof(value);
}

bool covered_by_header(const std::string &header) {
    return lower(header).find("circular") != std::string::npos;
}

bool linear_by_header(const std::string &header) {
    return lower(header).find("linear") != std::string::npos;
}

std::string trim(const std::string &text) {
    const auto begin = text.find_first_not_of(" \t");
    if (begin == std::string::npos) return {};
    const auto end = text.find_last_not_of(" \t");
    return text.substr(begin, end-begin+1);
}

bool location_bounds(
    const std::string &location, int &start, int &end
) {
    // GenBank/EMBL locations may use complement(), join(), or order(). For
    // masking, the first coordinate through the last coordinate represents
    // the complete target span; this also preserves a circular wrapped span.
    std::vector<int> coordinates;
    for (size_t i = 0; i < location.size();) {
        if (!std::isdigit((unsigned char) location[i])) {
            i ++;
            continue;
        }
        size_t finish = i;
        while (finish < location.size() &&
                std::isdigit((unsigned char) location[finish])) finish ++;
        try {
            coordinates.push_back(std::stoi(location.substr(i, finish-i)));
        } catch (const std::exception &) {
            return false;
        }
        i = finish;
    }
    if (coordinates.empty()) return false;
    start = coordinates.front()-1;
    end = coordinates.back();
    return start >= 0 && end > 0;
}

}  // namespace

feature::feature(
    std::string host, std::string type, int start, int end,
    float score, char strand
): host_(std::move(host)), type_(std::move(type)), start_(start),
   end_(end), score_(score), strand_(strand) {}

const std::string &feature::host() const { return host_; }
const std::string &feature::type() const { return type_; }
int feature::start() const { return start_; }
int feature::end() const { return end_; }
float feature::score() const { return score_; }
char feature::strand() const { return strand_; }
int feature::length() const { return end_-start_; }
bool feature::wrapped() const { return start_ > end_; }

feature_loader::feature_loader(
    const std::unordered_set<std::string> &types, float threshold
): types_(types), threshold_(threshold) {}

feature feature_loader::parse(
    const std::string &line, const std::string &format
) {
    const auto items = fields(line);
    const auto kind = lower(format);
    try {
        if (kind == "gff") {
            if (items.size() != 9)
                logger("Error: incorrect GFF file format (n_cols != 9).");
            if (items[3].empty() || items[4].empty() || items[6].empty())
                logger("Error: missing required GFF field (" + line + ")");
            return feature(items[0], items[2], std::stoi(items[3])-1,
                           std::stoi(items[4]), score_of(items[5]), items[6][0]);
        }
        if (kind == "bed") {
            if (items.size() < 3)
                logger("Error: incorrect BED file format (n_cols < 3).");
            const auto score = items.size() >= 5 ? score_of(items[4]) : 1.0F;
            const auto strand = items.size() >= 6 && !items[5].empty() ?
                                items[5][0] : '.';
            return feature(items[0], "CDS", std::stoi(items[1]),
                           std::stoi(items[2]), score, strand);
        }
    } catch (const std::exception &) {
        logger("Error: invalid " + format + " line (" + line + ")");
    }
    logger("Error: unsupported annotation file format");
}

void feature_loader::load(
    std::istream &handle, const std::string &format
) {
    features_.clear();
    provided_.clear();
    std::string line;
    while (std::getline(handle, line)) {
        if (!line.empty() && line.back() == '\r') line.pop_back();
        if (line.empty() || line[0] == '#') continue;
        const auto check = lower(line);
        if (startswith(check, "track") || startswith(check, "browser"))
            continue;

        auto item = parse(line, format);
        provided_.insert(item.host());
        if (accepts(item))
            features_[item.host()].push_back(std::move(item));
    }
}

bool feature_loader::provides(const std::string &host) const {
    return provided_.count(host) != 0;
}

bool feature_loader::accepts(const feature &item) const {
    return types_.count(item.type()) && item.score() >= threshold_;
}

const std::vector<feature> &feature_loader::features(
    const std::string &host
) const {
    static const std::vector<feature> empty;
    const auto found = features_.find(host);
    return found == features_.end() ? empty : found->second;
}

std::vector<feature> feature_loader::resolve(
    const std::string &host, const std::vector<feature> &embedded
) const {
    if (provides(host)) return features(host);

    std::vector<feature> masks;
    for (const auto &item : embedded) {
        if (accepts(item)) masks.push_back(item);
    }
    return masks;
}

region::region(
    std::string host, int start, int end,
    std::string sequence, std::string adjacent
): host_(std::move(host)), sequence_(std::move(sequence)),
   adjacent_(std::move(adjacent)), start_(start), end_(end) {}

const std::string &region::host() const { return host_; }
int region::start() const { return start_; }
int region::end() const { return end_; }
const std::string &region::sequence() const { return sequence_; }
const std::string &region::adjacent() const { return adjacent_; }
int region::length() const { return (int) sequence_.size(); }
bool region::wrapped() const { return start_ > end_; }

bool scaffold::next_line(std::istream &handle, std::string &line) {
    if (!pending_.empty()) {
        line.swap(pending_);
        return true;
    }
    if (!std::getline(handle, line)) return false;
    if (!line.empty() && line.back() == '\r') line.pop_back();
    return true;
}

void scaffold::normalize(std::string &sequence) {
    std::string result;
    result.reserve(sequence.size());
    for (unsigned char c : sequence) {
        if (std::isalpha(c) || c == '-' || c == '.')
            result.push_back((char) std::toupper(c));
    }
    sequence.swap(result);
}

std::string scaffold::record_name(const std::string &header) {
    auto begin = startswith(header, "LOCUS") ? 5 : 2;
    while (begin < (int) header.size() && std::isspace(header[begin])) begin ++;
    const auto finish = header.find_first_of(" \t;", begin);
    return header.substr(begin, finish == std::string::npos ?
                         std::string::npos : finish-begin);
}

void scaffold::read_header(
    const std::string &header, bool default_circular
) {
    auto begin = header[0] == '>' || header[0] == '@' ? 1 : 0;
    while (begin < (int) header.size() && std::isspace(header[begin])) begin ++;
    const auto finish = header.find_first_of(" \t", begin);
    name_ = header.substr(begin, finish == std::string::npos ?
                          std::string::npos : finish-begin);
    if (name_.empty()) logger("Error: invalid sequence header line");

    if (covered_by_header(header)) circular_ = true;
    else if (linear_by_header(header)) circular_ = false;
    else circular_ = default_circular;
}

void scaffold::parse_feature_line(const std::string &line) {
    if (line.size() <= 21 || line.compare(0, 5, "     ") != 0 ||
            std::isspace((unsigned char) line[5])) return;

    const auto type = trim(line.substr(5, 15));
    const auto location = trim(line.substr(21));
    int start = 0, end = 0;
    if (!type.empty() && !location.empty() &&
            location_bounds(location, start, end)) {
        embedded_features_.emplace_back(name_, type, start, end, 1.0F, '.');
    }
}

bool scaffold::is_adjacent(
    int start, int end, const feature &item
) const {
    // "Adjacent" means an endpoint is at most 30 bp from either IGR boundary.
    // The circular comparison is performed with wrapped coordinate copies so
    // both IGRs and features may cross the origin.
    const auto length = (int) sequence_.size();
    auto feature_start = item.start();
    auto feature_end = item.end();
    if (feature_end < feature_start) feature_end += length;

    if (!circular_) {
        const auto zone_start = std::max(0, start-ADJACENT_DISTANCE);
        const auto zone_end = std::min(length, end+ADJACENT_DISTANCE);
        return feature_start <= zone_end && feature_end >= zone_start;
    }

    auto region_end = end;
    if (region_end < start) region_end += length;
    const auto zone_start = start-ADJACENT_DISTANCE;
    const auto zone_end = region_end+ADJACENT_DISTANCE;
    for (auto shift = -length; shift <= length; shift += length) {
        if (feature_start+shift <= zone_end &&
                feature_end+shift >= zone_start) return true;
    }
    return false;
}

std::string scaffold::adjacent_features(
    int start, int end, const std::vector<feature> &masks
) const {
    // Emit complete, user-readable target coordinates in stable genomic order.
    // Output coordinates are 1-based and inclusive.
    using coordinate = std::tuple<std::string, int, int, char>;
    std::vector<coordinate> coordinates;
    for (const auto &item : masks) {
        if (is_adjacent(start, end, item))
            coordinates.emplace_back(
                item.host(), item.start(), item.end(), item.strand()
            );
    }
    std::sort(coordinates.begin(), coordinates.end());
    coordinates.erase(
        std::unique(coordinates.begin(), coordinates.end()),
        coordinates.end()
    );
    if (coordinates.empty()) return ".";

    std::ostringstream result;
    for (size_t i = 0; i < coordinates.size(); i ++) {
        if (i) result << ',';
        const auto &[host, feature_start, feature_end, strand] =
            coordinates[i];
        result << host << ':' << feature_start+1 << ".."
               << feature_end << '(' << strand << ')';
    }
    return result.str();
}

void scaffold::mask(
    std::vector<char> &covered, const feature &item
) const {
    // Clamp annotations to the sequence, and mask both pieces of a wrapped
    // feature on a circular scaffold.
    if (item.length() == 0) return;
    const auto length = (int) sequence_.size();
    if (circular_ && item.wrapped()) {
        const auto begin = std::clamp(item.start(), 0, length);
        const auto finish = std::clamp(item.end(), 0, length);
        std::fill(covered.begin()+begin, covered.end(), true);
        std::fill(covered.begin(), covered.begin()+finish, true);
        return;
    }

    const auto begin = std::clamp(item.start(), 0, length);
    const auto finish = std::clamp(item.end(), 0, length);
    if (begin < finish)
        std::fill(covered.begin()+begin, covered.begin()+finish, true);
}

std::string scaffold::subseq(int start, int end) const {
    if (start <= end) return sequence_.substr(start, end-start);
    return sequence_.substr(start) + sequence_.substr(0, end);
}

std::vector<region> scaffold::extract(
    const std::vector<feature> &masks, int minlen
) const {
    // Build the unmasked runs first. On a circular sequence, runs touching the
    // two sequence ends are one region crossing the origin.
    const auto length = (int) sequence_.size();
    if (length <= 0) return {};

    std::vector<char> covered(length, false);
    for (const auto &item : masks) mask(covered, item);

    std::vector<std::pair<int, int>> runs;
    for (int i = 0; i < length;) {
        while (i < length && covered[i]) i ++;
        if (i == length) break;
        const auto start = i;
        while (i < length && !covered[i]) i ++;
        runs.emplace_back(start, i);
    }

    if (circular_ && runs.size() >= 2 && runs.front().first == 0 &&
            runs.back().second == length) {
        runs.front().first = runs.back().first;
        runs.pop_back();
    }

    std::vector<region> results;
    for (const auto &[start, end] : runs) {
        // Apply the minimum length after merging circular end runs so the
        // reported length matches the complete biological IGR.
        const auto size = start <= end ? end-start : end+length-start;
        if (size >= minlen)
            results.emplace_back(
                name_, start, end, subseq(start, end),
                adjacent_features(start, end, masks)
            );
    }
    return results;
}

bool scaffold::load(std::istream &handle, bool default_circular) {
    sequence_.clear();
    name_.clear();
    embedded_features_.clear();
    circular_ = default_circular;

    std::string line;
    if (!next_line(handle, line)) return false;
    while (line.empty()) {
        if (!next_line(handle, line)) return false;
    }

    if (line[0] == '>') {
        // FASTA: sequence continues until the next header.
        read_header(line, default_circular);
        while (next_line(handle, line)) {
            if (!line.empty() && line[0] == '>') {
                pending_ = line;
                break;
            }
            sequence_ += line;
        }
        normalize(sequence_);
    } else if (line[0] == '@') {
        // FASTQ: ignore quality values after verifying their length.
        read_header(line, default_circular);
        bool separator_found = false;
        while (next_line(handle, line)) {
            if (!line.empty() && line[0] == '+') {
                separator_found = true;
                break;
            }
            sequence_ += line;
        }
        if (!separator_found) logger("Error: incomplete FASTQ record");
        normalize(sequence_);

        std::string quality;
        while (quality.size() < sequence_.size() && next_line(handle, line))
            quality += line;
        if (quality.size() != sequence_.size())
            logger("Error: FASTQ sequence and quality lengths differ");
    } else if (startswith(line, "LOCUS") || startswith(line, "ID")) {
        // GenBank/EMBL: parse embedded feature locations until ORIGIN/SQ,
        // then collect sequence lines until the record terminator.
        name_ = record_name(line);
        if (name_.empty()) logger("Error: invalid LOCUS/ID line");
        if (covered_by_header(line)) circular_ = true;
        else if (linear_by_header(line)) circular_ = false;
        else circular_ = default_circular;

        bool sequence_started = false;
        bool in_features = false;
        while (next_line(handle, line)) {
            if (!sequence_started && startswith(line, "FT")) {
                parse_feature_line("     " + line.substr(5));
                continue;
            }
            if (!sequence_started && startswith(line, "FEATURES")) {
                in_features = true;
                continue;
            }
            if (!sequence_started &&
                    (startswith(line, "ORIGIN") || startswith(line, "SQ"))) {
                sequence_started = true;
                in_features = false;
            } else if (line == "//") {
                break;
            } else if (!sequence_started && in_features) {
                parse_feature_line(line);
            } else if (sequence_started) {
                sequence_ += line;
            }
        }
        if (!sequence_started)
            logger("Error: no sequence section in GenBank/EMBL record");
        normalize(sequence_);
    } else {
        logger("Error: unsupported input sequence format.");
    }
    return true;
}

const std::string &scaffold::name() const { return name_; }

const std::vector<feature> &scaffold::embedded_features() const {
    return embedded_features_;
}
