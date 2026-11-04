/*
 * @Description: Files that carry results between the steps of the GxE workflow
 * @Author: Chao Ning
 *
 * .null.bin layout (native little-endian):
 *   char[8]   magic "FGXENULL"
 *   uint32    format version
 *   uint32    flags (bit 0: no_noisebye)
 *   uint64    fingerprint
 *   int64     number of samples n, environments q, variance components c
 *   string    trait name
 *   string[n] sample IDs
 *   string[q] environment names
 *   double[c] variance components
 *   double[n] V^{-1} y
 *   double[n * q] environment matrix, column-major
 *   int64     number of non-zeros nnz of V^{-1}
 *   int32[n + 1], int32[nnz], double[nnz]  V^{-1} in compressed-column form
 * A string is an int64 length followed by its bytes.
 */
#include "null_model_io.hpp"

#include <cerrno>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <unordered_map>

#include "../utils/fatal_error.hpp"
#include "../utils/string_utils.hpp"

namespace {

constexpr char kNullMagic[8] = {'F', 'G', 'X', 'E', 'N', 'U', 'L', 'L'};
constexpr std::uint32_t kNullFormatVersion = 1;
constexpr std::uint32_t kFlagNoNoiseBye = 1u;

// A printed double read back with 17 significant digits equals the original.
constexpr int kExactDoubleDigits = 17;

class BinaryWriter {
public:
    explicit BinaryWriter(const std::string& path) : _path(path), _out(path, std::ios::binary) {
        if (!_out.is_open()) {
            fatal_error("Fail to open output file: {}", _path);
        }
    }

    void raw(const void* data, std::size_t num_bytes) {
        _out.write(static_cast<const char*>(data), static_cast<std::streamsize>(num_bytes));
    }

    template <typename T>
    void value(T item) { raw(&item, sizeof(T)); }

    void text(const std::string& item) {
        value<std::int64_t>(static_cast<std::int64_t>(item.size()));
        raw(item.data(), item.size());
    }

    void finish() {
        _out.close();
        if (_out.fail()) {
            fatal_error("Failed to write {}. Is the disk full?", _path);
        }
    }

private:
    std::string _path;
    std::ofstream _out;
};

class BinaryReader {
public:
    explicit BinaryReader(const std::string& path) : _path(path), _in(path, std::ios::binary) {
        if (!_in.is_open()) {
            fatal_error("Fail to open file: {}", _path);
        }
    }

    void raw(void* data, std::size_t num_bytes) {
        _in.read(static_cast<char*>(data), static_cast<std::streamsize>(num_bytes));
        if (static_cast<std::size_t>(_in.gcount()) != num_bytes) {
            fatal_error("{} is truncated or corrupt.", _path);
        }
    }

    template <typename T>
    T value() {
        T item{};
        raw(&item, sizeof(T));
        return item;
    }

    std::string text() {
        const std::int64_t size = value<std::int64_t>();
        if (size < 0) {
            fatal_error("{} is corrupt.", _path);
        }
        std::string item(static_cast<std::size_t>(size), '\0');
        raw(item.data(), item.size());
        return item;
    }

    bool at_end() { return _in.peek() == std::ifstream::traits_type::eof(); }

private:
    std::string _path;
    std::ifstream _in;
};

// 64-bit FNV-1a.
class Fingerprint {
public:
    void add(const void* data, std::size_t num_bytes) {
        const unsigned char* bytes = static_cast<const unsigned char*>(data);
        for (std::size_t i = 0; i < num_bytes; ++i) {
            _hash = (_hash ^ bytes[i]) * 1099511628211ULL;
        }
    }

    void add(const std::string& item) {
        const std::int64_t size = static_cast<std::int64_t>(item.size());
        add(&size, sizeof(size));
        add(item.data(), item.size());
    }

    std::uint64_t digest() const { return _hash; }

private:
    std::uint64_t _hash = 14695981039346656037ULL;
};

std::string null_model_path(const std::string& prefix) { return prefix + ".null.bin"; }
std::string gamma_path(const std::string& prefix) { return prefix + ".gamma"; }

std::string to_hex(std::uint64_t value) {
    std::ostringstream stream;
    stream << std::hex << std::setw(16) << std::setfill('0') << value;
    return stream.str();
}

double parse_double(const std::string& text, const std::string& term, const std::string& path) {
    char* end = nullptr;
    const double value = std::strtod(text.c_str(), &end);
    if (end == text.c_str() || *end != '\0') {
        fatal_error("Invalid value for {} in {}: {}", term, path, text);
    }
    return value;
}

}  // namespace

std::uint64_t compute_null_fingerprint(const std::vector<std::string>& sample_ids,
                                       const std::vector<std::string>& env_names,
                                       const Eigen::VectorXd& varcom,
                                       const Eigen::VectorXd& vinv_y,
                                       const Eigen::MatrixXd& env_mat) {
    Fingerprint fingerprint;
    for (const std::string& id : sample_ids) fingerprint.add(id);
    for (const std::string& name : env_names) fingerprint.add(name);
    fingerprint.add(varcom.data(), sizeof(double) * static_cast<std::size_t>(varcom.size()));
    fingerprint.add(vinv_y.data(), sizeof(double) * static_cast<std::size_t>(vinv_y.size()));
    fingerprint.add(env_mat.data(), sizeof(double) * static_cast<std::size_t>(env_mat.size()));
    return fingerprint.digest();
}

void write_null_model(const std::string& prefix, const NullModel& model) {
    const std::int64_t num_id = static_cast<std::int64_t>(model.sample_ids.size());
    const std::int64_t num_env = static_cast<std::int64_t>(model.env_names.size());
    if (model.vinv_y.size() != num_id || model.env_mat.rows() != num_id || model.env_mat.cols() != num_env ||
            model.vinv.rows() != num_id || model.vinv.cols() != num_id || !model.vinv.isCompressed()) {
        fatal_error("Inconsistent null model; refusing to write {}", null_model_path(prefix));
    }

    // Write to a temporary name first, so that an interrupted run cannot leave a
    // partial file under the final name.
    const std::string final_path = null_model_path(prefix);
    const std::string temporary_path = final_path + ".tmp";
    BinaryWriter writer(temporary_path);

    writer.raw(kNullMagic, sizeof(kNullMagic));
    writer.value<std::uint32_t>(kNullFormatVersion);
    writer.value<std::uint32_t>(model.no_noisebye ? kFlagNoNoiseBye : 0u);
    writer.value<std::uint64_t>(model.fingerprint);
    writer.value<std::int64_t>(num_id);
    writer.value<std::int64_t>(num_env);
    writer.value<std::int64_t>(model.varcom.size());

    writer.text(model.trait_name);
    for (const std::string& id : model.sample_ids) writer.text(id);
    for (const std::string& name : model.env_names) writer.text(name);

    writer.raw(model.varcom.data(), sizeof(double) * static_cast<std::size_t>(model.varcom.size()));
    writer.raw(model.vinv_y.data(), sizeof(double) * static_cast<std::size_t>(num_id));
    writer.raw(model.env_mat.data(), sizeof(double) * static_cast<std::size_t>(num_id * num_env));

    const std::int64_t nnz = model.vinv.nonZeros();
    writer.value<std::int64_t>(nnz);
    writer.raw(model.vinv.outerIndexPtr(), sizeof(int) * static_cast<std::size_t>(num_id + 1));
    writer.raw(model.vinv.innerIndexPtr(), sizeof(int) * static_cast<std::size_t>(nnz));
    writer.raw(model.vinv.valuePtr(), sizeof(double) * static_cast<std::size_t>(nnz));
    writer.finish();

    if (std::rename(temporary_path.c_str(), final_path.c_str()) != 0) {
        fatal_error("Failed to rename {} to {}: {}", temporary_path, final_path, std::strerror(errno));
    }
    spdlog::info("Null model written to {} ({} samples, {} environments, {} non-zeros in V^-1)",
                 final_path, num_id, num_env, nnz);
}

NullModel read_null_model(const std::string& prefix, bool load_vinv) {
    const std::string path = null_model_path(prefix);
    spdlog::info("Reading null model: {}", path);
    BinaryReader reader(path);

    char magic[sizeof(kNullMagic)];
    reader.raw(magic, sizeof(magic));
    if (std::memcmp(magic, kNullMagic, sizeof(magic)) != 0) {
        fatal_error("{} is not a fastGxE null-model file.", path);
    }
    const std::uint32_t format_version = reader.value<std::uint32_t>();
    if (format_version != kNullFormatVersion) {
        fatal_error("{} has format version {}, but this program reads version {}.",
                    path, format_version, kNullFormatVersion);
    }

    NullModel model;
    const std::uint32_t flags = reader.value<std::uint32_t>();
    model.no_noisebye = (flags & kFlagNoNoiseBye) != 0;
    model.fingerprint = reader.value<std::uint64_t>();
    const std::int64_t num_id = reader.value<std::int64_t>();
    const std::int64_t num_env = reader.value<std::int64_t>();
    const std::int64_t num_varcom = reader.value<std::int64_t>();
    if (num_id <= 0 || num_env <= 0 || num_varcom <= 0) {
        fatal_error("{} is corrupt.", path);
    }

    model.trait_name = reader.text();
    model.sample_ids.reserve(static_cast<std::size_t>(num_id));
    for (std::int64_t i = 0; i < num_id; ++i) model.sample_ids.push_back(reader.text());
    model.env_names.reserve(static_cast<std::size_t>(num_env));
    for (std::int64_t i = 0; i < num_env; ++i) model.env_names.push_back(reader.text());

    model.varcom.resize(num_varcom);
    reader.raw(model.varcom.data(), sizeof(double) * static_cast<std::size_t>(num_varcom));
    model.vinv_y.resize(num_id);
    reader.raw(model.vinv_y.data(), sizeof(double) * static_cast<std::size_t>(num_id));
    model.env_mat.resize(num_id, num_env);
    reader.raw(model.env_mat.data(), sizeof(double) * static_cast<std::size_t>(num_id * num_env));

    if (load_vinv) {
        const std::int64_t nnz = reader.value<std::int64_t>();
        if (nnz < 0 || nnz > num_id * num_id) {
            fatal_error("{} is corrupt.", path);
        }
        model.vinv.resize(num_id, num_id);
        model.vinv.resizeNonZeros(nnz);
        reader.raw(model.vinv.outerIndexPtr(), sizeof(int) * static_cast<std::size_t>(num_id + 1));
        reader.raw(model.vinv.innerIndexPtr(), sizeof(int) * static_cast<std::size_t>(nnz));
        reader.raw(model.vinv.valuePtr(), sizeof(double) * static_cast<std::size_t>(nnz));
        if (!reader.at_end()) {
            fatal_error("{} has unexpected data after V^-1.", path);
        }
    }

    const std::uint64_t expected = compute_null_fingerprint(
        model.sample_ids, model.env_names, model.varcom, model.vinv_y, model.env_mat);
    if (expected != model.fingerprint) {
        fatal_error("{} failed its integrity check.", path);
    }

    spdlog::info("Null model: trait {}, {} samples, {} environments", model.trait_name, num_id, num_env);
    return model;
}

void write_gamma_file(const std::string& prefix, const GxeGamma& gamma,
                      const std::vector<std::string>& env_names) {
    const std::string path = gamma_path(prefix);
    std::ofstream fout(path);
    if (!fout.is_open()) {
        fatal_error("Fail to open output file: {}", path);
    }

    fout << std::setprecision(kExactDoubleDigits);
    fout << "term\tvalue\n";
    fout << "gamma_main\t" << gamma.main << "\n";
    for (std::size_t j = 0; j < env_names.size(); ++j) {
        fout << "gamma_noMain_" << env_names[j] << "\t" << gamma.no_main(j) << "\n";
    }
    for (std::size_t j = 0; j < env_names.size(); ++j) {
        fout << "gamma_withMain_" << env_names[j] << "\t" << gamma.with_main[j](1, 1) << "\n";
    }
    for (std::size_t j = 0; j < env_names.size(); ++j) {
        fout << "gamma_withMain00_" << env_names[j] << "\t" << gamma.with_main[j](0, 0) << "\n";
        fout << "gamma_withMain01_" << env_names[j] << "\t" << gamma.with_main[j](0, 1) << "\n";
        fout << "gamma_withMain10_" << env_names[j] << "\t" << gamma.with_main[j](1, 0) << "\n";
    }
    fout << "n_random_used\t" << gamma.num_random_used << "\n";
    fout << "n_random_requested\t" << gamma.num_random_requested << "\n";
    fout << "seed\t" << gamma.seed << "\n";
    fout << "null_fingerprint\t" << to_hex(gamma.null_fingerprint) << "\n";
    fout.close();
    if (fout.fail()) {
        fatal_error("Failed to write {}. Is the disk full?", path);
    }
}

GxeGamma read_gamma_file(const std::string& prefix, const std::vector<std::string>& env_names) {
    const std::string path = gamma_path(prefix);
    spdlog::info("Reading gamma factors: {}", path);
    std::ifstream fin(path);
    if (!fin.is_open()) {
        fatal_error("Fail to open file: {}", path);
    }

    std::unordered_map<std::string, std::string> values;
    std::string one_line;
    while (std::getline(fin, one_line)) {
        process_line(one_line);
        if (one_line.empty()) continue;
        const std::vector<std::string> fields = split_string(one_line);
        if (fields.size() == 2) values[fields[0]] = fields[1];
    }

    const auto text_of = [&](const std::string& term) -> const std::string& {
        const auto entry = values.find(term);
        if (entry == values.end()) {
            fatal_error("{} has no row for {}. Was it written by step 2 of this null model?", path, term);
        }
        return entry->second;
    };
    const auto number_of = [&](const std::string& term) {
        return parse_double(text_of(term), term, path);
    };

    GxeGamma gamma;
    gamma.main = number_of("gamma_main");
    gamma.no_main.resize(static_cast<Eigen::Index>(env_names.size()));
    gamma.with_main.reserve(env_names.size());
    for (std::size_t j = 0; j < env_names.size(); ++j) {
        gamma.no_main(j) = number_of("gamma_noMain_" + env_names[j]);
        Eigen::MatrixXd with_main(2, 2);
        with_main(0, 0) = number_of("gamma_withMain00_" + env_names[j]);
        with_main(0, 1) = number_of("gamma_withMain01_" + env_names[j]);
        with_main(1, 0) = number_of("gamma_withMain10_" + env_names[j]);
        with_main(1, 1) = number_of("gamma_withMain_" + env_names[j]);
        gamma.with_main.push_back(with_main);
    }
    gamma.num_random_used = static_cast<std::int64_t>(number_of("n_random_used"));
    gamma.num_random_requested = static_cast<std::int64_t>(number_of("n_random_requested"));
    gamma.seed = std::strtoull(text_of("seed").c_str(), nullptr, 10);
    gamma.null_fingerprint = std::strtoull(text_of("null_fingerprint").c_str(), nullptr, 16);
    return gamma;
}
