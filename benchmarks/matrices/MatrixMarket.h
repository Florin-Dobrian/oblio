#pragma once

// MatrixMarket.h -- read a Matrix Market file into a SparseMatrix.
//
// WHAT THIS IS FOR. The tree's matrices have all been generated until now, two families of grid
// plus a handful of hand-built examples, and every claim drawn on them is a claim about them.
// This is how real ones get in. The SuiteSparse Matrix Collection publishes Matrix Market as one
// of its three formats, and it is the one worth reading: text, one page of specification, and the
// header states the field type and the symmetry outright, which is what decides whether we can
// take the file at all.
//
// WHAT IS ACCEPTED, and everything else is REFUSED with the reason rather than worked around:
//
//   object      matrix          the format also covers dense arrays under `array`, not this
//   format      coordinate      the sparse spelling
//   field       real, integer   both become double here, the distinction being how the file
//                               spells its numbers rather than what they mean. Gset's weights
//                               are +1 and -1 and Trefethen's entries are primes, so several
//                               matrices in the collection are legitimately `integer`
//   symmetry    symmetric       A = A^T, and the file stores ONE TRIANGLE
//
// `complex` and `hermitian` are refused because this reader produces a SparseMatrix<double>,
// not because Oblio cannot factor them: it can, and complex Hermitian dynamic LDL is the one
// part of the numeric code with no 0.9 reference behind it, so real complex input is evidence
// nothing else can give it. That is a later step, and it is a second instantiation rather than
// a parse change. `general` is refused because a matrix stored in full may still be unsymmetric,
// and nothing here checks.
//
// `pattern` FILES ARE ACCEPTED ONLY ON REQUEST, through `acceptPattern`. They carry no values at
// all, two indices per line, so every entry becomes 1.0 and the diagonal stays the 0.0 the
// conversion inserts. That is useless for a residual and exactly right for a caller that reads
// only the structure, which is what an ordering does: `mmdmatrices.cpp` in experiments/ordering
// compares two orderings of the same pattern and never looks at a value. The default is false so
// that a caller who needs numbers cannot get a matrix of ones by accident.
//
// THIS READER MIRRORS. A `symmetric` file holds the lower triangle and Oblio stores both, so an
// off-diagonal entry is emitted twice, once each way. The diagonal is emitted once: fromTriplets
// ACCUMULATES on collision, which is what an assembly loop needs and what would double every
// diagonal value if a diagonal entry were mirrored onto itself.
//
// AND IT CHECKS THE TRIANGLE. The format says a `symmetric` file stores the lower triangle, and
// the collection's files do, but nothing in a file enforces it. An entry above the diagonal in a
// file we then mirror produces a duplicate that fromTriplets silently SUMS, so a malformed file
// would turn into a wrong number rather than into an error. The check is one comparison per
// entry and it is the only validation anywhere in this tree today.
//
// WHAT IS NOT CHECKED, and it is worth knowing rather than fixing here: whether the values make
// a well posed system. A graph's adjacency matrix carrying weights is a legitimate symmetric
// matrix and the factorizations will do something sensible with it, but a residual measured on
// one is a statement about our arithmetic rather than about anybody's problem.

#include "oblio/SparseMatrix.h"

#include <cctype>
#include <cstdlib>
#include <fstream>
#include <map>
#include <string>
#include <utility>
#include <vector>

namespace MatrixBenchmark {

// One nonzero as the file gives it: row, column, value, zero based by the time it is here.
struct Triplet {
    std::int32_t row;
    std::int32_t col;
    double       val;
};

// What a read produced. A refusal carries its reason, since the driver reports rather than
// stops: a file we cannot take is information about the collection, not an error in the run.
struct ReadResult {
    bool                        ok = false;
    std::string                 reason;       // empty when ok
    Oblio::SparseMatrix<double> matrix;
    std::size_t                 fileEntries = 0;   // entries as the file stored them, one triangle
};

namespace detail {

// Coordinate form to CSC. Taken from examples/example_matrix.cpp, whose comment calls it "the
// conversion worth stealing" and states the three preconditions a SparseMatrix has: full storage,
// a structurally present diagonal, and row indices sorted ascending with no duplicates. A map
// keyed by (column, row) gives the sort and the deduplication together, and the diagonal is
// inserted where nothing landed on it, which is what the symbolic factorization needs and is
// ordinary input for LDL even when the value is zero.
//
// It is the readable form rather than the fast one, one map node per nonzero, and that is the
// example's own stated choice. At the sizes this folder starts with it costs a fraction of a
// second; a counting sort is what it becomes if the matrices grow.
inline Oblio::SparseMatrix<double> fromTriplets(std::size_t size,
                                                const std::vector<Triplet>& triplets) {
    std::map<std::pair<std::int32_t, std::int32_t>, double> entries;

    for (const Triplet& t : triplets)
        entries[{t.col, t.row}] += t.val;

    for (std::int32_t j = 0; j < static_cast<std::int32_t>(size); ++j)
        entries.emplace(std::make_pair(j, j), 0.0);

    std::vector<std::size_t>  colPtr(size + 1, 0);
    std::vector<std::int32_t> rowIdx;
    std::vector<double>       val;
    rowIdx.reserve(entries.size());
    val.reserve(entries.size());

    auto entry = entries.begin();
    for (std::size_t j = 0; j < size; ++j) {
        while (entry != entries.end() && entry->first.first == static_cast<std::int32_t>(j)) {
            rowIdx.push_back(entry->first.second);
            val.push_back(entry->second);
            ++entry;
        }
        colPtr[j + 1] = rowIdx.size();
    }

    return Oblio::SparseMatrix<double>(size, colPtr, rowIdx, val);
}

inline std::string lowered(const std::string& text) {
    std::string result = text;
    for (char& c : result)
        c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
    return result;
}

// The banner, split on whitespace. Its five words are fixed by the format:
// %%MatrixMarket <object> <format> <field> <symmetry>
inline std::vector<std::string> bannerWords(const std::string& line) {
    std::vector<std::string> words;
    std::size_t at = 0;
    while (at < line.size()) {
        while (at < line.size() && std::isspace(static_cast<unsigned char>(line[at])))
            ++at;
        const std::size_t start = at;
        while (at < line.size() && !std::isspace(static_cast<unsigned char>(line[at])))
            ++at;
        if (at > start)
            words.push_back(lowered(line.substr(start, at - start)));
    }
    return words;
}

inline bool blankOrComment(const std::string& line) {
    for (char c : line) {
        if (std::isspace(static_cast<unsigned char>(c)))
            continue;
        return c == '%';
    }
    return true;
}

inline ReadResult refuse(const std::string& reason) {
    ReadResult result;
    result.ok = false;
    result.reason = reason;
    return result;
}

} // namespace detail

// Read one file. Never throws for a file we decline to take; a refusal comes back in the result.
// `acceptPattern` takes files that carry structure and no values; see the banner note above.
inline ReadResult readMatrixMarket(const std::string& path, bool acceptPattern = false) {
    std::ifstream file(path);
    if (!file)
        return detail::refuse("cannot open");

    std::string line;
    if (!std::getline(file, line))
        return detail::refuse("empty file");

    const std::vector<std::string> words = detail::bannerWords(line);
    if (words.size() < 5 || words[0] != "%%matrixmarket")
        return detail::refuse("no Matrix Market banner");

    const std::string& object = words[1];
    const std::string& format = words[2];
    const std::string& field = words[3];
    const std::string& symmetry = words[4];

    if (object != "matrix")
        return detail::refuse("object is " + object);
    if (format != "coordinate")
        return detail::refuse("format is " + format);
    const bool pattern = (field == "pattern");
    if (field != "real" && field != "integer" && !(pattern && acceptPattern))
        return detail::refuse("field is " + field);
    if (symmetry != "symmetric")
        return detail::refuse("symmetry is " + symmetry);

    // The size line, past any comments.
    std::size_t lineNumber = 1;
    while (std::getline(file, line)) {
        ++lineNumber;
        if (!detail::blankOrComment(line))
            break;
    }
    if (file.eof() && detail::blankOrComment(line))
        return detail::refuse("no size line");

    long long rows = 0;
    long long cols = 0;
    long long entries = 0;
    {
        const char* cursor = line.c_str();
        char* end = nullptr;
        rows = std::strtoll(cursor, &end, 10);
        cols = std::strtoll(end, &end, 10);
        entries = std::strtoll(end, &end, 10);
        if (end == cursor || rows <= 0 || cols <= 0 || entries < 0)
            return detail::refuse("malformed size line");
    }
    if (rows != cols)
        return detail::refuse("not square");

    const std::size_t size = static_cast<std::size_t>(rows);

    std::vector<Triplet> triplets;
    triplets.reserve(static_cast<std::size_t>(entries) * 2);

    long long read = 0;
    while (read < entries && std::getline(file, line)) {
        ++lineNumber;
        if (detail::blankOrComment(line))
            continue;

        const char* cursor = line.c_str();
        char* end = nullptr;
        const long long row = std::strtoll(cursor, &end, 10);
        const char* afterRow = end;
        const long long col = std::strtoll(afterRow, &end, 10);
        if (end == afterRow)
            return detail::refuse("malformed entry at line " + std::to_string(lineNumber));
        double val = 1.0;                     // a pattern file's every entry; see the banner note
        if (!pattern) {
            const char* afterCol = end;
            val = std::strtod(afterCol, &end);
            if (end == afterCol)
                return detail::refuse("no value at line " + std::to_string(lineNumber));
        }

        if (row < 1 || row > rows || col < 1 || col > cols)
            return detail::refuse("index out of range at line " + std::to_string(lineNumber));

        // The format puts a symmetric matrix's entries in the LOWER triangle. Nothing in the file
        // enforces it, and mirroring an upper entry would silently double a value.
        if (col > row)
            return detail::refuse("upper-triangle entry at line " + std::to_string(lineNumber));

        const std::int32_t r = static_cast<std::int32_t>(row - 1);
        const std::int32_t c = static_cast<std::int32_t>(col - 1);

        triplets.push_back({r, c, val});
        if (r != c)
            triplets.push_back({c, r, val});   // the mirror; the diagonal is emitted once

        ++read;
    }

    if (read != entries)
        return detail::refuse("file ends after " + std::to_string(read) + " of "
                              + std::to_string(entries) + " entries");

    ReadResult result;
    result.ok = true;
    result.fileEntries = static_cast<std::size_t>(entries);
    result.matrix = detail::fromTriplets(size, triplets);
    return result;
}

} // namespace MatrixBenchmark
