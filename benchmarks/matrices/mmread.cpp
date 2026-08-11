// mmread.cpp -- read every matrix named on the command line and say what came back.
//
// The reader before the pipeline, deliberately. Every surprise real input holds lands here
// first, and a program that only parses says what it found without a factorization in the way.
// It is also the first thing that will run over a whole downloaded set, so it doubles as the
// inventory: what is usable, what is refused and why.
//
//   ./mmread_cpp ../../data/*/*.mtx
//
// One line per file, sorted by whatever the shell hands over. A refusal prints its reason and
// the run continues, which is this folder's rule: a file we cannot take is information about
// the collection rather than an error in the run.
//
// Build, from this directory (no BLAS: nothing here factors anything):
//   g++ -std=c++17 -O3 -DNDEBUG -Wall -Wextra -I../../include
//       mmread.cpp ../../src/SparseMatrix.cpp ../../src/Types.cpp -o mmread_cpp

#include "MatrixMarket.h"

#include <cstdio>
#include <exception>
#include <string>

namespace {

// The file name without its directory, so a long path does not push the numbers off the line.
std::string shortName(const std::string& path) {
    const std::size_t slash = path.find_last_of('/');
    if (slash == std::string::npos)
        return path;
    const std::size_t before = path.find_last_of('/', slash - 1);
    if (before == std::string::npos)
        return path;
    return path.substr(before + 1);
}

// Structural zeros on the diagonal are inserted by the conversion, so a matrix whose file gave
// no diagonal at all is visible here and nowhere else. Graph patterns are the usual case.
std::size_t zeroDiagonal(const Oblio::SparseMatrix<double>& A) {
    std::size_t count = 0;
    const std::vector<std::size_t>&  colPtr = A.colPtr();
    const std::vector<std::int32_t>& rowIdx = A.rowIdx();
    const std::vector<double>&       val = A.val();

    for (std::size_t j = 0; j < A.size(); ++j)
        for (std::size_t cp = colPtr[j]; cp < colPtr[j + 1]; ++cp)
            if (rowIdx[cp] == static_cast<std::int32_t>(j) && val[cp] == 0.0)
                ++count;

    return count;
}

} // namespace

int main(int argc, char** argv) {
    if (argc < 2) {
        std::fprintf(stderr, "usage: %s <file.mtx> ...\n", argv[0]);
        return 1;
    }

    int taken = 0;
    int refused = 0;

    std::printf("%-44s %9s %11s %11s %9s\n",
                "matrix", "n", "file nnz", "nnz(A)", "zero diag");

    for (int i = 1; i < argc; ++i) {
        const std::string path = argv[i];
        std::fflush(stdout);   // the name is printed with the result, so flush before a crash

        MatrixBenchmark::ReadResult result;
        try {
            result = MatrixBenchmark::readMatrixMarket(path);
        } catch (const std::exception& error) {
            std::printf("%-44s SKIP  %s\n", shortName(path).c_str(), error.what());
            ++refused;
            continue;
        }

        if (!result.ok) {
            std::printf("%-44s SKIP  %s\n", shortName(path).c_str(), result.reason.c_str());
            ++refused;
            continue;
        }

        std::printf("%-44s %9zu %11zu %11zu %9zu\n",
                    shortName(path).c_str(),
                    result.matrix.size(),
                    result.fileEntries,
                    result.matrix.nnz(),
                    zeroDiagonal(result.matrix));
        ++taken;
    }

    std::printf("%d read, %d refused\n", taken, refused);
    return 0;
}
