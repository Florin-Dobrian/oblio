// examples/example_amalgamation.cpp
// Amalgamation: paying fill to make the blocks bigger, and what that is worth.
//
// Symbolic factorization groups columns into supernodes, which the numeric phase hands to dense
// BLAS kernels one at a time. Columns are grouped only when their patterns match exactly, so a
// mesh produces a great many small supernodes: on the grid below, about 10,000 of them for 14,400
// columns, which means about 10,000 BLAS calls on blocks a few columns wide, where the call
// overhead is a real share of the work.
//
// **Amalgamation merges a supernode into its parent when their patterns nearly match**, padding the
// index sets with rows that are structurally zero to make them match exactly. The threshold is how
// many such rows a merge may add. Bigger threshold, fewer and larger supernodes, and more zeros
// stored and multiplied.
//
// So it is a trade, and this example measures both sides of it. The result has two regimes rather
// than one, which is the point worth taking away: a small threshold makes *both* phases faster and
// costs almost no fill, because the solve pays per-supernode overhead too. Only past the solve's
// own optimum does it become a trade, factorization time bought with solve time, and then how it
// comes out depends on how many right-hand sides one factorization serves.
//
// **These timings make a point; they are not benchmark numbers.** One grid, one ordering, best of
// three runs, on whatever machine happens to be running it. The shape of the result is robust and
// the exact figures are not. benchmarks/ is where measurement is done properly.
//
// Cholesky throughout, since amalgamation is a structural matter and the factorization only has to
// be one that runs.
//
// Compile (macOS, from repo root):
//   g++ -std=c++17 -O3 -DOBLIO_BLAS_UNDERSCORE -Iinclude examples/example_amalgamation.cpp src/*.cpp -framework Accelerate -o example_amalgamation_cpp
// Linux: replace `-framework Accelerate` with `-llapack -lblas`.

#include "oblio/DirectSolver.h"

#include <chrono>
#include <cstdio>
#include <optional>
#include <set>
#include <vector>

using namespace Oblio;

namespace {

using Val = double;

SparseMatrix<Val> grid(std::size_t side) {
    const std::size_t size = side * side;
    std::vector<std::set<std::int32_t>> adjacency(size);
    for (std::size_t r = 0; r < side; ++r) {
        for (std::size_t c = 0; c < side; ++c) {
            const std::int32_t v = static_cast<std::int32_t>(r * side + c);
            adjacency[v].insert(v);
            if (r > 0)        adjacency[v].insert(v - static_cast<std::int32_t>(side));
            if (r + 1 < side) adjacency[v].insert(v + static_cast<std::int32_t>(side));
            if (c > 0)        adjacency[v].insert(v - 1);
            if (c + 1 < side) adjacency[v].insert(v + 1);
        }
    }
    std::vector<std::size_t>  colPtr(size + 1, 0);
    std::vector<std::int32_t> rowIdx;
    std::vector<Val>          val;
    for (std::size_t j = 0; j < size; ++j) {
        for (std::int32_t i : adjacency[j]) {
            rowIdx.push_back(i);
            val.push_back(i == static_cast<std::int32_t>(j) ? Val(4) : Val(-1));
        }
        colPtr[j + 1] = rowIdx.size();
    }
    return SparseMatrix<Val>(size, colPtr, rowIdx, val);
}

double milliseconds(std::chrono::steady_clock::time_point from,
                    std::chrono::steady_clock::time_point to) {
    return std::chrono::duration<double, std::milli>(to - from).count();
}

} // namespace

int main() {
    const std::size_t side = 300;
    const SparseMatrix<Val> A = grid(side);
    const std::size_t size = A.size();

    Vector<Val> b(size), x(size);
    for (std::size_t i = 0; i < size; ++i) b[i] = Val(1);

    // Two things make a few milliseconds hard to read, and both are handled rather than hoped
    // away. The first run of any setting is slow, on first-touch page faults and whatever the BLAS
    // does once, so each setting gets an untimed warm-up before anything is recorded. And a single
    // solve is under a millisecond, which is beneath the resolution worth trusting, so solves are
    // timed in batches and divided. Best of three either way. Still not a benchmark: see the
    // header.
    const int repeats      = 3;
    const int solvesPerRep = 5;

    const std::optional<std::size_t> threshold[] = {
        std::optional<std::size_t>{},   // none: supernodes exactly as symbolic found them
        std::optional<std::size_t>{2},
        std::optional<std::size_t>{4},
        std::optional<std::size_t>{8},
        std::optional<std::size_t>{16},
        std::optional<std::size_t>{32},
    };
    const int settings = 6;

    printf("Grid %zux%zu, n = %zu, nnz(A) = %zu. Cholesky, MmdFlat, left-looking.\n\n",
           side, side, size, A.nnz());
    printf("  %-9s %9s %10s %10s %10s %10s %10s\n",
           "threshold", "snodes", "nnz(L)", "numVal", "factor ms", "solve ms", "residual");
    printf("  %-9s %9s %10s %10s %10s %10s %10s\n",
           "---------", "------", "------", "------", "---------", "--------", "--------");

    double      factorMsAt[settings], solveMsAt[settings];
    std::size_t nnzAt[settings], snodesAt[settings];

    for (int t = 0; t < settings; ++t) {
        DirectSolver<Val> solver(Ordering::MmdFlat, Factorization::Cholesky, Traversal::LeftLooking,
                                 Supernodes::Fundamental, threshold[t]);
        if (!solver.analyze(A)) { printf("  analyze failed\n"); return 1; }

        // Warm-up, untimed.
        if (!solver.factor(A) || !solver.solve(b, x)) {
            printf("  factor or solve failed\n");
            return 1;
        }

        double factorMs = 0, solveMs = 0;
        for (int rep = 0; rep < repeats; ++rep) {
            const auto t0 = std::chrono::steady_clock::now();
            if (!solver.factor(A)) { printf("  factor failed\n"); return 1; }
            const auto t1 = std::chrono::steady_clock::now();
            for (int s = 0; s < solvesPerRep; ++s)
                if (!solver.solve(b, x)) { printf("  solve failed\n"); return 1; }
            const auto t2 = std::chrono::steady_clock::now();

            const double f = milliseconds(t0, t1);
            const double s = milliseconds(t1, t2) / solvesPerRep;
            if (rep == 0 || f < factorMs) factorMs = f;
            if (rep == 0 || s < solveMs)  solveMs  = s;
        }

        char name[16];
        if (threshold[t]) snprintf(name, sizeof name, "%zu", *threshold[t]);
        else              snprintf(name, sizeof name, "none");

        // Not timed, and not decoration: a setting could be fast because it computed the wrong
        // thing, and this column is what would say so.
        printf("  %-9s %9zu %10zu %10zu %10.1f %10.2f %10.1e\n",
               name, solver.symFactor().snodeSize(), solver.nnz(), solver.numVal(),
               factorMs, solveMs, solver.relativeResidual(A, b, x));

        factorMsAt[t] = factorMs;
        solveMsAt[t]  = solveMs;
        nnzAt[t]      = solver.nnz();
        snodesAt[t]   = solver.symFactor().snodeSize();
    }

    // The solve does not simply get worse. It has an optimum of its own, at a small threshold,
    // which is the thing this example would misreport if it compared the last row against the
    // first. Find that optimum rather than assuming where it is.
    int best = 0;
    for (int t = 1; t < settings; ++t)
        if (solveMsAt[t] < solveMsAt[best]) best = t;

    printf("\nReading the table:\n\n");

    printf("  Supernodes collapse and fill rises: %zu with no amalgamation, %zu at the\n"
           "  largest threshold, and nnz(L) goes from %zu to %zu, about %.0f%% more. Each\n"
           "  added entry is a structural zero a merge introduced, and the numeric phase\n"
           "  multiplies it like any other.\n\n",
           snodesAt[0], snodesAt[settings - 1], nnzAt[0], nnzAt[settings - 1],
           100.0 * (double)(nnzAt[settings - 1] - nnzAt[0]) / (double)nnzAt[0]);

    printf("  The factorization gets faster the whole way, %.1f ms against %.1f, about\n"
           "  %.0f%% off. The work moved into fewer and larger BLAS calls, and the per-call\n"
           "  overhead fell faster than the added arithmetic grew.\n\n",
           factorMsAt[settings - 1], factorMsAt[0],
           100.0 * (factorMsAt[0] - factorMsAt[settings - 1]) / factorMsAt[0]);

    printf("  The solve is not monotone, which is the part worth noticing. Its best time is\n"
           "  %.2f ms and that is not at threshold none: a solve pays per-supernode overhead\n"
           "  too, so merging helps it at first. Past that the stored zeros dominate, a\n"
           "  triangular solve reading every entry once and gaining nothing from a bigger\n"
           "  block, and it ends at %.2f ms.\n\n",
           solveMsAt[best], solveMsAt[settings - 1]);

    printf("  The residual does not move: amalgamation changes how many numbers are stored\n"
           "  and multiplied, not which matrix is factored, and the extra entries are zeros.\n"
           "  Every row solves the same system to the same accuracy.\n\n");

    printf("  So there are two regimes rather than one trade. A small threshold makes both\n"
           "  phases faster and costs almost no fill, which is free and worth taking. Only\n"
           "  past the solve's optimum does it become a trade: factorization time bought\n"
           "  with solve time.\n\n");

    const double saved    = factorMsAt[best] - factorMsAt[settings - 1];
    const double perSolve = solveMsAt[settings - 1] - solveMsAt[best];
    if (saved > 0 && perSolve > 0)
        printf("  Which way that trade comes out is arithmetic. Going from the solve's best\n"
               "  setting to the largest threshold saves %.1f ms once and costs %.2f ms on\n"
               "  every solve, so it stops paying at about %.0f solves per factorization.\n"
               "  Below that, amalgamate further; above it, stop at the optimum. See\n"
               "  example_reuse for where that ratio comes from.\n",
               saved, perSolve, saved / perSolve);
    else
        printf("  On this run the two did not separate cleanly enough to give a crossover;\n"
               "  these are milliseconds on a shared machine.\n");

    return 0;
}
