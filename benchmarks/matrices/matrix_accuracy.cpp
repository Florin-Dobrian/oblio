// matrix_accuracy.cpp -- solve every matrix named on the command line and report how good the
// answer is: the backward error, the residual, and what each factorization had to do to get
// there. The report drawn from it is ACCURACY.md beside this file.
//
// The first pass over real input. Ordering is MmdCompacted and the traversal is left-looking; both
// are held fixed on purpose, because the question here is whether the pipeline computes correctly
// on matrices nobody generated, not which ordering is best. The other axes come after this one has
// been read.
//
// NOT THE TREE'S DEFAULT, AND THAT IS DELIBERATE, 2026-08-21. The default is `AmdCompacted` and
// this pass ran under mmd from the beginning, so holding mmd keeps every published figure in
// ACCURACY.md comparable: `MmdCompacted` returns `MmdFlat`'s permutation exactly, which is
// genmmd's, so switching between those two cannot move a number here.
//
// Moving to `AmdCompacted` was tried and it COST A MATRIX. `Oberwolfach/LFAT5000` is killed by the
// OOM killer under amd where mmd factors it, at a 232-fold fill increase from delayed pivots, and
// amd's PREDICTED fill on that matrix is lower than mmd's. So the loss is the delayed-pivot
// cascade rather than the ordering's fill, and it is an open question rather than a reason to
// prefer mmd generally. See docs/NEXT.md.
//
//   ./matrix_accuracy_cpp ../../data/*/*.mtx
//
// SKIP AND REPORT, everywhere. A file we cannot read, an analysis that refuses, a factorization
// that refuses: each prints its reason in its own cell and the run continues. Over a set nobody
// has looked at yet, a refusal is information about the collection rather than an error in the
// run, and one run that always completes is worth more than one that stops at the first surprise.
//
// THREE FACTORIZATIONS, and the choice of which three is not obvious. Over the reals LDL^T and
// LDL^H are the same computation, the conjugate being the identity there, so running both would
// repeat one column rather than cover a new case. What actually separates on real input is
// static against dynamic pivoting:
//
//   Cholesky      succeeds only on a positive definite matrix, so it is a detector as much as a
//                 check. A quarter of this set has a structurally absent diagonal and it must
//                 refuse on all of those
//   StaticLDLT    runs on anything, and PERTURBS a pivot too small to divide by rather than
//                 failing. A poor residual here on an indefinite matrix is the documented
//                 behavior of the method, which is why the perturbation count sits beside it
//   DynamicLDLT   chooses its pivots as it goes, taking 2x2 blocks and delaying what it cannot
//                 use. This is the one expected to give a good residual on an indefinite matrix,
//                 and the delayed column count is what it paid for it
//
// Adding the two H variants is one line each if the consistency check is ever wanted; they should
// reproduce their T counterparts exactly on real input.
//
// THE RIGHT-HAND SIDE IS ALL ONES, which is a choice worth stating because it interacts with the
// paragraph below. The obvious alternative is to take a known solution and form b = A x, which
// puts b in the range of A by construction and so is solvable whatever A's rank; it also makes
// the forward error available. It is not the default here because singular matrices are labeled
// and skipped rather than solved, which removes the case that needs it. If they ever come back
// into the table, this is the line that has to change with them, and the README says what to
// change and what to watch for, under "What to change when singular matrices come back".
//
// SINGULAR MATRICES ARE LABELED, NOT HIDDEN AND NOT SOLVED. The test is the STRUCTURAL RANK, the
// size of a maximum matching between columns and rows: below n, every term of the determinant
// expansion vanishes and the matrix is singular for almost any values, established from the
// pattern with nothing factored. Those rows print their rank in place of the residual columns.
// What it still cannot see is a matrix that is numerically singular while structurally full,
// which nothing short of factoring will. See Matching.h and the README.
//
// WHAT THE RESIDUAL DOES AND DOES NOT SAY. ||A x - b|| / ||b|| is bounded by roughly the machine
// epsilon times ||A|| ||x|| / ||b|| for a stable solve, and that factor reaches the condition
// number when x is large relative to b. Real matrices span a far wider range of conditioning
// than the grids this tree has measured on, so a large figure in these columns is not by itself
// a defect, and the scale-free measure is the backward error rather than this. This driver
// reports and does not judge: there are no thresholds here and nothing fails.
//
// Build: `make`, or see the Makefile for the compiler line.

#include "MatrixMarket.h"
#include "Matching.h"

#include "oblio/DirectSolver.h"
#include "oblio/MultiplyEngine.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>

#include <sys/resource.h>
#include <sys/wait.h>
#include <unistd.h>
#include <exception>
#include <string>
#include <vector>

using namespace Oblio;

namespace {

// A numeric option's value, or -1 with a complaint on stderr.
//
// Parsed with strtod rather than strtoull so that both 200000000 and 2e8 work, and CHECKED, which
// is the point: an earlier version took strtoull, silently read 2e8 as 2, and skipped every matrix
// in the run without saying why. A wrong number that announces itself costs a second; one that
// does not costs an afternoon.
double optionValue(const std::string& arg, std::size_t offset) {
    const char* start = arg.c_str() + offset;
    char*       end = nullptr;
    const double value = std::strtod(start, &end);

    if (end == start || *end != '\0') {
        std::fprintf(stderr, "%s: not a number\n", arg.c_str());
        return -1;
    }
    if (value <= 0) {
        std::fprintf(stderr, "%s: must be positive\n", arg.c_str());
        return -1;
    }
    return value;
}

// The group and file name, which is the collection's own identifier for a matrix. A bare name
// is not unique across groups, and the full path pushes the numbers off the line.
std::string shortName(const std::string& path) {
    const std::size_t slash = path.find_last_of('/');
    if (slash == std::string::npos)
        return path;
    const std::size_t before = path.find_last_of('/', slash - 1);
    if (before == std::string::npos)
        return path;
    return path.substr(before + 1);
}

// nnz(L) from the symbolic factor, a supernode's own triangle plus its update rows. The same
// computation the other two benchmark folders make, duplicated rather than shared on the same
// grounds: a benchmark should stand alone.
std::size_t fill(const SymFactor& sf) {
    std::size_t nnz = 0;
    for (std::int32_t kk = 0; kk < static_cast<std::int32_t>(sf.snodeSize()); ++kk) {
        const std::size_t f = sf.frontSize(kk);
        const std::size_t u = sf.updateSize(kk);
        nnz += f * (f + 1) / 2 + f * u;
    }
    return nnz;
}

// Columns holding no nonzero value at all. Each is an exactly zero row and column, so each is a
// zero eigenvalue. Kept beside the structural rank because it names the COMMON case in this
// collection, the isolated vertex, and because it is what a reader can check by eye against the
// reader's own zero-diagonal column. The rank subsumes it: an empty column can never be matched.
std::size_t emptyColumns(const SparseMatrix<double>& A) {
    std::size_t count = 0;

    for (std::size_t j = 0; j < A.size(); ++j) {
        bool nonzero = false;
        for (std::size_t cp = A.colPtr()[j]; cp < A.colPtr()[j + 1]; ++cp)
            if (A.val()[cp] != 0.0)
                nonzero = true;
        if (!nonzero)
            ++count;
    }

    return count;
}

// The induced infinity norm of A, the largest absolute row sum. Computed from the columns, which
// is the same thing here: A is symmetric and stored fully.
double infNorm(const SparseMatrix<double>& A) {
    std::vector<double> rowSum(A.size(), 0.0);

    for (std::size_t j = 0; j < A.size(); ++j)
        for (std::size_t cp = A.colPtr()[j]; cp < A.colPtr()[j + 1]; ++cp)
            rowSum[static_cast<std::size_t>(A.rowIdx()[cp])] += std::abs(A.val()[cp]);

    double norm = 0.0;
    for (double sum : rowSum)
        norm = std::max(norm, sum);
    return norm;
}

double infNorm(const Vector<double>& v) {
    double norm = 0.0;
    for (std::size_t i = 0; i < v.size(); ++i)
        norm = std::max(norm, std::abs(v[i]));
    return norm;
}

// A matrix we solved, kept so the set can be classified at the end.
struct Solved {
    std::string name;
    std::size_t size = 0;
    bool        hasInertia = false;
    Inertia     inertia;
    bool        choleskyWorked = false;
    std::size_t perturbed = 0;

    // The best of the three on EACH MEASURE SEPARATELY, which may be two different
    // factorizations, and that is deliberate. They answer different questions: the smallest
    // residual says whether ANY factorization got the residual down, the smallest backward error
    // says whether ANY was backward stable, and neither answer is improved by insisting one
    // factorization supply both.
    //
    // Taking them as a pair went wrong twice, in opposite directions. Choosing by residual, on
    // GHS_indef/sit100, static LDL's 2.6e+03 edged dynamic's 3.2e+03 and dragged its backward
    // error of 1.3e-02 along, against dynamic's 9.6e-15. Choosing by backward error, on
    // HB/plat1919, Cholesky's 7.3e-17 edged dynamic's 9.9e-17 and dragged its residual of 1.0e+03
    // along, against dynamic's 1.5e+02. Both times a rounding-level difference in one column
    // decided the other one. Decoupling them removes the tie-break entirely.
    double bestResidual = 0;
    double bestBackward = 0;
};

// A matrix we declined to solve, kept so the set can be listed together at the end.
struct Labeled {
    std::string name;
    std::size_t size = 0;
    std::size_t rank = 0;
    std::size_t empty = 0;
};

// One factorization's outcome. `note` carries what a number cannot: a refusal, or what the
// method had to do to get its answer.
struct Outcome {
    bool        solved = false;
    double      backward = 0;   // ||r|| / (||A|| ||x|| + ||b||), which conditioning does not inflate
    double      residual = 0;   // ||r|| / ||b||, which it does
    bool        hasInertia = false;
    Inertia     inertia;
    std::size_t perturbed = 0;  // static LDL only: pivots too small to divide by
    std::size_t actualFill = 0; // nnz(L) as factored, which delays can push above the prediction
    std::string note;
};

// Factor and solve with the analysis already done, so the three share one ordering and one
// symbolic factorization, which is also the property being relied on: the analysis depends on
// the pattern and not on the method.
Outcome run(DirectSolver<double>& solver, Factorization factorization,
            const SparseMatrix<double>& A, const Vector<double>& b, double aNorm) {
    Outcome outcome;
    solver.setFactorization(factorization);

    Vector<double> x(A.size());

    try {
        if (!solver.factor(A)) {
            outcome.note = "refused";
            return outcome;
        }
        if (!solver.solve(b, x)) {
            outcome.note = "no solve";
            return outcome;
        }

        // Both measures come off ONE residual vector and one set of norms, so the pair cannot
        // disagree about anything but its denominator, which is the whole point of showing both.
        Vector<double>       r(A.size());
        const MultiplyEngine multiply;
        multiply.residual(A, x, b, r);

        const double rNorm = infNorm(r);
        const double xNorm = infNorm(x);
        const double bNorm = infNorm(b);

        outcome.backward = rNorm / (aNorm * xNorm + bNorm);
        outcome.residual = (bNorm > 0.0) ? rNorm / bNorm : 0.0;
    } catch (const std::exception& error) {
        outcome.note = error.what();
        return outcome;
    }

    outcome.solved = true;

    if (factorization == Factorization::StaticLDLT) {
        outcome.perturbed = solver.numPerturbations();
        if (outcome.perturbed > 0)
            outcome.note = std::to_string(outcome.perturbed) + "p";
    } else if (factorization == Factorization::DynamicLDLT) {
        const std::size_t delayed = solver.numDelayedColumns();
        if (delayed > 0)
            outcome.note = std::to_string(delayed) + "d";
    }

    // The inertia of A, read off the signs of D: A = L D L^H is a congruence and a congruence
    // preserves those signs, which is Sylvester's law, so counting them in D counts them in A
    // without forming an eigenvalue. Recorded for every factorization and PRINTED from the
    // dynamic one only, because it is a property of A rather than of the method and the three
    // should agree. They do not always, and the exception is documented: a static LDL that
    // perturbed reports the inertia of the matrix it FACTORED, so on eurqsa it reads 83 zero
    // eigenvalues where dynamic reads 8, and on saylr3 it reads none where dynamic reads 2.
    outcome.hasInertia = solver.inertia(outcome.inertia);
    outcome.actualFill = solver.nnz();

    return outcome;
}

// A cell is the pair: backward error first, because it is the verdict, then the relative
// residual, which is the same numerator over a denominator that conditioning inflates. The note
// is a perturbation count (p) or a delayed-column count (d).
void printCell(const Outcome& outcome, bool /*last*/) {
    if (outcome.solved)
        std::printf(" %8.1e %8.1e", outcome.backward, outcome.residual);
    else
        std::printf(" %17s", "-");

    std::printf(" %-10s", outcome.note.c_str());
}

// One matrix, start to finish, run in a CHILD PROCESS so that nothing it does can take the
// run down. Everything it learns goes back over `out` as ONE LINE, since a child cannot
// append to the parent's vectors:
//
//   R  solved      size hasInertia pos neg zero choleskyOk perturbed bestRes bestBwd
//   G  structurally singular      size rank empty
//   C  over the fill cap
//   P  refused by the reader      isPattern
//
// The human-readable row goes to stdout, which the child inherits, so it lands exactly
// where it would have without the fork. If the child DIES before writing its line, the
// parent sees an empty pipe and finishes the row itself: that is the whole point, and
// GHS_indef/bloweya is the case that forced it. See the README.
void processMatrix(const std::string& path, std::size_t maxFill, std::FILE* out) {
    const std::string name = shortName(path);

    MatrixBenchmark::ReadResult file;
    try {
        file = MatrixBenchmark::readMatrixMarket(path);
    } catch (const std::exception& error) {
        std::printf("%-34s SKIP  %s\n", name.c_str(), error.what());
        std::fprintf(out, "P\t0\n");
        return;
    }
    if (!file.ok) {
        std::printf("%-34s SKIP  %s\n", name.c_str(), file.reason.c_str());
        std::fprintf(out, "P\t%d\n",
                     file.reason == "field is pattern" ? 1 : 0);
        return;
    }

    const SparseMatrix<double>& A = file.matrix;

    Vector<double> b(A.size());
    for (std::size_t k = 0; k < A.size(); ++k)
        b[k] = 1.0;

    // Once per matrix: it is the same for all three factorizations.
    const double aNorm = infNorm(A);

    // The name goes out before the work, so an abort names the matrix that caused it.
    std::printf("%-34s %7zu %9zu", name.c_str(), A.size(), A.nnz());
    std::fflush(stdout);

    // A structurally singular matrix is LABELED rather than hidden, and not solved. With b
    // chosen rather than formed from a known solution, the system is inconsistent whenever A
    // is singular, no x satisfies it, and what comes out is decided by arbitrary tiny pivots:
    // netscience read 2.19 on one machine and 1.25e+33 on another, on the same code. Printing
    // that beside honest residuals would invite reading it as a defect, which it is not. See
    // the README for the worked case and for what these rows are owed: a consistent
    // right-hand side, and the rank the factorization found.
    const std::size_t rank = MatrixBenchmark::structuralRank(A);
    if (rank < A.size()) {
        const std::size_t empty = emptyColumns(A);
        std::printf("  SINGULAR, structural rank %zu of %zu\n", rank, A.size());
        std::fprintf(out, "G\t%zu\t%zu\t%zu\n", A.size(), rank, empty);
        return;
    }

    DirectSolver<double> solver(Ordering::MmdCompacted, Factorization::Cholesky,
                                Traversal::LeftLooking);

    bool analyzed = false;
    try {
        analyzed = solver.analyze(A);
    } catch (const std::exception& error) {
        std::printf("  ANALYZE  %s\n", error.what());
        std::fprintf(out, "P\t0\n");
        return;
    }
    if (!analyzed) {
        std::printf("  ANALYZE  refused\n");
        std::fprintf(out, "P\t0\n");
        return;
    }

    const std::size_t predicted = fill(solver.symFactor());
    std::printf(" %10zu", predicted);
    std::fflush(stdout);   // so a death after this point still names the size

    // Eight bytes a value is the floor, before indices, the update matrices and
    // whatever a delayed column widens. If the prediction alone will not fit,
    // nothing downstream will.
    if (predicted > maxFill) {
        const double bytes = double(predicted) * 8.0;
        if (bytes >= 1073741824.0)
            std::printf("  TOO LARGE, %.1f GB of values at least\n", bytes / 1073741824.0);
        else
            std::printf("  TOO LARGE, %.0f MB of values at least\n", bytes / 1048576.0);
        std::fprintf(out, "C\n");
        return;
    }



    const Outcome cholesky = run(solver, Factorization::Cholesky, A, b, aNorm);
    const Outcome staticLdl = run(solver, Factorization::StaticLDLT, A, b, aNorm);
    const Outcome dynamicLdl = run(solver, Factorization::DynamicLDLT, A, b, aNorm);

    // What the factorization actually held, against what the analysis predicted, and the
    // ratio between them, which is THE PRICE OF THE PIVOTING. The two can differ ONLY under
    // dynamic pivoting: a delayed column widens its parent's front, so the gap is exactly
    // what the delays cost. A statically pivoted factor moves nothing and reads 1.0x, which
    // is every row with no Nd in its note.
    if (dynamicLdl.actualFill > 0)
        std::printf(" %10zu %6.1fx", dynamicLdl.actualFill,
                    double(dynamicLdl.actualFill) / double(predicted));
    else
        std::printf(" %10s %7s", "-", "-");

    printCell(cholesky, false);
    printCell(staticLdl, false);
    printCell(dynamicLdl, false);

    // The inertia sits at the end because it is read from a factorization that has to run
    // first, and it is taken from the DYNAMIC one, which does not perturb and so describes A
    // rather than a nearby matrix.
    if (dynamicLdl.hasInertia) {
        std::printf(" %7zu %7zu %6zu", dynamicLdl.inertia.positive,
                    dynamicLdl.inertia.negative, dynamicLdl.inertia.zero);

        // The same word the structurally singular rows carry, so the two kinds read alike at
        // a glance. They are not established alike: that one is proved from the pattern
        // before anything runs, this one is found by the factorization, and the row is
        // solved and reported rather than declined.
        if (dynamicLdl.inertia.zero > 0)
            std::printf("  SINGULAR");
    }

    std::printf("\n");

    Solved record;
    record.size = A.size();
    record.hasInertia = dynamicLdl.hasInertia;
    record.inertia = dynamicLdl.inertia;
    record.choleskyWorked = cholesky.solved;
    record.perturbed = staticLdl.perturbed;

    bool first = true;
    for (const Outcome* o : {&cholesky, &staticLdl, &dynamicLdl})
        if (o->solved) {
            if (first || o->residual < record.bestResidual)
                record.bestResidual = o->residual;
            if (first || o->backward < record.bestBackward)
                record.bestBackward = o->backward;
            first = false;
        }

    std::fprintf(out, "R\t%zu\t%d\t%zu\t%zu\t%zu\t%d\t%zu\t%.17g\t%.17g\n",
                 record.size, record.hasInertia ? 1 : 0,
                 record.inertia.positive, record.inertia.negative,
                 record.inertia.zero, record.choleskyWorked ? 1 : 0,
                 record.perturbed, record.bestResidual, record.bestBackward);
}


} // namespace

int main(int argc, char** argv) {
    // A FILL CAP, because nothing before the analysis predicts the cost. --max-nnz in ssget.py
    // bounds nnz(A) and the memory a factorization wants is nnz(L), which is a different number
    // by orders of magnitude: GHS_indef/bloweya is n = 30004 and took the whole machine down. The
    // predicted fill is known exactly once analyze() has run and costs nothing to read, so the
    // check goes there and the row is skipped with its number rather than the run dying.
    //
    // The default is generous enough to pass everything this folder has fetched so far, the
    // largest being PARSEC/benzene at 13.9 million. Raise it with --max-fill=N.
    std::size_t maxFill = 50000000;

    // An address-space limit for each child, so that where the kernel enforces it the failure
    // arrives as std::bad_alloc and reads as a skip rather than as a signal. Zero disables it.
    // 0 by default because macOS generally ignores RLIMIT_AS and a limit that is silently not
    // applied is worse than none: it would suggest a guard that is not there.
    std::size_t memoryLimit = 0;

    std::vector<const char*> files;
    for (int i = 1; i < argc; ++i) {
        const std::string arg = argv[i];
        if (arg.rfind("--max-fill=", 0) == 0) {
            const double value = optionValue(arg, 11);
            if (value < 0)
                return 1;
            maxFill = static_cast<std::size_t>(value);
        } else if (arg.rfind("--max-memory-gb=", 0) == 0) {
            const double value = optionValue(arg, 16);
            if (value < 0)
                return 1;
            memoryLimit = static_cast<std::size_t>(value * 1073741824.0);
        } else if (arg.rfind("--", 0) == 0) {
            std::fprintf(stderr, "unknown option %s\n", arg.c_str());
            return 1;
        } else
            files.push_back(argv[i]);
    }

    if (files.empty()) {
        std::fprintf(stderr,
                     "usage: %s [--max-fill=N] [--max-memory-gb=G] <file.mtx> ...\n",
                     argv[0]);
        return 1;
    }

    std::printf("ordering MmdCompacted, traversal left-looking, b all ones\n");
    std::printf("bwd = ||Ax - b|| / (||A|| ||x|| + ||b||), res = ||Ax - b|| / ||b||, "
                "infinity norms throughout\n");
    std::printf("note: Np, the static LDL perturbed N pivots it could not divide by, having no way\n      to pivot. Nd, the dynamic LDL delayed N columns to a parent rather than\n      pivot on them.\n");
    std::printf("inertia: eigenvalues of A by sign, read from D and taken from the "
                "DYNAMIC factorization,\n         which does not perturb. A nonzero "
                "`zero` is marked SINGULAR: that matrix is\n         NUMERICALLY singular, which the structural test before the run cannot see.\n\n");
    std::printf("%-34s %7s %9s %26s  %-28s %-28s %-28s %s\n",
                "", "", "", "nnz(L) ------------------",
                "Cholesky", "StaticLDLT", "DynamicLDLT", "inertia of A");
    std::printf("%-34s %7s %9s %10s %10s %7s  %8s %8s %-10s %8s %8s %-10s %8s %8s %-10s %7s %7s %6s\n",
                "matrix", "n", "nnz(A)", "predicted", "actual", "ratio",
                "bwd", "res", "", "bwd", "res", "", "bwd", "res", "",
                "pos", "neg", "zero");

    int read = 0;
    int refused = 0;         // every file we could not take
    int refusedPattern = 0;  // of those, the ones with no values at all
    int capped = 0;          // analyzed, then skipped for predicted fill
    int died = 0;            // the child process did not survive the matrix

    // The labeled rows are collected as well as printed in place. Four among sixty is easy to
    // lose while scanning a table, and they are a set worth reading as one: what we have declined
    // to answer, and how much of each matrix is the reason.
    std::vector<Labeled> singular;

    // And the solved rows, for the classification at the end.
    std::vector<Solved> solved;

    // ONE CHILD PER MATRIX. The parent forks, waits, and reads the child's one-line record off a
    // pipe. A child that exits normally reported everything it found; a child killed by a signal
    // reported nothing, and the parent finishes its half-written row and carries on.
    //
    // This exists because SIGKILL cannot be caught. The kernel destroys the process without
    // unwinding, so no handler runs, and `std::bad_alloc` is thrown only when an allocation FAILS
    // rather than when the pages are later touched under overcommit. GHS_indef/bloweya took two
    // runs down that way. Where RLIMIT_AS is enforced the child gets an exception instead and the
    // existing catch turns it into a clean skip, which is a much better message; where it is not,
    // as on macOS, the fork is what holds.
    for (const char* argument : files) {
        const std::string path = argument;
        const std::string name = shortName(path);

        std::fflush(stdout);   // or the child inherits the parent's buffer and duplicates it

        int channel[2];
        if (::pipe(channel) != 0) {
            std::fprintf(stderr, "pipe failed\n");
            return 1;
        }

        const pid_t child = ::fork();
        if (child < 0) {
            // No fork, so no guard: do the work here and accept the risk rather than stop.
            ::close(channel[0]);
            std::FILE* out = ::fdopen(channel[1], "w");
            processMatrix(path, maxFill, out);
            std::fclose(out);
            continue;
        }

        if (child == 0) {
            ::close(channel[0]);
            std::FILE* out = ::fdopen(channel[1], "w");
            if (memoryLimit > 0) {
                // Where this is honored the allocation fails and the child throws, which reads far
                // better than a signal. macOS generally ignores it, hence the fork above.
                struct rlimit limit;
                limit.rlim_cur = memoryLimit;
                limit.rlim_max = memoryLimit;
                ::setrlimit(RLIMIT_AS, &limit);
            }
            processMatrix(path, maxFill, out);
            std::fflush(out);
            std::fflush(stdout);
            std::_Exit(0);
        }

        ::close(channel[1]);

        std::string reply;
        {
            char buffer[512];
            std::FILE* in = ::fdopen(channel[0], "r");
            if (in != nullptr) {
                if (std::fgets(buffer, sizeof(buffer), in) != nullptr)
                    reply = buffer;
                std::fclose(in);
            }
        }

        int status = 0;
        ::waitpid(child, &status, 0);

        if (reply.empty()) {
            // The child said nothing, so it did not finish. Its row is half written.
            if (WIFSIGNALED(status))
                std::printf("  KILLED by signal %d, out of memory\n", WTERMSIG(status));
            else
                std::printf("  DIED without a result\n");
            ++died;
            continue;
        }

        std::vector<std::string> field;
        {
            std::size_t at = 0;
            while (at <= reply.size()) {
                const std::size_t tab = reply.find_first_of("\t\n", at);
                if (tab == std::string::npos)
                    break;
                field.push_back(reply.substr(at, tab - at));
                at = tab + 1;
            }
        }
        if (field.empty())
            continue;

        if (field[0] == "P") {
            ++refused;
            if (field.size() > 1 && field[1] == "1")
                ++refusedPattern;
        } else if (field[0] == "C") {
            ++capped;
        } else if (field[0] == "G" && field.size() >= 4) {
            Labeled row;
            row.name = name;
            row.size = std::strtoull(field[1].c_str(), nullptr, 10);
            row.rank = std::strtoull(field[2].c_str(), nullptr, 10);
            row.empty = std::strtoull(field[3].c_str(), nullptr, 10);
            singular.push_back(row);
        } else if (field[0] == "R" && field.size() >= 10) {
            Solved row;
            row.name = name;
            row.size = std::strtoull(field[1].c_str(), nullptr, 10);
            row.hasInertia = field[2] == "1";
            row.inertia.positive = std::strtoull(field[3].c_str(), nullptr, 10);
            row.inertia.negative = std::strtoull(field[4].c_str(), nullptr, 10);
            row.inertia.zero = std::strtoull(field[5].c_str(), nullptr, 10);
            row.choleskyWorked = field[6] == "1";
            row.perturbed = std::strtoull(field[7].c_str(), nullptr, 10);
            row.bestResidual = std::strtod(field[8].c_str(), nullptr);
            row.bestBackward = std::strtod(field[9].c_str(), nullptr);
            solved.push_back(row);
            ++read;
        }
    }

    if (!singular.empty()) {
        std::printf("\nlabeled structurally singular and not solved:\n");
        std::printf("%-34s %8s %8s %8s %8s\n", "matrix", "n", "s-rank", "deficit", "empty");
        for (const Labeled& row : singular)
            std::printf("%-34s %8zu %8zu %8zu %8zu\n", row.name.c_str(), row.size, row.rank,
                        row.size - row.rank, row.empty);
        std::printf("\nStructural rank is the size of a maximum matching between columns and rows,\n"
                    "so a deficit means EVERY term of the determinant vanishes and the matrix is\n"
                    "singular for almost any values. `empty` counts columns with no nonzero at\n"
                    "all, which is the common case here and a subset of the deficit.\n");
    }

    // ------------------------------------------------------------------------------------------
    // The classification, which is exact where it can be and says so where it cannot.
    // ------------------------------------------------------------------------------------------

    if (!solved.empty()) {
        std::size_t positiveDefinite = 0;
        std::size_t negativeDefinite = 0;
        std::size_t indefinite = 0;
        std::size_t numericallySingular = 0;
        std::size_t unknown = 0;

        // Order matters: a matrix with a zero eigenvalue is singular whatever the other two
        // counts say, so that test comes first and the definiteness tests only see the rest.
        for (const Solved& row : solved) {
            if (!row.hasInertia)
                ++unknown;
            else if (row.inertia.zero > 0)
                ++numericallySingular;
            else if (row.inertia.positive == row.size)
                ++positiveDefinite;
            else if (row.inertia.negative == row.size)
                ++negativeDefinite;
            else
                ++indefinite;
        }

        std::printf("\nof the %d files read:\n",
                    read + refused + (int)singular.size() + capped + died);
        std::printf("  pattern, no values    %4d   nothing to solve, and nothing lost: an\n"
                    "                               ordering study would want exactly these\n",
                    refusedPattern);
        if (refused - refusedPattern > 0)
            std::printf("  refused for other     %4d   see the reason on the row\n",
                        refused - refusedPattern);
        std::printf("  with values           %4d\n",
                    read + (int)singular.size() + capped + died);
        std::printf("    structurally singular %3zu   declined before factoring, listed above\n",
                    singular.size());
        if (capped > 0)
            std::printf("    over the fill cap     %3d   analyzed, then skipped; raise it with\n"
                        "                               --max-fill=N\n", capped);
        if (died > 0)
            std::printf("    killed               %4d   the child did not survive the matrix, so\n"
                        "                               the row above says how far it got\n", died);
        std::printf("    solved                %3zu\n", solved.size());

        std::printf("\nof the %zu solved, by the inertia of A:\n", solved.size());
        std::printf("  positive definite     %4zu   every eigenvalue strictly positive\n",
                    positiveDefinite);
        std::printf("  negative definite     %4zu   every eigenvalue strictly negative; Cholesky\n"
                    "                               must refuse these, and -A would factor\n",
                    negativeDefinite);
        std::printf("  indefinite            %4zu   both signs present, none zero\n", indefinite);
        std::printf("  numerically singular  %4zu   at least one exactly zero, whatever the other\n"
                    "                               signs. Counted ONLY here, not in the three\n"
                    "                               rows above, the zero test coming first\n",
                    numericallySingular);
        if (unknown > 0)
            std::printf("  no inertia            %4zu\n", unknown);

        // A free consistency check with no threshold in it: Cholesky succeeds exactly on a
        // positive definite matrix, so the two accounts must agree on every row.
        std::size_t disagreements = 0;
        for (const Solved& row : solved)
            if (row.hasInertia
                && row.choleskyWorked != (row.inertia.positive == row.size && row.inertia.zero == 0))
                ++disagreements;

        if (disagreements == 0)
            std::printf("\nCholesky agrees with the inertia on all %zu rows: it succeeded on every\n"
                        "positive definite matrix and on no other, which is what it is for.\n",
                        solved.size());
        else
            std::printf("\nCHOLESKY DISAGREES WITH THE INERTIA ON %zu ROWS. One of the two is wrong.\n",
                        disagreements);

        // WHERE WE FAILED TO PRODUCE A USABLE ANSWER, ranked rather than thresholded. The
        // question this answers is not "is the matrix singular" but "did the best of our three
        // factorizations bring the residual down", and a ranking needs no cutoff to be invented
        // for it: the worst rows are the worst rows whatever the rest of the table looks like.
        //
        // An earlier version of this block selected on `res / bwd` and returned 21 rows of 35,
        // which is not a pointer. That quantity is roughly ||A|| ||x|| / ||b||, and it grows with
        // the SCALING of A as much as with its conditioning, so a matrix whose entries are large
        // ranks high for reasons that have nothing to do with difficulty. The tell was `bwd`
        // figures below machine epsilon, which no solve achieves: their denominator was simply
        // large.
        std::vector<const Solved*> worst;
        for (const Solved& row : solved)
            worst.push_back(&row);
        std::sort(worst.begin(), worst.end(), [](const Solved* a, const Solved* b) {
            return a->bestResidual > b->bestResidual;
        });

        const std::size_t show = std::min<std::size_t>(worst.size(), 8);
        std::printf("\nthe %zu largest residuals, taking the best of the three factorizations on\n"
                    "EACH MEASURE SEPARATELY, which may be two different factorizations:\n", show);
        std::printf("%-34s %10s %10s %10s %7s\n",
                    "matrix", "best res", "best bwd", "perturbed", "zero");
        for (std::size_t i = 0; i < show; ++i) {
            const Solved* row = worst[i];
            std::printf("%-34s %10.1e %10.1e %10zu %7zu\n", row->name.c_str(),
                        row->bestResidual, row->bestBackward, row->perturbed,
                        row->hasInertia ? row->inertia.zero : 0);
        }
        std::printf("\nSo `best res` says whether ANY factorization got the residual down and\n"
                    "`best bwd` whether ANY was backward stable. A large residual beside a `bwd`\n"
                    "at machine precision is the matrix and not us: no x reproduces b, either\n"
                    "because none exists or because ||x|| is so large that forming Ax loses b\n"
                    "entirely, and `zero` says which of those is proved. A large residual beside a\n"
                    "large `bwd` and no perturbations would be ours, and there is no such row\n"
                    "today.\n");
    }

    std::printf("\n%d solved, %zu labeled singular, %d over the fill cap, %d killed, %d skipped\n",
                read, singular.size(), capped, died, refused);
    return 0;
}
