// matrix_performance.cpp -- where a solve's time goes, per ordering and per traversal, on real
// matrices.
//
// `../pipeline` asks this question on grid Laplacians and says in its own header that a break-even
// computed there may be a property of grids rather than of the orderings. This is the same
// question with the input replaced, which is the whole of what is new: the phase split, the
// timing protocol and the loop nesting are all taken from that folder deliberately, so a row here
// can be read beside a row there.
//
//   make performance                        every matrix in performance_candidates.txt
//   ./matrix_performance_cpp ../../data/HB/*.mtx      any subset, the shell choosing it
//   ./matrix_performance_cpp --repeats=1 ...          a quick pass
//
// CHOLESKY ONLY, and that is not a limitation here but the point. Cholesky never pivots, so:
//
//   - the factor's structure is exactly what the analysis predicted, which makes nnz(L) a property
//     of the ordering alone and therefore comparable across the four;
//   - no column is delayed, so no row can surprise the memory budget the way a dynamically
//     pivoted one can, and the fill cap below is exact rather than a lower bound;
//   - and the numeric phase does the same arithmetic under every traversal, so a difference
//     between left-looking, right-looking and multifrontal is a difference in how the work is
//     scheduled rather than in how much of it there is.
//
// Everything here is therefore about the ordering and the traversal, with the factorization held
// still. The accuracy of these factorizations is the subject of ACCURACY.md and is not measured
// again.
//
// THE SET IS THE POSITIVE DEFINITE ONE. Cholesky refuses anything else, so `ssget.py --posdef`
// selects the matrices this driver can use, and `performance_candidates.txt` is the list. A matrix
// that refuses is reported and stepped over, as everywhere in this folder.
//
// FOUR ORDERINGS, TWO NAMED. MMD and AMD are the vendored codes and need `../../private`; MMD3 and
// AMD3 are Oblio's. The vendored pair is the reference: it says whether our implementations cost
// what they should and fill what they should. Where private/ is absent those two rows report a
// refusal and the other two still run.

// ASK FOR A PERFORMANCE CORE, on the one platform where cores differ. Apple Silicon runs a
// command-line process at QOS_CLASS_DEFAULT, which prefers a performance core but permits the
// scheduler to park the thread on an efficiency one, and that placement is STICKY over long
// stretches rather than jittering per iteration, so a minimum over repeats does not filter it.
// The same call in ../ordering turned a scattered null into a clean result on 2026-08-10.
#ifdef __APPLE__
#include <pthread.h>
#endif

#include "MatrixMarket.h"

#include "oblio/DirectSolver.h"
#include "oblio/ElmForest.h"
#include "oblio/ElmForestEngine.h"
#include "oblio/MultiplyEngine.h"
#include "oblio/OrderEngine.h"
#include "oblio/Permutation.h"
#include "oblio/SymFactor.h"
#include "oblio/SymFactorEngine.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>

using namespace Oblio;

namespace {

// The group and file name, which is the collection's own identifier for a matrix.
std::string shortName(const std::string& path) {
    const std::size_t slash = path.find_last_of('/');
    if (slash == std::string::npos)
        return path;
    const std::size_t before = path.find_last_of('/', slash - 1);
    if (before == std::string::npos)
        return path;
    std::string name = path.substr(before + 1);
    if (name.size() > 4 && name.compare(name.size() - 4, 4, ".mtx") == 0)
        name.erase(name.size() - 4);
    return name;
}

// The induced infinity norm, for the backward error's denominator. Computed from the columns,
// which is the same thing here: A is symmetric and stored fully.
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

// nnz(L) from the symbolic factor, a supernode's own triangle plus its update rows. The same
// computation ../ordering and ../pipeline make, duplicated rather than shared on the same grounds:
// a benchmark should stand alone.
std::size_t fill(const SymFactor& sf) {
    std::size_t nnz = 0;
    for (std::int32_t kk = 0; kk < static_cast<std::int32_t>(sf.snodeSize()); ++kk) {
        const std::size_t f = sf.frontSize(kk);
        const std::size_t u = sf.updateSize(kk);
        nnz += f * (f + 1) / 2 + f * u;
    }
    return nnz;
}

// Best of N after a warm-up, the work a callable so every phase is timed the same way. Taken from
// ../pipeline, including the warm-up, which matters more here than there: a matrix read from disk
// has just been touched by the reader and its pages are cold in a way a generated grid's are not.
template <class Work>
double bestOf(Work work, int repeats) {
    work();
    double best = 1e30;
    for (int trial = 0; trial < repeats; ++trial) {
        const auto t0 = std::chrono::steady_clock::now();
        work();
        const auto t1 = std::chrono::steady_clock::now();
        best = std::min(best, std::chrono::duration<double, std::milli>(t1 - t0).count());
    }
    return best;
}

struct Timing {
    bool        ok = false;
    std::size_t nnzL = 0;
    double      order = 0;      // OrderEngine::compute alone
    double      analyze = 0;    // the whole facade call, left-looking
    double      analyzeMf = 0;  // the same, multifrontal: it sorts children and relabels
    double      factLl = 0;
    double      factRl = 0;
    double      factMf = 0;
    double      solve = 0;

    // THE ANSWER, beside the cost. A timing report that does not say whether the answers were any
    // good is half a report, and it costs nothing here: the solve runs for the timing anyway, so
    // the residual is one extra multiply. Both measures are computed as in ../matrices and mean
    // the same there: `bwd` is the verdict, since it does not inherit the conditioning, and `res`
    // is reported beside it. Every matrix in this set is positive definite and factored by
    // Cholesky, so nothing perturbs and nothing is delayed; these columns should be at machine
    // precision on every row and the run says so at the end if they are not.
    double backward = 0;
    double residual = 0;

    // PER TRAVERSAL, because they do not always agree. Cholesky's refusal is a numerical event,
    // not a structural one: the three traversals accumulate the same arithmetic in different
    // orders, so a matrix sitting at the definiteness boundary can factor under one and refuse
    // under another. HB/plat1919 is the case, right-looking refusing under three orderings of the
    // four, and it is the matrix the accuracy run put at a residual of 1.0e+03 with a backward
    // error of 7.3e-17. Without these flags a refusal printed as a time of 0.00, which reads as
    // the fastest cell in the table.
    bool okLl = false;
    bool okRl = false;
    bool okMf = false;
};

// The fill an ordering produces, and whether it can be produced at all. Cheap: analysis only, no
// arithmetic. Run for every ordering before any timing, so a matrix too large for the cap is
// skipped whole rather than half measured.
std::size_t fillFor(const SparseMatrix<double>& A, Ordering method, bool& ok) {
    const OrderEngine oe(method);
    Permutation       P;
    ok = false;
    if (!oe.compute(A, P))
        return 0;

    const ElmForestEngine fe;
    ElmForest             ef;
    if (!fe.compute(A, P, ef))
        return 0;

    const SymFactorEngine se;
    SymFactor             sf;
    if (!se.compute(A, P, ef, sf))
        return 0;

    ok = true;
    return fill(sf);
}

// One traversal's factor time, with the analysis already done and outside the timed region.
// Reported separately from analyze because analyze depends on the traversal a little, multifrontal
// sorting children and relabeling supernodes, and factor depends on it a great deal.
double factorTime(const SparseMatrix<double>& A, Ordering method, Traversal traversal, int repeats,
                  bool& ok) {
    DirectSolver<double> solver(method, Factorization::Cholesky, traversal);
    ok = solver.analyze(A);
    if (!ok)
        return 0;
    if (!solver.factor(A)) {   // not positive definite, which Cholesky is entitled to refuse
        ok = false;
        return 0;
    }
    return bestOf([&] { solver.factor(A); }, repeats);
}

Timing measure(const SparseMatrix<double>& A, Ordering method, int repeats, double aNorm) {
    Timing t;

    bool ok = false;
    t.nnzL = fillFor(A, method, ok);
    if (!ok)
        return t;

    const OrderEngine oe(method);
    t.order = bestOf([&] { Permutation P; oe.compute(A, P); }, repeats);

    t.analyze = bestOf([&] {
        DirectSolver<double> s(method, Factorization::Cholesky, Traversal::LeftLooking);
        s.analyze(A);
    }, repeats);

    t.analyzeMf = bestOf([&] {
        DirectSolver<double> s(method, Factorization::Cholesky, Traversal::Multifrontal);
        s.analyze(A);
    }, repeats);

    t.factLl = factorTime(A, method, Traversal::LeftLooking, repeats, t.okLl);
    if (!t.okLl)
        return t;   // Cholesky refused: not positive definite, and the row says so
    t.factRl = factorTime(A, method, Traversal::RightLooking, repeats, t.okRl);
    t.factMf = factorTime(A, method, Traversal::Multifrontal, repeats, t.okMf);

    DirectSolver<double> solver(method, Factorization::Cholesky, Traversal::LeftLooking);
    if (solver.analyze(A) && solver.factor(A)) {
        Vector<double> b(A.size()), x(A.size()), r(A.size());
        for (std::size_t i = 0; i < A.size(); ++i)
            b[i] = 1.0;
        t.solve = bestOf([&] { solver.solve(b, x); }, repeats);

        const MultiplyEngine multiply;
        multiply.residual(A, x, b, r);

        const double rNorm = infNorm(r);
        const double xNorm = infNorm(x);
        const double bNorm = infNorm(b);

        t.backward = rNorm / (aNorm * xNorm + bNorm);
        t.residual = (bNorm > 0.0) ? rNorm / bNorm : 0.0;
    }

    t.ok = true;
    return t;
}

// The four orderings, in the order the table prints them: each of ours beside the vendored code it
// reproduces, so the pair can be read across rather than hunted for.
struct Method {
    Ordering    method;
    const char* name;
};

const Method kMethods[] = {
    {Ordering::MMD,  "MMD"},
    {Ordering::MMD3, "MMD3"},
    {Ordering::AMD,  "AMD"},
    {Ordering::AMD3, "AMD3"},
};
const int kNumMethods = 4;

// Everything one matrix produced, kept for the summary.
struct Record {
    std::string name;
    Timing      timing[kNumMethods];
};

// The median, which is the right average for a SHARE. Shares are bounded in [0, 1] and their
// distribution across a heterogeneous set is skewed by the small matrices, where the analysis is
// everything and the factorization is nothing. A median says what a typical matrix looks like; a
// mean would be dragged by the tail.
double median(std::vector<double> values) {
    if (values.empty())
        return 0;
    std::sort(values.begin(), values.end());
    const std::size_t half = values.size() / 2;
    if (values.size() % 2 == 1)
        return values[half];
    return 0.5 * (values[half - 1] + values[half]);
}

// The geometric mean, which is the right average for a ratio: it treats 2x and 0.5x as equal and
// opposite where the arithmetic mean does not, and it is what makes "on average 4 percent more
// fill" mean the same thing read in either direction.
double geoMean(const std::vector<double>& values) {
    if (values.empty())
        return 0;
    double sum = 0;
    for (double v : values)
        sum += std::log(v);
    return std::exp(sum / double(values.size()));
}

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

} // namespace

int main(int argc, char** argv) {
#ifdef __APPLE__
    // Before anything is timed. See the note beside the include.
    pthread_set_qos_class_self_np(QOS_CLASS_USER_INTERACTIVE, 0);
#endif

    // The cap is on nnz(L), which for Cholesky is exactly what the factor will hold. It is checked
    // for EVERY ordering before any timing, so a matrix one ordering cannot fit is skipped whole
    // rather than reported with three rows and a gap.
    std::size_t maxFill = 50000000;
    int         repeats = 3;

    std::vector<const char*> files;
    for (int i = 1; i < argc; ++i) {
        const std::string arg = argv[i];
        if (arg.rfind("--max-fill=", 0) == 0) {
            const double value = optionValue(arg, 11);
            if (value < 0)
                return 1;
            maxFill = static_cast<std::size_t>(value);
        } else if (arg.rfind("--repeats=", 0) == 0) {
            const double value = optionValue(arg, 10);
            if (value < 0)
                return 1;
            repeats = static_cast<int>(value);
        } else if (arg.rfind("--", 0) == 0) {
            std::fprintf(stderr, "unknown option %s\n", arg.c_str());
            return 1;
        } else
            files.push_back(argv[i]);
    }

    if (files.empty()) {
        std::fprintf(stderr, "usage: %s [--max-fill=N] [--repeats=N] <file.mtx> ...\n", argv[0]);
        return 1;
    }

    std::printf("Cholesky, real, best of %d after a warm-up, milliseconds\n", repeats);
    std::printf("order is OrderEngine::compute alone; analyze is the whole facade call, so the\n"
                "forest and symbolic cost is the difference and no phase is double counted\n\n");

    std::vector<Record> records;
    int refused = 0;
    int capped = 0;
    int notDefinite = 0;
    double worstBackward = 0;   // over every ordering of every matrix
    double worstResidual = 0;   // read and analyzed, but Cholesky refused under every ordering

    for (const char* argument : files) {
        const std::string path = argument;
        const std::string name = shortName(path);

        MatrixBenchmark::ReadResult file;
        try {
            file = MatrixBenchmark::readMatrixMarket(path);
        } catch (const std::exception& error) {
            std::printf("%-36s SKIP  %s\n", name.c_str(), error.what());
            ++refused;
            continue;
        }
        if (!file.ok) {
            std::printf("%-36s SKIP  %s\n", name.c_str(), file.reason.c_str());
            ++refused;
            continue;
        }

        const SparseMatrix<double>& A = file.matrix;

        // The fill of every ordering first, so the cap decides before anything is timed.
        std::size_t largest = 0;
        for (int m = 0; m < kNumMethods; ++m) {
            bool ok = false;
            largest = std::max(largest, fillFor(A, kMethods[m].method, ok));
        }
        if (largest > maxFill) {
            std::printf("%-36s %8zu %10zu  TOO LARGE, nnz(L) %zu\n",
                        name.c_str(), A.size(), A.nnz(), largest);
            ++capped;
            continue;
        }

        const double aNorm = infNorm(A);

        std::printf("%-36s %8zu %10zu\n", name.c_str(), A.size(), A.nnz());
        std::printf("  %-6s %11s %8s %8s %8s %8s %8s %8s %8s %9s %9s\n",
                    "order", "nnz(L)", "order", "anlzLL", "anlzMF", "factLL", "factRL", "factMF",
                    "solve", "bwd", "res");
        std::fflush(stdout);

        Record record;
        record.name = name;

        for (int m = 0; m < kNumMethods; ++m) {
            const Timing t = measure(A, kMethods[m].method, repeats, aNorm);
            record.timing[m] = t;

            if (!t.ok) {
                std::printf("  %-6s %11s\n", kMethods[m].name,
                            t.nnzL == 0 ? "refused" : "not SPD");
                continue;
            }
            std::printf("  %-6s %11zu %8.2f %8.2f %8.2f %8.2f", kMethods[m].name, t.nnzL,
                        t.order, t.analyze, t.analyzeMf, t.factLl);
            if (t.okRl) std::printf(" %8.2f", t.factRl); else std::printf(" %8s", "refused");
            if (t.okMf) std::printf(" %8.2f", t.factMf); else std::printf(" %8s", "refused");
            std::printf(" %8.2f %9.1e %9.1e\n", t.solve, t.backward, t.residual);

            worstBackward = std::max(worstBackward, t.backward);
            worstResidual = std::max(worstResidual, t.residual);
        }
        std::fflush(stdout);

        // A matrix no ordering could factor is not a measurement. Cholesky refusing means the
        // matrix is not positive definite, which the index said it was, so the count is worth
        // keeping separate rather than folding into the total.
        bool any = false;
        for (int m = 0; m < kNumMethods; ++m)
            if (record.timing[m].ok)
                any = true;

        if (any)
            records.push_back(record);
        else
            ++notDefinite;
    }

    // ------------------------------------------------------------------------------------------
    // The summary, which is the part a report quotes.
    // ------------------------------------------------------------------------------------------

    if (!records.empty()) {
        // Only the matrices where every ordering succeeded, so the averages compare like with
        // like. A matrix where the vendored pair refused would otherwise flatter our two.
        std::vector<const Record*> complete;
        for (const Record& r : records) {
            bool all = true;
            for (int m = 0; m < kNumMethods; ++m)
                if (!r.timing[m].ok)
                    all = false;
            if (all)
                complete.push_back(&r);
        }

        std::printf("\nall four orderings succeeded on %zu of the %zu matrices measured\n",
                    complete.size(), records.size());

        if (!complete.empty()) {
            std::printf("\ngeometric mean relative to the best ordering on each matrix:\n");
            std::printf("  %-6s %10s %10s %10s %10s\n",
                        "order", "nnz(L)", "order", "analyze", "factLL");

            for (int m = 0; m < kNumMethods; ++m) {
                std::vector<double> fillRatio, orderRatio, analyzeRatio, factRatio;

                for (const Record* r : complete) {
                    double bestFill = 1e300, bestOrder = 1e300;
                    double bestAnalyze = 1e300, bestFact = 1e300;
                    for (int k = 0; k < kNumMethods; ++k) {
                        bestFill = std::min(bestFill, double(r->timing[k].nnzL));
                        bestOrder = std::min(bestOrder, r->timing[k].order);
                        bestAnalyze = std::min(bestAnalyze, r->timing[k].analyze);
                        bestFact = std::min(bestFact, r->timing[k].factLl);
                    }
                    if (bestFill > 0)    fillRatio.push_back(double(r->timing[m].nnzL) / bestFill);
                    if (bestOrder > 0)   orderRatio.push_back(r->timing[m].order / bestOrder);
                    if (bestAnalyze > 0) analyzeRatio.push_back(r->timing[m].analyze / bestAnalyze);
                    if (bestFact > 0)    factRatio.push_back(r->timing[m].factLl / bestFact);
                }

                std::printf("  %-6s %10.3f %10.3f %10.3f %10.3f\n", kMethods[m].name,
                            geoMean(fillRatio), geoMean(orderRatio), geoMean(analyzeRatio),
                            geoMean(factRatio));
            }

            // WHERE A ONE-SHOT SOLVE'S TIME GOES, at a fixed ordering and the multifrontal
            // traversal. Four shares that sum to one:
            //
            //   ordering            OrderEngine::compute alone
            //   rest of analysis    the elimination forest and the symbolic factorization, which
            //                       is the whole analysis minus the ordering
            //   factorization       the numeric phase
            //   solve               the triangular solves
            //
            // Reported at the median and again over the largest quartile by nnz(L), on the
            // expectation that the split would move with size. IT LARGELY DOES NOT: on 2026-08-11
            // the two columns came out at 21.9 against 21.7 percent for the ordering and 46.7
            // against 47.9 for the factorization. Individual large matrices do lean further
            // toward the factorization, up to three quarters, but the analysis grows with the
            // factor and the medians stay put. The second column is kept because that was worth
            // establishing rather than assuming.
            //
            // Each median is taken over its own column, so the four DO NOT sum to 100 and are not
            // meant to: a median is not additive. They say what a typical matrix's share of each
            // phase looks like, not how one matrix divides.
            //
            // The analysis is also the part that DISAPPEARS on reuse. A caller factoring the same
            // pattern repeatedly pays it once, so these shares are the one-shot case and the
            // upper bound on what the analysis can cost.
            {
                std::vector<const Record*> bySize(complete);
                std::sort(bySize.begin(), bySize.end(), [](const Record* a, const Record* b) {
                    return a->timing[0].nnzL > b->timing[0].nnzL;
                });
                const std::size_t quartile = std::max<std::size_t>(1, bySize.size() / 4);

                std::vector<double> ordAll, restAll, factAll, solveAll;
                std::vector<double> ordBig, restBig, factBig, solveBig;

                for (std::size_t i = 0; i < bySize.size(); ++i) {
                    const Timing& t = bySize[i]->timing[0];
                    if (!t.okMf)
                        continue;
                    const double total = t.analyzeMf + t.factMf + t.solve;
                    if (total <= 0)
                        continue;

                    const double ord = t.order / total;
                    const double rest = (t.analyzeMf - t.order) / total;
                    const double fact = t.factMf / total;
                    const double sol = t.solve / total;

                    ordAll.push_back(ord);   restAll.push_back(rest);
                    factAll.push_back(fact); solveAll.push_back(sol);
                    if (i < quartile) {
                        ordBig.push_back(ord);   restBig.push_back(rest);
                        factBig.push_back(fact); solveBig.push_back(sol);
                    }
                }

                std::printf("\nwhere a one-shot solve's time goes, at %s and multifrontal,\n"
                            "median share of the total:\n", kMethods[0].name);
                std::printf("  %-20s %10s %14s\n", "phase", "all", "largest 25%");
                std::printf("  %-20s %9.1f%% %13.1f%%\n", "ordering",
                            100.0 * median(ordAll), 100.0 * median(ordBig));
                std::printf("  %-20s %9.1f%% %13.1f%%\n", "rest of analysis",
                            100.0 * median(restAll), 100.0 * median(restBig));
                std::printf("  %-20s %9.1f%% %13.1f%%\n", "factorization",
                            100.0 * median(factAll), 100.0 * median(factBig));
                std::printf("  %-20s %9.1f%% %13.1f%%\n", "solve",
                            100.0 * median(solveAll), 100.0 * median(solveBig));
                std::printf("  (each median is taken over its own column, so the four need not\n"
                            "   sum to 100. The analysis is paid once and reused on a repeated\n"
                            "   factorization, so these are the one-shot shares.)\n");
            }

            // The traversals, at a fixed ordering, since that is the axis they vary.
            std::vector<double> ll, rl, mf;
            std::size_t traversalRefusals = 0;
            for (const Record* r : complete) {
                const Timing& t = r->timing[0];
                if (!t.okLl || !t.okRl || !t.okMf) {
                    ++traversalRefusals;   // a matrix all three cannot factor is not a comparison
                    continue;
                }
                const double best = std::min({t.factLl, t.factRl, t.factMf});
                if (best > 0) {
                    ll.push_back(t.factLl / best);
                    rl.push_back(t.factRl / best);
                    mf.push_back(t.factMf / best);
                }
            }
            std::printf("\ntraversals at %s, geometric mean of the factor time relative to\n"
                        "the best traversal on each matrix:\n", kMethods[0].name);
            std::printf("  left-looking  %6.3f\n  right-looking %6.3f\n  multifrontal  %6.3f\n",
                        geoMean(ll), geoMean(rl), geoMean(mf));
            if (traversalRefusals > 0)
                std::printf("  (%zu %s left out, a traversal having refused where another\n"
                            "   succeeded; see the refused cells in the table)\n",
                            traversalRefusals,
                            traversalRefusals == 1 ? "matrix" : "matrices");
        }
    }

    // The answers, in one line. Cholesky on a positive definite matrix perturbs nothing and delays
    // nothing, so these should be at machine precision on every row; a report of timings that does
    // not say whether the answers were good is half a report.
    if (!records.empty())
        std::printf("\nworst backward error over every ordering of every matrix %.1e. Cholesky on\n"
                    "a positive definite matrix neither perturbs nor delays, so that is the\n"
                    "arithmetic alone. The worst relative residual is %.1e, which is conditioning\n"
                    "rather than error: see ACCURACY.md on why the two are reported together.\n",
                    worstBackward, worstResidual);

    std::printf("\n%zu measured, %d not positive definite, %d over the fill cap, %d skipped\n",
                records.size(), notDefinite, capped, refused);
    return 0;
}
