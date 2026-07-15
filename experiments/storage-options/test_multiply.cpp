// test_multiply.cpp - one multiply function, two storages.
//
// Builds one sparse matrix, stores it twice (CSC and vector-of-vectors) from the same
// entries, and runs the *same* multiply source over both, each reading its matrix through
// that storage's own per-column lookups (rowIdxPtr / valPtr / colLen).
//
// Correctness is the first result and the point of the experiment: the two runs must agree
// bit for bit, because they are the same additions in the same order over the same values,
// reached through pointers that differ only in where they point.
//
// Timing is the second result, and it says something the design needs: how much the layout
// costs, and how much the pointer indirection costs, which turn out to be very different
// numbers.

#include "MultiplyEngine.h"
#include "SparseMatrixStatic.h"
#include "SparseMatrixDynamic.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <random>
#include <vector>

using namespace StorageOptions;

namespace {

// The matrix, as per-column lists. Both storages are built from this, so their content is
// identical by construction and the experiment measures layout, not two different matrices.
struct Columns {
    std::size_t size;
    std::vector<std::vector<std::int32_t>> rowIdx;
    std::vector<std::vector<double>>       val;
};

Columns build(std::size_t size, std::size_t bandwidth, std::mt19937& rng) {
    Columns c;
    c.size = size;
    c.rowIdx.resize(size);
    c.val.resize(size);
    std::uniform_real_distribution<double> uniform(-1.0, 1.0);

    for (std::int32_t j = 0; j < static_cast<std::int32_t>(size); ++j) {
        const std::int32_t hi = std::min(static_cast<std::int32_t>(size),
                                         j + static_cast<std::int32_t>(bandwidth));
        for (std::int32_t i = j; i < hi; ++i) {
            c.rowIdx[j].push_back(i);
            c.val[j].push_back(uniform(rng));
        }
    }
    return c;
}

SparseMatrixStatic toCsc(const Columns& c) {
    std::vector<std::size_t>  colPtr(c.size + 1, 0);
    std::vector<std::int32_t> rowIdx;
    std::vector<double>       val;
    for (std::int32_t j = 0; j < static_cast<std::int32_t>(c.size); ++j) {
        rowIdx.insert(rowIdx.end(), c.rowIdx[j].begin(), c.rowIdx[j].end());
        val.insert(val.end(), c.val[j].begin(), c.val[j].end());
        colPtr[j + 1] = rowIdx.size();
    }
    return SparseMatrixStatic(c.size, std::move(colPtr), std::move(rowIdx), std::move(val));
}

// Columns allocated in order. The allocator lays them out nearly contiguously, so this is the
// most favourable case a vector of vectors can have: a flat layout with an extra hop.
SparseMatrixDynamic toVv(const Columns& c) {
    return SparseMatrixDynamic(c.size, c.rowIdx, c.val);
}

// The same content, but the columns allocated in random order with other allocations
// interleaved, so they land wherever the allocator has room. This is what a vector of vectors
// looks like once it has been *used*: in dynamic factorization, fronts grow and reallocate at
// different times, and end up scattered. The version above flatters VV; this one does not.
SparseMatrixDynamic toVvScattered(const Columns& c, std::mt19937& rng) {
    std::vector<std::int32_t> order(c.size);
    for (std::int32_t j = 0; j < static_cast<std::int32_t>(c.size); ++j) order[j] = j;
    std::shuffle(order.begin(), order.end(), rng);

    std::vector<std::vector<std::int32_t>> rowIdx(c.size);
    std::vector<std::vector<double>>       val(c.size);
    std::vector<std::vector<double>>       spacer;   // pushes the next column elsewhere
    spacer.reserve(c.size);

    for (std::int32_t k = 0; k < static_cast<std::int32_t>(c.size); ++k) {
        const std::int32_t j = order[k];
        rowIdx[j] = c.rowIdx[j];
        val[j]    = c.val[j];
        spacer.emplace_back(8, 0.0);
    }
    return SparseMatrixDynamic(c.size, std::move(rowIdx), std::move(val));
}

double maxDiff(const std::vector<double>& a, const std::vector<double>& b) {
    double d = 0.0;
    for (std::size_t i = 0; i < a.size(); ++i)
        d = std::max(d, std::abs(a[i] - b[i]));
    return d;
}

template<class Fn>
double timeIt(Fn&& fn, int repeats) {
    using Clock = std::chrono::steady_clock;
    fn();   // warm up
    const auto t0 = Clock::now();
    for (int r = 0; r < repeats; ++r) fn();
    const auto t1 = Clock::now();
    return std::chrono::duration<double, std::milli>(t1 - t0).count() / repeats;
}

} // namespace

// Exercise the asymmetry: setValues on both, setColumn on the dynamic one only.
//
// This is the part of the experiment that is about the *interface*, not the speed. An object
// offers what its storage makes cheap; the static matrix therefore has no setColumn, and that
// absence is not an omission to be fixed.
void testMutation() {
    std::cout << "\nMutation: what each storage is willing to promise\n";

    // A 3x3, so it can be checked by hand.  A = [[2,1,0],[1,3,1],[0,1,4]]
    std::vector<std::size_t>  colPtr = {0, 2, 5, 7};
    std::vector<std::int32_t> rowIdx = {0,1, 0,1,2, 1,2};
    std::vector<double>       val    = {2,1, 1,3,1, 1,4};

    SparseMatrixStatic  s(3, colPtr, rowIdx, val);
    SparseMatrixDynamic d(3, {{0,1}, {0,1,2}, {1,2}}, {{2,1}, {1,3,1}, {1,4}});

    MultiplyEngine eng;
    std::vector<double> x = {1, 1, 1}, y(3);

    // multiply() ACCUMULATES into y (the BLAS convention: y += A x), so y must be zeroed before
    // each call. The benchmark above does the same. Forgetting it is exactly the kind of quiet
    // wrong answer that a test with no expected value would never catch.
    auto run = [&](auto& A) {
        std::fill(y.begin(), y.end(), 0.0);
        eng.multiply(A, x.data(), y.data());
    };

    // setValues: BOTH support it, and it costs nothing in either. Same structure, new numbers, is
    // the mutation the solver actually does most (a Newton step, a time step, refactorize). The
    // calls are now identical on both classes, which is the point of dropping to column granularity.
    const bool sOk = s.setValues(0, {20,10}) && s.setValues(1, {10,30,10}) && s.setValues(2, {10,40});
    const bool dOk = d.setValues(0, {20,10}) && d.setValues(1, {10,30,10}) && d.setValues(2, {10,40});

    run(s);
    const bool sVal = sOk && y[0]==30 && y[1]==50 && y[2]==50;

    run(d);
    const bool dVal = dOk && y[0]==30 && y[1]==50 && y[2]==50;

    std::cout << "  setValues, static   " << (sVal ? "ok" : "FAILED")
              << "   (cheap: nothing moves)\n";
    std::cout << "  setValues, dynamic  " << (dVal ? "ok" : "FAILED")
              << "   (cheap: nothing moves)\n";

    // setColumn: the DYNAMIC one only. Give column 1 a different pattern entirely.
    const bool cOk = d.setColumn(1, {0, 2}, {7.0, 9.0});
    run(d);
    // now A = [[20,7,0],[10,0,10],[0,9,40]],  A*(1,1,1) = (27, 20, 49)
    const bool dCol = cOk && y[0]==27 && y[1]==20 && y[2]==49;
    std::cout << "  setColumn, dynamic  " << (dCol ? "ok" : "FAILED")
              << "   (cheap: the column owns its buffer, neighbours untouched)\n";

    // And rejection: out of range, unsorted, mismatched lengths.
    const bool rej = !d.setColumn(9, {0}, {1.0})          // column out of range
                  && !d.setColumn(0, {1, 0}, {1.0, 2.0})  // rows not sorted
                  && !d.setColumn(0, {0, 1}, {1.0});      // lengths disagree
    std::cout << "  setColumn rejects   " << (rej ? "ok" : "FAILED")
              << "   (bad column, unsorted rows, mismatched lengths)\n";

    std::cout << "  setColumn, static   absent by design"
              << "   (it would shift every later column: O(nnz), not O(column))\n";
}

// The invalidation rule, demonstrated rather than merely asserted.
//
// **Structural mutation invalidates every pointer previously extracted. Value mutation does not.**
// That is the whole rule, it holds in both storages, and it is the one the solver's dynamic factor
// will live by: delayed pivoting grows a front, which reallocates its buffer, which dangles every
// pointer into it. A use-after-free there is silent.
//
// Nothing in C++ enforces this. So the least we can do is show it: extract the pointers, mutate,
// and observe which pointers still point where they did.
void testInvalidation() {
    std::cout << "\nPointer validity: which mutations move the buffers\n";

    SparseMatrixDynamic d(3, {{0,1}, {0,1,2}, {1,2}}, {{2,1}, {1,3,1}, {1,4}});

    const double* before = d.valPtr(1);       // where column 1's values lived

    // setValues: the buffers stay put. The old pointer still points at column 1, and now sees the
    // new numbers. Safe to keep, though few callers should want to.
    d.setValues(0, {20,10}); d.setValues(1, {10,30,10}); d.setValues(2, {10,40});
    const bool valueKept = (d.valPtr(1) == before);
    std::cout << "  setValues   pointer " << (valueKept ? "UNCHANGED" : "moved    ")
              << "   (buffer reused; contents overwritten in place)\n";

    // setColumn: the buffer is replaced. The old pointer is dangling, and we must not read it.
    // We can only observe that it is no longer where the column lives.
    d.setColumn(1, {0, 2}, {7.0, 9.0});
    const bool structMoved = (d.valPtr(1) != before);
    std::cout << "  setColumn   pointer " << (structMoved ? "MOVED    " : "unchanged")
              << "   (buffer replaced; anything held from before now dangles)\n";

    std::cout << "  rule                " << ((valueKept && structMoved) ? "holds" : "BROKEN")
              << "        structural mutation invalidates, value mutation does not\n";
}

int main() {
    const std::size_t size      = 200000;
    const std::size_t bandwidth = 16;
    const int         repeats   = 20;

    std::mt19937 rng(20260712);
    const Columns c = build(size, bandwidth, rng);

    const SparseMatrixStatic csc  = toCsc(c);
    const SparseMatrixDynamic  vv   = toVv(c);
    const SparseMatrixDynamic  vvS  = toVvScattered(c, rng);

    MultiplyEngine eng;

    std::vector<double> x(size);
    for (std::int32_t j = 0; j < static_cast<std::int32_t>(size); ++j)
        x[j] = 1.0 + 0.001 * static_cast<double>(j % 97);

    std::vector<double> yBase(size, 0.0), yCsc(size, 0.0), yVv(size, 0.0), yVvS(size, 0.0);

    std::cout << "storage-options: one multiply, two storages (static and dynamic)\n"
              << "  size = " << size << ", nnz = " << csc.nnz()
              << ", repeats = " << repeats << "\n\n";

    const double tBase = timeIt([&]{ std::fill(yBase.begin(), yBase.end(), 0.0);
                                     eng.multiplyStatic(csc, x.data(), yBase.data()); }, repeats);
    const double tCsc  = timeIt([&]{ std::fill(yCsc.begin(), yCsc.end(), 0.0);
                                     eng.multiply(csc, x.data(), yCsc.data()); }, repeats);
    const double tVv   = timeIt([&]{ std::fill(yVv.begin(), yVv.end(), 0.0);
                                     eng.multiply(vv, x.data(), yVv.data()); }, repeats);
    const double tVvS  = timeIt([&]{ std::fill(yVvS.begin(), yVvS.end(), 0.0);
                                     eng.multiply(vvS, x.data(), yVvS.data()); }, repeats);

    const double dCsc = maxDiff(yBase, yCsc);
    const double dVv  = maxDiff(yBase, yVv);
    const double dVvS = maxDiff(yBase, yVvS);
    const bool   ok   = (dCsc == 0.0 && dVv == 0.0 && dVvS == 0.0);

    std::cout << std::fixed << std::setprecision(3);
    std::cout << "  hand-written flat (baseline) " << std::setw(8) << tBase << " ms\n";
    std::cout << "  multiply(), static (direct)  " << std::setw(8) << tCsc << " ms   "
              << std::setprecision(2) << tCsc / tBase << "x\n" << std::setprecision(3);
    std::cout << "  multiply(), dynamic, packed  " << std::setw(8) << tVv << " ms   "
              << std::setprecision(2) << tVv / tBase << "x\n" << std::setprecision(3);
    std::cout << "  multiply(), dynamic, scatter " << std::setw(8) << tVvS << " ms   "
              << std::setprecision(2) << tVvS / tBase << "x\n\n";

    std::cout << "  all results identical: " << (ok ? "yes" : "NO") << "\n";
    if (!ok)
        std::cout << "    max diff  csc " << dCsc << "  vv " << dVv
                  << "  vv-scattered " << dVvS << "\n";

    testMutation();
    testInvalidation();

    std::cout << "\nWhat this shows:\n"
              << "  The three timed multiply() rows are one multiply SOURCE, monomorphized per\n"
              << "  storage. Each instantiation reads its matrix through that storage's own\n"
              << "  lookups (rowIdxPtr / valPtr / colLen) and names no member, no buffer, no\n"
              << "  layout, so the source is written once and the compiler specializes it. This\n"
              << "  is direct access, the consumer templated on the storage, which is exactly the\n"
              << "  shape the numeric engine uses on the factor through blockPtr. There is no\n"
              << "  extractor and no pointer-array view, so nothing is owned and nothing dangles.\n"
              << "\n"
              << "  Reaching a column through the lookup is free. multiply() over the static\n"
              << "  matrix matches the hand-written flat baseline, which calls no lookup at all\n"
              << "  and walks colPtr raw, so the abstraction costs nothing on the layout we care\n"
              << "  most about.\n"
              << "\n"
              << "  What costs is the layout. The two dynamic rows are the same class holding\n"
              << "  the same content, differing only in the order their inner vectors were\n"
              << "  allocated. A flat buffer guarantees consecutive columns are adjacent; a\n"
              << "  vector of vectors only ever borrows that from the allocator.\n"
              << "\n"
              << "  Two caveats on the scattered number. It is constructed, not observed: we\n"
              << "  shuffle the allocation order and interleave spacers, to remove the\n"
              << "  allocator's help entirely. And this kernel is the harshest setting a cache\n"
              << "  miss can have, roughly two flops per element loaded, so a miss shows at\n"
              << "  full price. A factorization front handed to BLAS-3 does O(n^3) work on\n"
              << "  O(n^2) data, and the same miss amortizes. Read the packed row as the\n"
              << "  dynamic layout's structural cost and the scattered row as an upper bound,\n"
              << "  not a forecast.\n";

    return ok ? 0 : 1;
}
