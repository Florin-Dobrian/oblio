# Notes for Oblio, from the polyglot repo

Written 2026-08-02, in a session on the `data-structures` playground, while implementing dense and
sparse matvec and dense LU factor and solve in C++, Python, Scala and Rust. Two questions about
Oblio surfaced from the cross-language contrast. Neither is a change to make now; both are worth
recording before the reasoning is lost.

Nothing here transfers a convention from the playground. That repo holds four languages in parallel
and can only adopt a rule all four can carry; Oblio is one language chosen deliberately and can tune
a rule to what C++ makes easy. Different constraints, different answers. What travels is the
reasoning, not the conclusion.

## 1. The index-type rule, and whether it should narrow to storage

**What.** CODING_RULES currently says a loop counter takes the type of what it counts: a loop over
entities counts indices, so the counter is `std::int32_t` and the bound is cast once in the
condition. The playground arrived at a narrower rule instead, `std::int32_t` for indices that are
*stored* in an array and `std::size_t` for loop counters, with the cast at the point a stored index
becomes a subscript.

**Why the narrower rule is attractive.** The 32-bit cap earns its keep on storage, where an index
array is as long as the matrix has nonzeros and its element width is a first-order memory cost. A
loop counter is ephemeral and costs nothing either way, so making it signed buys no memory and adds
a cast per loop condition.

**Why Oblio's rule is probably still right, and this is likely to close with no change.** The
narrower rule was forced in the playground by Rust, which subscripts only with `usize`: an `i32`
counter there would give up `for i in 0..n` and turn every loop into a hand-driven `while`. That
pressure does not exist in C++, and Oblio's rule buys something the narrow one does not:

- **Symmetry between ascending and descending entity loops.** With `std::int32_t` the two directions
  have the same shape, `for (i = 0; i < n; ++i)` and `for (i = n - 1; i >= 0; --i)`, and the loop
  exits at `-1`, which the type represents. With `std::size_t` the descending form needs
  `for (std::size_t i = n; i-- > 0;)`, which is a real trick: it exploits the return value of
  post-decrement to get a termination test the type cannot express, and a reader has to verify that
  it enters at `n - 1` and exits after `0`.
- The dense LU solve in the playground has four such loops, two of them descending, and the
  asymmetry is visible on the page.

Oblio is full of descending sweeps, so the uniformity is worth more here than the casts cost.

**What to do.** Record the reasoning as a DESIGN_DECISIONS entry rather than change anything. The
current rule is a considered position and the entry should say why, including what the alternative
buys and why it does not apply. A rule that has survived a stated challenge is stronger than one
that has never been questioned.

## 2. Return by value where the result is one-shot

**What.** Oblio's engines are consistently out-parameter APIs: `compute(A, Permutation&)`,
`MultiplyEngine` filling a `y` rather than returning one. Some of that is a C++98 artifact, since
0.9 was developed between 1998 and 2005, starting the year C++98 was standardized and ending six
years before move semantics arrived, throughout which returning a container risked a full copy.

**Why C++11 changed this.** Before move constructors, `return y;` on a named local had two
outcomes: NRVO applied and the object was built directly in the caller's storage, or NRVO did not
apply and the *copy* constructor ran, allocating a second buffer and copying every element. NRVO
was optional and could not handle every shape, so returning containers by value was a habit to
avoid. C++11 made a returned local an rvalue, so the fallback became a move rather than a copy:
stealing a pointer, never copying a buffer. Measured in the playground, `return y;` produces exactly
one constructor call at `-O0` and `-O3` alike, and with `-fno-elide-constructors` the fallback is a
move. The bad case stopped being expensive, which is what made return-by-value safe to write without
thinking.

The companion idiom is a by-value parameter plus `std::move` into the member:

```
Vector(std::size_t size, std::vector<double> val)
    : mSize(size), mVal(std::move(val)) {}
```

This costs one copy for an lvalue argument and zero for an rvalue, where a `const&` parameter costs
one copy always. The `std::move` is never a pessimization and should be there by default. The one
thing to watch is a constructor body that reads the parameter after the initializer list has moved
from it, which sees an emptied object; the fix is to hoist the read into the initializer list ahead
of the move, which is where the size guard already sits.

**Why this is not a blanket change.** Two distinctions matter, and separating them is the actual
work:

- **Some out-params are correct in any era.** A numeric kernel filling a caller-supplied buffer is
  how we factor once and solve many times without reallocating per call. BLAS and LAPACK are
  out-parameter APIs today for exactly that reason. These should stay as they are.
- **Some are access design rather than history.** An engine fills an object it is a `friend` of, and
  that is the deliberate write-grant rule in DESIGN_DECISIONS, not a workaround for missing moves.

So the candidates are the one-shot results: a function whose output is constructed, returned once
and not written into again. Those would read better by value and now cost nothing to return.

**Why not now.** Changing a signature from out-param to return-by-value changes what a caller does,
which puts it on the rewrite track by CLAUDE.md's own definition of the port invariant. Nothing is
`verified` buffer-for-buffer against 0.9 yet, so an API change would land ahead of the comparison it
is supposed to survive.

**What to do.** A TODO entry naming the specific candidates, one line each, split into the three
groups above: one-shot results worth returning, hot-loop buffers that stay out-params, and
friend-filled objects that stay as they are. Writing the inventory is cheap and can happen now; the
change waits for the port to reach `verified`.
