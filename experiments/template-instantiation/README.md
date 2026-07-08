# Template-Instantiation Example

A minimal, self-contained study of three ways to instantiate a `Val`-templated
class for `double` and `std::complex<double>`. Reference / teaching only — **not**
part of the main Oblio build.

## The question

Oblio's value-dependent classes are templated on the scalar type but only ever used
with a small, fixed set (`double`, `complex<double>`). How should those templates be
compiled? This example builds the same tiny stack — `Matrix`, `Vector`, and a
`MultiplyEngine` doing dense `y = A·x` — three ways, and shows they produce identical
results while differing only in *how* the code is compiled.

## The three variants

Files are named `<Class><Variant>` (e.g. `MatrixPlainExplicit.h`):

- **Implicit** (`*Implicit.h`) — full template body in the header; every translation
  unit instantiates it implicitly. This is what Oblio 0.9 effectively did.
- **PlainExplicit** (`*PlainExplicit.{h,cpp}`) — declaration-only header; bodies plus
  explicit instantiation (`template class Matrix<double>;`) in the `.cpp`.
- **GuardedExplicit** (`*GuardedExplicit.{h,cpp}`) — same as PlainExplicit, plus
  `extern template` lines in the header.

## Build & run

```
make test     # builds and runs all three; they must print identical results
make clean
```

One shared source, `test_multiply.cpp`, compiles against each variant via a macro
(`-DOBLIO_TI_IMPLICIT` / `-DOBLIO_TI_PLAIN_EXPLICIT` / `-DOBLIO_TI_GUARDED_EXPLICIT`).
Executables carry the `_cpp` language suffix and are gitignored.

## What it demonstrates

- All three compute the same thing — the variant is purely a compilation strategy.
- Implicit re-instantiates the template in *every* translation unit (the build cost
  0.9 paid). The explicit variants compile each specialization once, in the `.cpp`.
- PlainExplicit and GuardedExplicit behave identically, at runtime and at link time.
  `extern template` suppresses implicit instantiation of *visible* header bodies — but
  a declaration-only header has none, so here it suppresses nothing: it is
  documentation, not mechanism.
- Bonus: build either explicit variant *without* its `.cpp` files and it fails at link
  with undefined references — proof the bodies live only in the compiled objects.

## Decision it informed

The real tree uses **guarded explicit**: declaration-only headers, bodies + explicit
instantiation in `.cpp`, `extern template` in headers as an intent annotation. Full
rationale — the two-axis framing (implicit vs explicit; plain vs guarded) and the
C++98/C++11 history — is in `../../docs/DESIGN_DECISIONS.md` (the explicit-instantiation
entry).
