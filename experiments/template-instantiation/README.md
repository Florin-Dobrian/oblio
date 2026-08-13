# Template-Instantiation Example

A minimal, self-contained study of three ways to instantiate a `Val`-templated
class for `double` and `std::complex<double>`. Reference / teaching only, **not**
part of the main Oblio build.

## The question

Oblio's value-dependent classes are templated on the scalar type but only ever used
with a small, fixed set (`double`, `complex<double>`). How should those templates be
compiled? This example builds the same tiny stack, `Matrix`, `Vector`, and a
`MultiplyEngine` doing dense `y = A*x`, three ways, and shows they produce identical
results while differing only in *how* the code is compiled.

## The three variants

Files are named `<Class><Variant>` (e.g. `MatrixPlainExplicit.h`):

- **Implicit** (`*Implicit.h`), full template body in the header; every translation
  unit instantiates it implicitly. This is what Oblio 0.9 effectively did.
- **PlainExplicit** (`*PlainExplicit.{h,cpp}`), declaration-only header apart from one
  defaulted default constructor; bodies plus explicit instantiation
  (`template class Matrix<double>;`) in the `.cpp`.
- **GuardedExplicit** (`*GuardedExplicit.{h,cpp}`), same as PlainExplicit, plus
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

- All three compute the same thing, the variant is purely a compilation strategy.
- Implicit re-instantiates the template in *every* translation unit (the build cost
  0.9 paid). The explicit variants compile each specialization once, in the `.cpp`.
- PlainExplicit and GuardedExplicit behave identically, at runtime and at link time, and
  the guard is documentation rather than mechanism here. Each header carries exactly one
  visible body, the defaulted default constructor, and it was worth checking whether that
  gives `extern template` something to act on. **It does not, and the check is the
  interesting part.** Linking a program that default-constructs both classes *without* the
  `.cpp` files leaves the same three symbols undefined under each variant,
  `Matrix<double>::rows()`, `Matrix<double>::cols()` and `Vector<double>::size()`;
  `Matrix<double>::Matrix()` is in neither list. A
  defaulted default constructor over scalars and `std::vector` members is trivial, so no
  out-of-line function is emitted for it at all, and there is no symbol for the guard to
  suppress or for the linker to resolve. The prediction that the guard would gain mechanism
  on that member was wrong, and measuring it took ten minutes.

  **Checked on both toolchains**, since "no symbol is emitted" is an implementation matter
  rather than a guarantee, and this tree has been caught before by a negative result on one
  compiler being read as a result about another: Apple clang on arm64 and GCC 13 on x86-64
  name the same three symbols. Count the symbols and not the lines if this is ever re-run,
  because the two linkers report differently, GNU `ld` one line per use site and Apple's
  grouping the uses under each symbol.
- So the shape to take from this is the *rule*, not a mechanism: a body in the header is
  safe **provided** the class stays explicitly instantiated with the guard present. Whether
  the guard has work to do depends on the body, and for a trivial one it has none. The real
  tree is this arrangement, which is why it is worth modelling: `include/oblio/Vector.h`
  carries `Vector() = default;` in an otherwise declaration-only header beside its
  `extern template` lines, and CLAUDE.md's definitions-in-cpp invariant names exactly this
  as its deliberate exception. An accidental inline without the guard is the bug; a chosen
  one under it is fine.
- Bonus: build either explicit variant *without* its `.cpp` files and it fails at link
  with undefined references, proof the bodies live only in the compiled objects. Unchanged
  by the above, and the same symbols go missing under both variants.

## Decision it informed

The real tree uses **guarded explicit**: declaration-only headers, bodies + explicit
instantiation in `.cpp`, `extern template` in headers as an intent annotation. Full
rationale, the two-axis framing (implicit vs explicit; plain vs guarded) and the
C++98/C++11 history, is in `../../docs/DESIGN_DECISIONS.md` (the explicit-instantiation
entry).
