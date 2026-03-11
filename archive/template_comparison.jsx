import { useState } from "react";

const STYLES = `
  @import url('https://fonts.googleapis.com/css2?family=JetBrains+Mono:ital,wght@0,300;0,400;0,500;0,600;1,300&family=Syne:wght@400;600;700&display=swap');

  * { box-sizing: border-box; margin: 0; padding: 0; }

  body { background: #0e0e10; }

  .root {
    font-family: 'Syne', sans-serif;
    background: #0e0e10;
    color: #e8e6e0;
    min-height: 100vh;
    padding: 2rem;
  }

  .header {
    text-align: center;
    margin-bottom: 2.5rem;
  }

  .header h1 {
    font-size: 1.5rem;
    font-weight: 700;
    letter-spacing: 0.04em;
    color: #f0ede6;
    margin-bottom: 0.4rem;
  }

  .header p {
    font-size: 0.82rem;
    color: #6e6b63;
    font-family: 'JetBrains Mono', monospace;
    font-weight: 300;
    letter-spacing: 0.02em;
  }

  .tabs {
    display: flex;
    justify-content: center;
    gap: 0.25rem;
    margin-bottom: 2rem;
    background: #16161a;
    border: 1px solid #252529;
    border-radius: 8px;
    padding: 0.3rem;
    width: fit-content;
    margin-left: auto;
    margin-right: auto;
  }

  .tab {
    font-family: 'JetBrains Mono', monospace;
    font-size: 0.75rem;
    font-weight: 500;
    padding: 0.45rem 1.1rem;
    border-radius: 5px;
    border: none;
    cursor: pointer;
    transition: all 0.15s ease;
    letter-spacing: 0.03em;
    background: transparent;
    color: #5a5750;
  }

  .tab:hover { color: #9e9a90; }

  .tab.active {
    background: #252529;
    color: #c8c4bc;
  }

  .columns {
    display: grid;
    grid-template-columns: 1fr 1fr;
    gap: 1.25rem;
    max-width: 1280px;
    margin: 0 auto;
  }

  .panel {
    background: #131316;
    border: 1px solid #222226;
    border-radius: 10px;
    overflow: hidden;
    display: flex;
    flex-direction: column;
  }

  .panel-header {
    display: flex;
    align-items: center;
    justify-content: space-between;
    padding: 0.75rem 1.1rem;
    border-bottom: 1px solid #1e1e22;
    background: #111114;
  }

  .panel-title {
    font-size: 0.78rem;
    font-weight: 600;
    letter-spacing: 0.06em;
    text-transform: uppercase;
  }

  .panel-title.tpl { color: #7a9ecf; }
  .panel-title.exp { color: #9ecf8a; }

  .panel-filename {
    font-family: 'JetBrains Mono', monospace;
    font-size: 0.68rem;
    color: #3e3e44;
    font-weight: 300;
  }

  .panel-body {
    padding: 1.1rem 1.25rem;
    overflow-x: auto;
    flex: 1;
  }

  pre {
    font-family: 'JetBrains Mono', monospace;
    font-size: 0.72rem;
    line-height: 1.75;
    font-weight: 300;
  }

  .note-box {
    margin: 0.6rem 1.25rem 1rem;
    background: #1a1a1e;
    border-left: 2px solid #333338;
    border-radius: 0 5px 5px 0;
    padding: 0.65rem 0.9rem;
  }

  .note-box p {
    font-size: 0.73rem;
    color: #5e5b55;
    font-family: 'JetBrains Mono', monospace;
    font-weight: 300;
    line-height: 1.6;
  }

  .note-box p strong {
    color: #7a7570;
    font-weight: 500;
  }

  /* Syntax colours */
  .kw  { color: #c792ea; }   /* keyword */
  .tp  { color: #7a9ecf; }   /* type / template param */
  .fn  { color: #ddc790; }   /* function name */
  .cm  { color: #3d3d44; font-style: italic; }   /* comment */
  .st  { color: #9ecf8a; }   /* string / literal */
  .pp  { color: #c86e6e; }   /* preprocessor */
  .op  { color: #6e9ecf; }   /* operator / punctuation accent */
  .nm  { color: #e8c88a; }   /* number */
  .pl  { color: #e8e6e0; }   /* plain */

  .divider {
    height: 1px;
    background: #1e1e22;
    margin: 0.5rem 1.25rem;
  }

  .section-label {
    font-family: 'JetBrains Mono', monospace;
    font-size: 0.62rem;
    font-weight: 500;
    letter-spacing: 0.1em;
    text-transform: uppercase;
    color: #35353c;
    padding: 0.5rem 1.25rem 0.25rem;
  }

  .summary {
    max-width: 1280px;
    margin: 1.5rem auto 0;
    display: grid;
    grid-template-columns: 1fr 1fr;
    gap: 1.25rem;
  }

  .summary-panel {
    background: #131316;
    border: 1px solid #222226;
    border-radius: 10px;
    padding: 1rem 1.25rem;
  }

  .summary-panel h3 {
    font-size: 0.72rem;
    font-weight: 600;
    letter-spacing: 0.08em;
    text-transform: uppercase;
    margin-bottom: 0.75rem;
  }

  .summary-panel.tpl h3 { color: #7a9ecf; }
  .summary-panel.exp h3 { color: #9ecf8a; }

  .summary-panel ul {
    list-style: none;
    display: flex;
    flex-direction: column;
    gap: 0.4rem;
  }

  .summary-panel li {
    font-family: 'JetBrains Mono', monospace;
    font-size: 0.69rem;
    color: #5a5750;
    font-weight: 300;
    line-height: 1.5;
    padding-left: 1rem;
    position: relative;
  }

  .summary-panel li::before {
    content: '·';
    position: absolute;
    left: 0;
  }

  .summary-panel.tpl li::before { color: #4a6e9f; }
  .summary-panel.exp li::before { color: #6e9f5a; }

  .summary-panel li em {
    font-style: normal;
    color: #8a8680;
  }
`;

// ─── Code fragments ──────────────────────────────────────────────────────────

const SECTIONS = [
  {
    id: "matrix",
    label: "Matrix",
    tpl: {
      file: "Matrix.h",
      note: "Full implementation in the header. Instantiated implicitly for every translation unit that includes it, for every Val type used.",
      code: [
        { t: "cm", v: "// Matrix.h — full implementation in header" },
        { t: "pp", v: "#pragma once" },
        { t: "pl", v: "" },
        { t: "kw", v: "template", s: " <" },
        { t: "kw", v: "class", s: " " },
        { t: "tp", v: "Val", s: ">" },
        { t: "kw", v: "class", s: " " },
        { t: "tp", v: "Matrix", s: " {" },
        { t: "kw", v: "public", s: ":" },
        { t: "pl", v: "  Matrix()" },
        { t: "pl", v: "    : mRows(0), mCols(0) {}" },
        { t: "pl", v: "" },
        { t: "pl", v: "  Matrix(" },
        { t: "tp", v: "  size_t", s: " rows, " },
        { t: "tp", v: "size_t", s: " cols," },
        { t: "kw", v: "    const", s: " std::vector<" },
        { t: "tp", v: "Val", s: ">& vals)" },
        { t: "pl", v: "    : mRows(rows), mCols(cols)," },
        { t: "pl", v: "      mVals(vals) {}" },
        { t: "pl", v: "" },
        { t: "tp", v: "  Val", s: " operator()(" },
        { t: "tp", v: "size_t", s: " i, " },
        { t: "tp", v: "size_t", s: " j) " },
        { t: "kw", v: "const", s: " {" },
        { t: "kw", v: "    return", s: " mVals[i * mCols + j];" },
        { t: "pl", v: "  }" },
        { t: "pl", v: "" },
        { t: "tp", v: "  size_t", s: " rows() " },
        { t: "kw", v: "const", s: " { " },
        { t: "kw", v: "return", s: " mRows; }" },
        { t: "tp", v: "  size_t", s: " cols() " },
        { t: "kw", v: "const", s: " { " },
        { t: "kw", v: "return", s: " mCols; }" },
        { t: "pl", v: "" },
        { t: "kw", v: "private", s: ":" },
        { t: "tp", v: "  size_t", s: " mRows, mCols;" },
        { t: "pl", v: "  std::vector<" },
        { t: "tp", v: "Val", s: "> mVals;" },
        { t: "pl", v: "};" },
        { t: "pl", v: "" },
        { t: "cm", v: "// No .cc file — everything is here." },
        { t: "cm", v: "// Recompiled for every Val in every .cc" },
        { t: "cm", v: "// that includes this header." },
      ],
    },
    exp: {
      file: "Matrix.h + Matrix.cc",
      note: "Header holds only the declaration. Implementation moves to a .cc file. Two explicit instantiations cover all needed types; adding float is one extra line.",
      hcode: [
        { t: "cm", v: "// Matrix.h — declaration only" },
        { t: "pp", v: "#pragma once" },
        { t: "pl", v: "" },
        { t: "kw", v: "template", s: " <" },
        { t: "kw", v: "class", s: " " },
        { t: "tp", v: "Val", s: ">" },
        { t: "kw", v: "class", s: " " },
        { t: "tp", v: "Matrix", s: " {" },
        { t: "kw", v: "public", s: ":" },
        { t: "pl", v: "  Matrix();" },
        { t: "pl", v: "  Matrix(" },
        { t: "tp", v: "    size_t", s: " rows, " },
        { t: "tp", v: "size_t", s: " cols," },
        { t: "kw", v: "    const", s: " std::vector<" },
        { t: "tp", v: "Val", s: ">& vals);" },
        { t: "tp", v: "  Val", s: " operator()(" },
        { t: "tp", v: "size_t", s: " i, " },
        { t: "tp", v: "size_t", s: " j) " },
        { t: "kw", v: "const", s: ";" },
        { t: "tp", v: "  size_t", s: " rows() " },
        { t: "kw", v: "const", s: ";" },
        { t: "tp", v: "  size_t", s: " cols() " },
        { t: "kw", v: "const", s: ";" },
        { t: "kw", v: "private", s: ":" },
        { t: "tp", v: "  size_t", s: " mRows, mCols;" },
        { t: "pl", v: "  std::vector<" },
        { t: "tp", v: "Val", s: "> mVals;" },
        { t: "pl", v: "};" },
        { t: "pl", v: "" },
        { t: "cm", v: "// Suppress implicit instantiation" },
        { t: "cm", v: "// in all other translation units." },
        { t: "kw", v: "extern template class", s: " " },
        { t: "tp", v: "Matrix", s: "<" },
        { t: "kw", v: "double", s: ">;" },
        { t: "kw", v: "extern template class", s: " " },
        { t: "tp", v: "Matrix", s: "<std::complex<" },
        { t: "kw", v: "double", s: ">>;" },
      ],
      ccode: [
        { t: "cm", v: "// Matrix.cc — implementation + instantiations" },
        { t: "pp", v: '#include "Matrix.h"' },
        { t: "pl", v: "" },
        { t: "kw", v: "template", s: " <" },
        { t: "kw", v: "class", s: " " },
        { t: "tp", v: "Val", s: ">" },
        { t: "tp", v: "Matrix", s: "<" },
        { t: "tp", v: "Val", s: ">::Matrix()" },
        { t: "pl", v: "  : mRows(0), mCols(0) {}" },
        { t: "pl", v: "" },
        { t: "kw", v: "template", s: " <" },
        { t: "kw", v: "class", s: " " },
        { t: "tp", v: "Val", s: ">" },
        { t: "tp", v: "Matrix", s: "<" },
        { t: "tp", v: "Val", s: ">::Matrix(" },
        { t: "tp", v: "    size_t", s: " rows, " },
        { t: "tp", v: "size_t", s: " cols," },
        { t: "kw", v: "    const", s: " std::vector<" },
        { t: "tp", v: "Val", s: ">& vals)" },
        { t: "pl", v: "  : mRows(rows), mCols(cols)," },
        { t: "pl", v: "    mVals(vals) {}" },
        { t: "pl", v: "" },
        { t: "kw", v: "template", s: " <" },
        { t: "kw", v: "class", s: " " },
        { t: "tp", v: "Val", s: ">" },
        { t: "tp", v: "Val", s: " " },
        { t: "tp", v: "Matrix", s: "<" },
        { t: "tp", v: "Val", s: ">::operator()(" },
        { t: "tp", v: "    size_t", s: " i, " },
        { t: "tp", v: "size_t", s: " j) " },
        { t: "kw", v: "const", s: " {" },
        { t: "kw", v: "  return", s: " mVals[i * mCols + j];" },
        { t: "pl", v: "}" },
        { t: "pl", v: "" },
        { t: "cm", v: "// Explicit instantiations — compiled once." },
        { t: "kw", v: "template class", s: " " },
        { t: "tp", v: "Matrix", s: "<" },
        { t: "kw", v: "double", s: ">;" },
        { t: "kw", v: "template class", s: " " },
        { t: "tp", v: "Matrix", s: "<std::complex<" },
        { t: "kw", v: "double", s: ">>;" },
        { t: "cm", v: "// Adding float later: one line here." },
      ],
    },
  },
  {
    id: "vector",
    label: "Vector",
    tpl: {
      file: "Vector.h",
      note: "Same pattern as Matrix — fully in the header, implicitly instantiated everywhere it is included.",
      code: [
        { t: "cm", v: "// Vector.h — full implementation in header" },
        { t: "pp", v: "#pragma once" },
        { t: "pl", v: "" },
        { t: "kw", v: "template", s: " <" },
        { t: "kw", v: "class", s: " " },
        { t: "tp", v: "Val", s: ">" },
        { t: "kw", v: "class", s: " " },
        { t: "tp", v: "Vector", s: " {" },
        { t: "kw", v: "public", s: ":" },
        { t: "pl", v: "  Vector();" },
        { t: "pl", v: "  " },
        { t: "kw", v: "explicit", s: " Vector(" },
        { t: "tp", v: "size_t", s: " size)" },
        { t: "pl", v: "    : mSize(size), mVals(size) {}" },
        { t: "pl", v: "" },
        { t: "tp", v: "  Val", s: "& operator[](" },
        { t: "tp", v: "size_t", s: " i) {" },
        { t: "kw", v: "    return", s: " mVals[i];" },
        { t: "pl", v: "  }" },
        { t: "kw", v: "  const", s: " " },
        { t: "tp", v: "Val", s: "& operator[](" },
        { t: "tp", v: "size_t", s: " i) " },
        { t: "kw", v: "const", s: " {" },
        { t: "kw", v: "    return", s: " mVals[i];" },
        { t: "pl", v: "  }" },
        { t: "pl", v: "" },
        { t: "tp", v: "  size_t", s: " size() " },
        { t: "kw", v: "const", s: " { " },
        { t: "kw", v: "return", s: " mSize; }" },
        { t: "pl", v: "" },
        { t: "kw", v: "private", s: ":" },
        { t: "tp", v: "  size_t", s: " mSize;" },
        { t: "pl", v: "  std::vector<" },
        { t: "tp", v: "Val", s: "> mVals;" },
        { t: "pl", v: "};" },
      ],
    },
    exp: {
      file: "Vector.h + Vector.cc",
      note: "Declaration in header with extern template to suppress implicit instantiation. Implementation and explicit instantiations in the .cc file.",
      hcode: [
        { t: "cm", v: "// Vector.h — declaration only" },
        { t: "pp", v: "#pragma once" },
        { t: "pl", v: "" },
        { t: "kw", v: "template", s: " <" },
        { t: "kw", v: "class", s: " " },
        { t: "tp", v: "Val", s: ">" },
        { t: "kw", v: "class", s: " " },
        { t: "tp", v: "Vector", s: " {" },
        { t: "kw", v: "public", s: ":" },
        { t: "pl", v: "  Vector();" },
        { t: "kw", v: "  explicit", s: " Vector(" },
        { t: "tp", v: "size_t", s: " size);" },
        { t: "tp", v: "  Val", s: "& operator[](" },
        { t: "tp", v: "size_t", s: " i);" },
        { t: "kw", v: "  const", s: " " },
        { t: "tp", v: "Val", s: "& operator[](" },
        { t: "tp", v: "size_t", s: " i) " },
        { t: "kw", v: "const", s: ";" },
        { t: "tp", v: "  size_t", s: " size() " },
        { t: "kw", v: "const", s: ";" },
        { t: "kw", v: "private", s: ":" },
        { t: "tp", v: "  size_t", s: " mSize;" },
        { t: "pl", v: "  std::vector<" },
        { t: "tp", v: "Val", s: "> mVals;" },
        { t: "pl", v: "};" },
        { t: "pl", v: "" },
        { t: "kw", v: "extern template class", s: " " },
        { t: "tp", v: "Vector", s: "<" },
        { t: "kw", v: "double", s: ">;" },
        { t: "kw", v: "extern template class", s: " " },
        { t: "tp", v: "Vector", s: "<std::complex<" },
        { t: "kw", v: "double", s: ">>;" },
      ],
      ccode: [
        { t: "cm", v: "// Vector.cc" },
        { t: "pp", v: '#include "Vector.h"' },
        { t: "pl", v: "" },
        { t: "kw", v: "template", s: " <" },
        { t: "kw", v: "class", s: " " },
        { t: "tp", v: "Val", s: ">" },
        { t: "tp", v: "Vector", s: "<" },
        { t: "tp", v: "Val", s: ">::Vector()" },
        { t: "pl", v: "  : mSize(0) {}" },
        { t: "pl", v: "" },
        { t: "kw", v: "template", s: " <" },
        { t: "kw", v: "class", s: " " },
        { t: "tp", v: "Val", s: ">" },
        { t: "tp", v: "Vector", s: "<" },
        { t: "tp", v: "Val", s: ">::Vector(" },
        { t: "tp", v: "    size_t", s: " size)" },
        { t: "pl", v: "  : mSize(size), mVals(size) {}" },
        { t: "pl", v: "" },
        { t: "cm", v: "// ... other method definitions ..." },
        { t: "pl", v: "" },
        { t: "kw", v: "template class", s: " " },
        { t: "tp", v: "Vector", s: "<" },
        { t: "kw", v: "double", s: ">;" },
        { t: "kw", v: "template class", s: " " },
        { t: "tp", v: "Vector", s: "<std::complex<" },
        { t: "kw", v: "double", s: ">>;" },
      ],
    },
  },
  {
    id: "matvec",
    label: "MatVec",
    tpl: {
      file: "MultiplyEngine.h",
      note: "The multiply function is fully templated in the header. It is recompiled for each Val type in every translation unit that includes MultiplyEngine.h.",
      code: [
        { t: "cm", v: "// MultiplyEngine.h — full impl in header" },
        { t: "pp", v: "#pragma once" },
        { t: "pp", v: '#include "Matrix.h"' },
        { t: "pp", v: '#include "Vector.h"' },
        { t: "pl", v: "" },
        { t: "kw", v: "class", s: " " },
        { t: "tp", v: "MultiplyEngine", s: " {" },
        { t: "kw", v: "public", s: ":" },
        { t: "kw", v: "  template", s: " <" },
        { t: "kw", v: "class", s: " " },
        { t: "tp", v: "Val", s: ">" },
        { t: "pl", v: "  " },
        { t: "tp", v: "Vector", s: "<" },
        { t: "tp", v: "Val", s: "> Multiply(" },
        { t: "kw", v: "    const", s: " " },
        { t: "tp", v: "Matrix", s: "<" },
        { t: "tp", v: "Val", s: ">& A," },
        { t: "kw", v: "    const", s: " " },
        { t: "tp", v: "Vector", s: "<" },
        { t: "tp", v: "Val", s: ">& x) " },
        { t: "kw", v: "const", s: " {" },
        { t: "tp", v: "    Vector", s: "<" },
        { t: "tp", v: "Val", s: "> y(A.rows());" },
        { t: "kw", v: "    for", s: " (" },
        { t: "tp", v: "size_t", s: " i = " },
        { t: "nm", v: "0", s: "; i < A.rows(); ++i) {" },
        { t: "tp", v: "      Val", s: " sum{" },
        { t: "nm", v: "0", s: "};" },
        { t: "kw", v: "      for", s: " (" },
        { t: "tp", v: "size_t", s: " j = " },
        { t: "nm", v: "0", s: "; j < A.cols(); ++j)" },
        { t: "pl", v: "        sum += A(i,j) * x[j];" },
        { t: "pl", v: "      y[i] = sum;" },
        { t: "pl", v: "    }" },
        { t: "kw", v: "    return", s: " y;" },
        { t: "pl", v: "  }" },
        { t: "pl", v: "};" },
        { t: "pl", v: "" },
        { t: "cm", v: "// Works for any Val automatically." },
        { t: "cm", v: "// But recompiled in every .cc that" },
        { t: "cm", v: "// includes this header." },
      ],
    },
    exp: {
      file: "MultiplyEngine.h + MultiplyEngine.cc",
      note: "The function signature is declared in the header with extern template suppression. The implementation and explicit instantiations live in the .cc, compiled exactly once.",
      hcode: [
        { t: "cm", v: "// MultiplyEngine.h — declaration only" },
        { t: "pp", v: "#pragma once" },
        { t: "pp", v: '#include "Matrix.h"' },
        { t: "pp", v: '#include "Vector.h"' },
        { t: "pl", v: "" },
        { t: "kw", v: "class", s: " " },
        { t: "tp", v: "MultiplyEngine", s: " {" },
        { t: "kw", v: "public", s: ":" },
        { t: "kw", v: "  template", s: " <" },
        { t: "kw", v: "class", s: " " },
        { t: "tp", v: "Val", s: ">" },
        { t: "pl", v: "  " },
        { t: "tp", v: "Vector", s: "<" },
        { t: "tp", v: "Val", s: "> Multiply(" },
        { t: "kw", v: "    const", s: " " },
        { t: "tp", v: "Matrix", s: "<" },
        { t: "tp", v: "Val", s: ">& A," },
        { t: "kw", v: "    const", s: " " },
        { t: "tp", v: "Vector", s: "<" },
        { t: "tp", v: "Val", s: ">& x) " },
        { t: "kw", v: "const", s: ";" },
        { t: "pl", v: "};" },
        { t: "pl", v: "" },
        { t: "cm", v: "// Suppress implicit instantiation." },
        { t: "kw", v: "extern template", s: " " },
        { t: "tp", v: "Vector", s: "<" },
        { t: "kw", v: "double", s: ">" },
        { t: "tp", v: "  MultiplyEngine", s: "::Multiply(" },
        { t: "kw", v: "    const", s: " " },
        { t: "tp", v: "Matrix", s: "<" },
        { t: "kw", v: "double", s: ">& A," },
        { t: "kw", v: "    const", s: " " },
        { t: "tp", v: "Vector", s: "<" },
        { t: "kw", v: "double", s: ">& x) " },
        { t: "kw", v: "const", s: ";" },
        { t: "kw", v: "extern template", s: " " },
        { t: "tp", v: "Vector", s: "<std::complex<" },
        { t: "kw", v: "double", s: ">>" },
        { t: "tp", v: "  MultiplyEngine", s: "::Multiply(" },
        { t: "kw", v: "    const", s: " " },
        { t: "tp", v: "Matrix", s: "<std::complex<" },
        { t: "kw", v: "double", s: ">>& A," },
        { t: "kw", v: "    const", s: " " },
        { t: "tp", v: "Vector", s: "<std::complex<" },
        { t: "kw", v: "double", s: ">>& x) " },
        { t: "kw", v: "const", s: ";" },
      ],
      ccode: [
        { t: "cm", v: "// MultiplyEngine.cc" },
        { t: "pp", v: '#include "MultiplyEngine.h"' },
        { t: "pl", v: "" },
        { t: "kw", v: "template", s: " <" },
        { t: "kw", v: "class", s: " " },
        { t: "tp", v: "Val", s: ">" },
        { t: "tp", v: "Vector", s: "<" },
        { t: "tp", v: "Val", s: ">" },
        { t: "tp", v: "MultiplyEngine", s: "::Multiply(" },
        { t: "kw", v: "    const", s: " " },
        { t: "tp", v: "Matrix", s: "<" },
        { t: "tp", v: "Val", s: ">& A," },
        { t: "kw", v: "    const", s: " " },
        { t: "tp", v: "Vector", s: "<" },
        { t: "tp", v: "Val", s: ">& x) " },
        { t: "kw", v: "const", s: " {" },
        { t: "tp", v: "  Vector", s: "<" },
        { t: "tp", v: "Val", s: "> y(A.rows());" },
        { t: "kw", v: "  for", s: " (" },
        { t: "tp", v: "size_t", s: " i = " },
        { t: "nm", v: "0", s: "; i < A.rows(); ++i) {" },
        { t: "tp", v: "    Val", s: " sum{" },
        { t: "nm", v: "0", s: "};" },
        { t: "kw", v: "    for", s: " (" },
        { t: "tp", v: "size_t", s: " j = " },
        { t: "nm", v: "0", s: "; j < A.cols(); ++j)" },
        { t: "pl", v: "      sum += A(i,j) * x[j];" },
        { t: "pl", v: "    y[i] = sum;" },
        { t: "pl", v: "  }" },
        { t: "kw", v: "  return", s: " y;" },
        { t: "pl", v: "}" },
        { t: "pl", v: "" },
        { t: "cm", v: "// Compiled once for each type." },
        { t: "kw", v: "template", s: " " },
        { t: "tp", v: "Vector", s: "<" },
        { t: "kw", v: "double", s: ">" },
        { t: "tp", v: "MultiplyEngine", s: "::Multiply(" },
        { t: "kw", v: "  const", s: " " },
        { t: "tp", v: "Matrix", s: "<" },
        { t: "kw", v: "double", s: ">& A," },
        { t: "kw", v: "  const", s: " " },
        { t: "tp", v: "Vector", s: "<" },
        { t: "kw", v: "double", s: ">& x) " },
        { t: "kw", v: "const", s: ";" },
        { t: "kw", v: "template", s: " " },
        { t: "tp", v: "Vector", s: "<std::complex<" },
        { t: "kw", v: "double", s: ">>" },
        { t: "tp", v: "MultiplyEngine", s: "::Multiply(" },
        { t: "kw", v: "  const", s: " " },
        { t: "tp", v: "Matrix", s: "<std::complex<" },
        { t: "kw", v: "double", s: ">>& A," },
        { t: "kw", v: "  const", s: " " },
        { t: "tp", v: "Vector", s: "<std::complex<" },
        { t: "kw", v: "double", s: ">>& x) " },
        { t: "kw", v: "const", s: ";" },
      ],
    },
  },
];

// ─── Render helpers ──────────────────────────────────────────────────────────

function renderTokens(tokens) {
  return tokens.map((tok, i) => (
    <span key={i}>
      <span className={tok.t}>{tok.v}</span>
      {tok.s !== undefined ? <span className="pl">{tok.s}</span> : null}
      {tok.nl !== false && i === tokens.length - 1 ? null : null}
    </span>
  ));
}

function CodeBlock({ tokens }) {
  // Group flat token array into lines by splitting on newline-value tokens
  const lines = [];
  let current = [];
  tokens.forEach((tok) => {
    if (tok.v === "" && tok.t === "pl" && tok.s === undefined) {
      lines.push(current);
      current = [];
      lines.push(null); // blank line
    } else {
      current.push(tok);
    }
  });
  if (current.length) lines.push(current);

  return (
    <pre>
      {lines.map((line, i) =>
        line === null ? (
          <div key={i}>&nbsp;</div>
        ) : (
          <div key={i}>{renderTokens(line)}</div>
        )
      )}
    </pre>
  );
}

// ─── Main component ──────────────────────────────────────────────────────────

export default function App() {
  const [activeTab, setActiveTab] = useState("matrix");
  const section = SECTIONS.find((s) => s.id === activeTab);

  return (
    <>
      <style>{STYLES}</style>
      <div className="root">
        <div className="header">
          <h1>Template Strategy Comparison</h1>
          <p>header-only templates &nbsp;vs&nbsp; explicit instantiation</p>
        </div>

        <div className="tabs">
          {SECTIONS.map((s) => (
            <button
              key={s.id}
              className={`tab ${activeTab === s.id ? "active" : ""}`}
              onClick={() => setActiveTab(s.id)}
            >
              {s.label}
            </button>
          ))}
        </div>

        <div className="columns">
          {/* LEFT — current template style */}
          <div className="panel">
            <div className="panel-header">
              <span className="panel-title tpl">Current — Header Template</span>
              <span className="panel-filename">{section.tpl.file}</span>
            </div>
            <div className="panel-body">
              <CodeBlock tokens={section.tpl.code} />
            </div>
            <div className="note-box">
              <p><strong>Note: </strong>{section.tpl.note}</p>
            </div>
          </div>

          {/* RIGHT — explicit instantiation */}
          <div className="panel">
            <div className="panel-header">
              <span className="panel-title exp">Proposed — Explicit Instantiation</span>
              <span className="panel-filename">{section.exp.file}</span>
            </div>

            <div className="section-label">header</div>
            <div className="panel-body">
              <CodeBlock tokens={section.exp.hcode} />
            </div>

            <div className="divider" />

            <div className="section-label">implementation</div>
            <div className="panel-body">
              <CodeBlock tokens={section.exp.ccode} />
            </div>

            <div className="note-box">
              <p><strong>Note: </strong>{section.exp.note}</p>
            </div>
          </div>
        </div>

        {/* Summary */}
        <div className="summary">
          <div className="summary-panel tpl">
            <h3>Header Template — tradeoffs</h3>
            <ul>
              <li>Implementation exposed in every header — <em>no .cc files</em></li>
              <li>Implicitly instantiated in every translation unit that includes the header</li>
              <li>Recompiled for every scalar type in every .cc — <em>slowest builds</em></li>
              <li>Adding a new scalar type requires <em>no changes</em> — just use it</li>
              <li>Implementation changes force recompilation of all dependents</li>
            </ul>
          </div>
          <div className="summary-panel exp">
            <h3>Explicit Instantiation — tradeoffs</h3>
            <ul>
              <li>Headers hold declarations only — <em>clean, fast to parse</em></li>
              <li><code>extern template</code> suppresses implicit instantiation everywhere else</li>
              <li>Each type compiled <em>exactly once</em> in its .cc file</li>
              <li>Adding a new scalar type: <em>one line per .cc file</em></li>
              <li>Implementation changes recompile only the affected .cc</li>
            </ul>
          </div>
        </div>
      </div>
    </>
  );
}
