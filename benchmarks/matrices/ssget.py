#!/usr/bin/env python3
"""Choose and fetch matrices from the SuiteSparse Matrix Collection.

Two modes, deliberately separate, because one of them chooses and the other pays:

    ./ssget.py list                      print the matrices that match the filter
    ./ssget.py list --max-n 200000       the same, with a different range
    ./ssget.py list --per-kind 2 --max-nnz 500000     a spread across kinds, cheap end only
    ./ssget.py list --per-kind 6 --values             only matrices that carry values
    ./ssget.py fetch candidates.txt      download the ones named in a file
    ./ssget.py fetch --all               download everything the filter matches

`list` reads one small index file and downloads no matrix at all, so it is cheap to re-run
until the candidate list looks right. `fetch` takes that list and populates data/, and it is
where the bytes are: the index carries nnz for every matrix, so the cost of a set is visible
before any of it moves. The output of `list` is valid input to `fetch` unedited, since fetch
reads the first whitespace token of each line and ignores the rest, so the loop is print,
delete the rows we do not want, fetch.

THE FILTER IS DELIBERATELY NARROW, and it matches what the reader beside this script accepts:
real values, square, and numerically symmetric, which for this collection means the Matrix
Market file stores one triangle under a `symmetric` header and is symmetric both structurally
and numerically by construction. A matrix stored `general` whose values happen to be symmetric
is excluded, and so is anything complex or pattern-only. Widening any of that is a change here
and a change in the reader together.

BINARY MEANS PATTERN HERE, which is what `--values` is for. A binary matrix has values that are
all 1, so the collection stores its pattern and nothing else, and the file's Matrix Market field
type comes out `pattern` with two indices per line and no third column. The index's `isReal` flag
does NOT predict this: it means "not complex". The index's `isBinary` flag does, exactly, on the
evidence of two independent samples, 19 of 19 and then 41 of 41 pattern files matching the binary
rows and no others. So `--values` excludes binary matrices and what comes back carries numbers.
Leave it off when a pattern is all a study needs, which is true of any ordering question.

WHAT IS NOT FILTERED, because the index cannot say it: whether the diagonal is structurally
present. It need not be. The conversion the driver uses inserts a structural zero on any
diagonal nothing landed on, which is what the symbolic factorization needs and is ordinary
input for LDL, though Cholesky will refuse it.

Layout: matrices land in data/<Group>/<name>.mtx, the collection's own layout, so a path says
where the file came from and two groups may share a matrix name without colliding. The index
itself is cached at data/ssstats.csv and re-read from there; --refresh downloads it again.

The collection lives at https://sparse.tamu.edu. Its own interfaces page notes that the site
is https while the data files sit on an http repository, which some browsers block; command
line downloads are not affected. If the collection ever moves, SITE below is the one constant
to change.
"""

import argparse
import io
import os
import sys
import tarfile
import urllib.request

SITE = "https://sparse.tamu.edu"
INDEX_URL = SITE + "/files/ssstats.csv"
MATRIX_URL = SITE + "/MM/{group}/{name}.tar.gz"

TIMEOUT = 60


def repo_data_dir():
    """data/ at the repo root, this script living two levels below it."""
    here = os.path.dirname(os.path.abspath(__file__))
    return os.path.join(here, "..", "..", "data")


# ----------------------------------------------------------------------------------------------
# The index
# ----------------------------------------------------------------------------------------------

# ssstats.csv is one line per matrix behind a two-line preamble, a matrix count and a version
# stamp. The fields are positional:
#
#   0 group  1 name  2 rows  3 cols  4 nnz  5 isReal  6 isBinary  7 isND  8 posdef
#   9 pattern symmetry  10 numerical symmetry  11 kind  [12 an extra field in newer versions]
#
# We take the first eleven positionally and treat everything after as the kind, dropping a
# trailing purely numeric field if one is there. That way a kind containing a comma, or another
# column appended at the end later, costs nothing. The preamble is skipped by trying to parse
# rather than by counting lines, so a change to its shape does not silently drop a matrix.


class Matrix:
    def __init__(self, group, name, rows, cols, nnz, is_real, is_binary, posdef, nsym, kind):
        self.group = group
        self.name = name
        self.rows = rows
        self.cols = cols
        self.nnz = nnz
        self.is_real = is_real
        self.is_binary = is_binary
        self.posdef = posdef
        self.nsym = nsym
        self.kind = kind

    @property
    def full_name(self):
        return self.group + "/" + self.name


def parse_index_line(line):
    """One matrix, or None if the line is preamble or unparsable."""
    fields = line.rstrip("\n").split(",")
    if len(fields) < 12:
        return None

    try:
        rows = int(fields[2])
        cols = int(fields[3])
        nnz = int(fields[4])
        is_real = int(fields[5])
        is_binary = int(fields[6])
        posdef = int(fields[8])
        nsym = float(fields[10])
    except ValueError:
        return None

    rest = fields[11:]
    if len(rest) > 1 and rest[-1].strip().isdigit():
        rest = rest[:-1]
    kind = ",".join(rest).strip()

    return Matrix(fields[0].strip(), fields[1].strip(), rows, cols, nnz,
                  is_real, is_binary, posdef, nsym, kind)


def load_index(data_dir, refresh):
    """The cached index, downloading it if it is missing or if refresh was asked for."""
    path = os.path.join(data_dir, "ssstats.csv")

    if refresh or not os.path.exists(path):
        os.makedirs(data_dir, exist_ok=True)
        sys.stderr.write("downloading index: " + INDEX_URL + "\n")
        with urllib.request.urlopen(INDEX_URL, timeout=TIMEOUT) as response:
            body = response.read()
        with open(path, "wb") as out:
            out.write(body)

    matrices = []
    with open(path, "r", encoding="utf-8", errors="replace") as index:
        for line in index:
            matrix = parse_index_line(line)
            if matrix is not None:
                matrices.append(matrix)

    if not matrices:
        sys.exit("index at " + path + " parsed to nothing; try --refresh")

    return matrices


def select(matrices, min_n, max_n, max_nnz, per_kind, values_only):
    """Real, square, numerically symmetric, in the size range, and optionally sampled by kind."""
    chosen = []
    for m in matrices:
        if not m.is_real:
            continue
        if m.rows != m.cols:
            continue
        if m.nsym != 1.0:
            continue
        if m.rows < min_n or m.rows > max_n:
            continue
        if max_nnz is not None and m.nnz > max_nnz:
            continue
        if values_only and m.is_binary:
            continue
        chosen.append(m)

    chosen.sort(key=lambda m: (m.nnz, m.full_name))

    if per_kind is not None:
        chosen = sample_by_kind(chosen, per_kind)

    return chosen


def sample_by_kind(matrices, per_kind):
    """The `per_kind` cheapest of each kind, the list arriving sorted by nnz.

    The kind column is the only account of structure the index gives us, and structural variety
    is the whole reason for going to real matrices at all: a hundred matrices that are all
    `structural problem` would be one family again in a new costume.

    Taking the CHEAPEST of each kind biases toward the small end within every kind, which is
    what a correctness survey wants and is not what a timing ladder will want later. Spreading
    within a kind is a change here when that day comes.
    """
    taken = {}
    chosen = []
    for m in matrices:
        count = taken.get(m.kind, 0)
        if count < per_kind:
            taken[m.kind] = count + 1
            chosen.append(m)
    return chosen


# ----------------------------------------------------------------------------------------------
# list
# ----------------------------------------------------------------------------------------------


def do_list(matrices):
    width = max([len(m.full_name) for m in matrices] + [len("matrix")])

    print("# %-*s %9s %11s %7s %7s  %s"
          % (width, "matrix", "n", "nnz", "posdef", "binary", "kind"))
    for m in matrices:
        print("  %-*s %9d %11d %7d %7d  %s"
              % (width, m.full_name, m.rows, m.nnz, m.posdef, m.is_binary, m.kind))

    total = sum(m.nnz for m in matrices)
    kinds = len({m.kind for m in matrices})
    print("# %d matrices, %d nonzeros in total, %d kinds" % (len(matrices), total, kinds))


# ----------------------------------------------------------------------------------------------
# fetch
# ----------------------------------------------------------------------------------------------


def read_wanted(path):
    """Group/Name from the first token of each line, so `list` output feeds straight back in."""
    wanted = []
    with open(path, "r", encoding="utf-8") as listing:
        for line in listing:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            wanted.append(line.split()[0])
    return wanted


def header_line(path):
    """The Matrix Market banner, reported so a header we cannot read shows up at fetch time."""
    with open(path, "r", encoding="utf-8", errors="replace") as source:
        return source.readline().strip()


def fetch_one(matrix, data_dir):
    """Download one tarball and write its .mtx to data/<Group>/<name>.mtx. Returns the path."""
    group_dir = os.path.join(data_dir, matrix.group)
    target = os.path.join(group_dir, matrix.name + ".mtx")

    if os.path.exists(target):
        return target, True

    url = MATRIX_URL.format(group=matrix.group, name=matrix.name)
    with urllib.request.urlopen(url, timeout=TIMEOUT) as response:
        payload = response.read()

    # The tarball holds <name>/<name>.mtx, and for some matrices a right-hand side or a set of
    # coordinates beside it. We want the one named for the matrix. Members are read and written
    # by hand rather than through extractall, which would honor whatever path the archive
    # carries.
    wanted = matrix.name + ".mtx"
    with tarfile.open(fileobj=io.BytesIO(payload), mode="r:gz") as archive:
        member = None
        for candidate in archive.getmembers():
            if candidate.isfile() and os.path.basename(candidate.name) == wanted:
                member = candidate
                break
        if member is None:
            raise RuntimeError("no " + wanted + " in the archive")

        os.makedirs(group_dir, exist_ok=True)
        source = archive.extractfile(member)
        with open(target, "wb") as out:
            while True:
                block = source.read(1 << 20)
                if not block:
                    break
                out.write(block)

    return target, False


def do_fetch(matrices, wanted, data_dir):
    by_name = {m.full_name: m for m in matrices}

    if wanted is None:
        selected = matrices
    else:
        selected = []
        for name in wanted:
            if name in by_name:
                selected.append(by_name[name])
            else:
                print("%-40s SKIP  not in the filtered index" % name)

    failed = 0
    cached_count = 0
    got_count = 0

    for m in selected:
        try:
            path, cached = fetch_one(m, data_dir)
        except Exception as error:  # a refusal is information, so report it and keep going
            print("%-40s FAIL  %s" % (m.full_name, error))
            failed += 1
            continue

        if cached:
            cached_count += 1
        else:
            got_count += 1

        print("%-40s %s  %s" % (m.full_name, "have" if cached else "got ", header_line(path)))

    # `have` against `got` is the number worth seeing on a second run: it says how much of a
    # widened filter was already on disk, so growing the set never looks like starting over.
    print("# %d of %d in %s: %d already there, %d downloaded, %d failed"
          % (cached_count + got_count, len(selected), data_dir,
             cached_count, got_count, failed))


# ----------------------------------------------------------------------------------------------


def main():
    parser = argparse.ArgumentParser(
        description="choose and fetch matrices from the SuiteSparse Matrix Collection")
    parser.add_argument("mode", choices=["list", "fetch"])
    parser.add_argument("listing", nargs="?",
                        help="for fetch: a file naming the matrices, as printed by list")
    parser.add_argument("--all", action="store_true",
                        help="for fetch: take everything the filter matches")
    parser.add_argument("--min-n", type=int, default=1000, help="smallest n (default 1000)")
    parser.add_argument("--max-n", type=int, default=100000, help="largest n (default 100000)")
    parser.add_argument("--max-nnz", type=int, default=None,
                        help="largest nnz(A), which bounds both the download and the "
                             "factorization where n bounds neither")
    parser.add_argument("--per-kind", type=int, default=None,
                        help="take only the N cheapest matrices of each kind, for structural "
                             "variety rather than more of one family")
    parser.add_argument("--values", action="store_true",
                        help="exclude binary matrices, which this collection stores WITHOUT "
                             "values; use it when the run needs numbers rather than a pattern")
    parser.add_argument("--data", default=None, help="where matrices go (default the repo's data/)")
    parser.add_argument("--refresh", action="store_true", help="download the index again")
    args = parser.parse_args()

    data_dir = args.data if args.data else repo_data_dir()
    index = load_index(data_dir, args.refresh)
    matrices = select(index, args.min_n, args.max_n, args.max_nnz, args.per_kind,
                      args.values)

    if not matrices:
        sys.exit("nothing matches the filter in that size range")

    if args.mode == "list":
        do_list(matrices)
        return

    if args.listing is None and not args.all:
        sys.exit("fetch needs a listing file, or --all")
    if args.listing is not None and args.all:
        sys.exit("fetch takes a listing file or --all, not both")

    wanted = None if args.all else read_wanted(args.listing)
    do_fetch(matrices, wanted, data_dir)


if __name__ == "__main__":
    main()
