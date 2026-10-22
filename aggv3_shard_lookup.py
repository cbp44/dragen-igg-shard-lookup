#!/usr/bin/env python3
"""
aggv3_shard_lookup.py
=====================

Convert genomic coordinates (GRCh38) into Genomics England aggV3 shard and
sub-shard numbers. This is a command-line version of the JavaScript used
by the aggV3 online shard lookup tool, and it uses the same overlap logic.

Accepted coordinate formats
---------------------------
    1:242434                  single position, without the 'chr' prefix
    chr1:242434               single position, with the 'chr' prefix
    chr1:242434-537588484     a range (may cover several shards)
    chrX:1,000,000-2,000,000  thousands separators are allowed
    MT:100 / chrM:100         mitochondrial genome

Chromosome names are case-insensitive. If a range is given backwards
(start > end), the two ends are swapped, as the JavaScript tool does.

Output
------
By default the script prints a comma-separated list of the shards that
overlap the query, with duplicates removed, e.g.

    $ ./aggv3_shard_lookup.py chr1:242434-537588484
    1,2,3,4,5,6,7,8

If you give several coordinates, each output line starts with the query
and a tab.

With --subshard the script prints a tab-separated table with one row per
overlapping sub-shard:

    query  chrom  subshard_start  subshard_end  shard  subshard

Shard table
-----------
The shard/sub-shard boundaries are read from a JSON file:
    * the path given with --shard-file, or, if that is not given,
    * 'aggv3_shards.json' in the same folder as this script.

The JSON file is the `shardsByChr` object from the JavaScript, i.e. the
text from the first '{' to the last '}':

    {"chr1": [[10061, 1111562, 1, 1], ...], "chr2": [...], ...}

Each entry is [start, end, shard, subshard], with inclusive coordinates.

Exit codes
----------
    0  success: every query overlapped at least one sub-shard
    1  at least one query did not overlap any sub-shard (e.g. it falls
       in a centromere/gap or between two sub-shards)
    2  invalid command-line arguments or badly formed coordinates
    3  the shard table file is missing or invalid
"""

from __future__ import annotations

import argparse
import bisect
import json
import re
import sys
from pathlib import Path
from typing import Dict, List, NamedTuple, Optional, Sequence

__version__ = "1.1.0"

EXIT_OK = 0
EXIT_NO_HITS = 1
EXIT_USAGE = 2
EXIT_DATA = 3

DEFAULT_SHARD_FILENAME = "aggv3_shards.json"
DEFAULT_SHARD_PATH = Path(__file__).resolve().parent / DEFAULT_SHARD_FILENAME


# --------------------------------------------------------------------------- #
# Data structures
# --------------------------------------------------------------------------- #
class Region(NamedTuple):
    """One sub-shard: [start, end] (inclusive) on a chromosome."""
    start: int
    end: int
    shard: int
    subshard: int


class Query(NamedTuple):
    """A parsed, validated coordinate query."""
    raw: str      # the string as the user typed it
    chrom: str    # normalised chromosome name, e.g. 'chr1', 'chrX', 'chrM'
    start: int
    end: int

    @property
    def is_range(self) -> bool:
        return self.start != self.end


class Hit(NamedTuple):
    """A sub-shard that overlaps a query."""
    chrom: str
    start: int
    end: int
    shard: int
    subshard: int


class CoordinateError(ValueError):
    """Raised when a coordinate string cannot be parsed or is invalid."""


class ShardTableError(Exception):
    """Raised when the shard table file is missing or malformed."""


# --------------------------------------------------------------------------- #
# Chromosome name handling
# --------------------------------------------------------------------------- #
def normalise_chrom(token: str) -> str:
    """
    Convert a chromosome name to the 'chrN' form used in the shard table.

    Examples: '1' -> 'chr1', 'CHR1' -> 'chr1', 'x' -> 'chrX',
              'chrMT' -> 'chrM', '01' -> 'chr1'.
    """
    name = re.sub(r"^chr", "", token.strip(), flags=re.IGNORECASE)
    if not name:
        raise CoordinateError("missing chromosome name")
    name = name.upper()
    if name == "MT":
        name = "M"
    if name.isdigit():
        name = str(int(name))  # strip leading zeros: '01' -> '1'
    return "chr" + name


# --------------------------------------------------------------------------- #
# Shard table
# --------------------------------------------------------------------------- #
class ShardTable:
    """Sub-shard intervals for each chromosome, sorted by start position."""

    def __init__(self, regions_by_chrom: Dict[str, List[Region]]) -> None:
        self._regions = regions_by_chrom
        # Precomputed start positions so binary search is fast (bisect).
        self._starts = {c: [r.start for r in regs]
                        for c, regs in regions_by_chrom.items()}

    @property
    def chromosomes(self) -> List[str]:
        return list(self._regions)

    def __contains__(self, chrom: str) -> bool:
        return chrom in self._regions

    def find_overlaps(self, chrom: str, q_start: int, q_end: int) -> List[Hit]:
        """
        Return all sub-shards that overlap [q_start, q_end] (inclusive).

        Same steps as the JavaScript `findOverlaps`:
          1. binary-search for the first region whose start >= q_start
          2. step back one, because the region before may span q_start
          3. skip regions that end before q_start
          4. collect regions until one starts after q_end
        """
        regions = self._regions.get(chrom)
        if not regions:
            return []

        i = bisect.bisect_left(self._starts[chrom], q_start)  # lowerBoundByStart
        if i > 0:
            i -= 1
        n = len(regions)
        while i < n and regions[i].end < q_start:
            i += 1

        hits: List[Hit] = []
        while i < n and regions[i].start <= q_end:
            r = regions[i]
            if r.end >= q_start:
                hits.append(Hit(chrom, r.start, r.end, r.shard, r.subshard))
            i += 1
        return hits


def resolve_shard_file(cli_value: Optional[str]) -> Path:
    """Use --shard-file if given, otherwise aggv3_shards.json next to the script."""
    if cli_value:
        return Path(cli_value).expanduser()
    return DEFAULT_SHARD_PATH


def load_shard_table(path: Path) -> ShardTable:
    """Load and validate the shard table from a JSON file."""
    if not path.is_file():
        raise ShardTableError(
            f"shard table JSON file not found: {path}\n"
            f"Either create '{DEFAULT_SHARD_FILENAME}' in the same folder as "
            f"this script by pasting the JavaScript `shardsByChr` object "
            f"into it, or give its location with --shard-file PATH."
        )

    try:
        text = path.read_text(encoding="utf-8")
    except OSError as exc:
        raise ShardTableError(f"cannot read {path}: {exc}") from exc

    try:
        raw = json.loads(text)
    except json.JSONDecodeError as exc:
        hint = ""
        if "shardsByChr" in text:
            hint = ("\nHint: the file seems to contain JavaScript. Keep only "
                    "the object itself: delete 'const shardsByChr =' at the "
                    "start and the ';' at the end.")
        raise ShardTableError(f"{path}: invalid JSON ({exc}){hint}") from exc

    if not isinstance(raw, dict) or not raw:
        raise ShardTableError(
            f"{path}: expected a non-empty JSON object of chromosomes, "
            f'e.g. {{"chr1": [[start, end, shard, subshard], ...], ...}}'
        )

    regions_by_chrom: Dict[str, List[Region]] = {}
    for chrom_key, entries in raw.items():
        if not isinstance(entries, list):
            raise ShardTableError(f"{path}: value for {chrom_key!r} is not a list")
        regions: List[Region] = []
        for n, entry in enumerate(entries):
            if (not isinstance(entry, list) or len(entry) != 4
                    or not all(isinstance(v, int) and not isinstance(v, bool)
                               for v in entry)):
                raise ShardTableError(
                    f"{path}: entry {n} of {chrom_key!r} is not "
                    f"[start, end, shard, subshard]: {entry!r}"
                )
            region = Region(*entry)
            if region.start > region.end:
                raise ShardTableError(
                    f"{path}: entry {n} of {chrom_key!r} has start > end: {entry!r}"
                )
            regions.append(region)
        regions.sort()  # binary search needs regions sorted by start
        try:
            regions_by_chrom[normalise_chrom(chrom_key)] = regions
        except CoordinateError:
            raise ShardTableError(
                f"{path}: invalid chromosome key {chrom_key!r}"
            ) from None

    return ShardTable(regions_by_chrom)


# --------------------------------------------------------------------------- #
# Coordinate parsing
# --------------------------------------------------------------------------- #
_POSITION_RE = re.compile(r"(\d[\d,]*)(?:\s*-\s*(\d[\d,]*))?")


def _to_int(token: str) -> int:
    digits = token.replace(",", "")
    if not digits.isdigit():
        raise CoordinateError(f"invalid number {token!r}")
    return int(digits)


def parse_coordinate(raw: str, table: ShardTable,
                     warn: bool = True) -> Query:
    """
    Parse 'chr1:15000', '1:15000' or 'chr1:15000-20000' into a Query.

    Raises CoordinateError with a clear message if the input is invalid.
    """
    text = raw.strip()
    if not text:
        raise CoordinateError("empty coordinate")

    if text.count(":") != 1:
        raise CoordinateError(
            f"{raw!r}: expected the form CHROM:POS or CHROM:START-END "
            f"(e.g. chr1:15000 or chr1:15000-20000)"
        )

    chrom_part, pos_part = text.split(":")
    try:
        chrom = normalise_chrom(chrom_part)
    except CoordinateError as exc:
        raise CoordinateError(f"{raw!r}: {exc}") from None

    if chrom not in table:
        valid = ", ".join(c.replace("chr", "") for c in table.chromosomes)
        raise CoordinateError(
            f"{raw!r}: unknown chromosome {chrom_part.strip()!r} "
            f"(valid chromosomes: {valid})"
        )

    pos_part = pos_part.strip()
    m = _POSITION_RE.fullmatch(pos_part)
    if not m:
        raise CoordinateError(
            f"{raw!r}: invalid position {pos_part!r}; positions must be "
            f"non-negative integers, e.g. 15000 or 15000-20000"
        )

    start = _to_int(m.group(1))
    end = _to_int(m.group(2)) if m.group(2) is not None else start

    if start > end:  # the JavaScript tool also swaps the ends
        if warn:
            print(f"warning: {raw!r}: start > end, swapping to "
                  f"{chrom}:{end}-{start}", file=sys.stderr)
        start, end = end, start

    return Query(raw=text, chrom=chrom, start=start, end=end)


# --------------------------------------------------------------------------- #
# Output formatting
# --------------------------------------------------------------------------- #
def unique_shards(hits: Sequence[Hit]) -> List[int]:
    """Distinct shard numbers, in ascending order."""
    return sorted({h.shard for h in hits})


def format_shard_list(hits: Sequence[Hit]) -> str:
    return ",".join(str(s) for s in unique_shards(hits))


# --------------------------------------------------------------------------- #
# Command-line interface
# --------------------------------------------------------------------------- #
def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="aggv3_shard_lookup.py",
        description="Find the Genomics England aggV3 shard (and optionally "
                    "sub-shard) numbers for GRCh38 genomic coordinates.",
        epilog="examples:\n"
               "  %(prog)s chr1:242434\n"
               "  %(prog)s 1:242434\n"
               "  %(prog)s chr1:242434-537588484\n"
               "  %(prog)s chr1:242434-537588484 --subshard\n"
               "  %(prog)s chr2:100000 chrX:5,000,000 --shard-file /path/to/aggv3_shards.json",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "coordinates", nargs="+", metavar="COORD",
        help="genomic coordinate(s): CHROM:POS or CHROM:START-END; the "
             "'chr' prefix is optional (e.g. 1:242434, chr1:242434, "
             "chr1:242434-537588484)",
    )
    parser.add_argument(
        "-s", "--subshard", action="store_true",
        help="also report sub-shard numbers and boundaries "
             "(tab-separated table, one row per sub-shard)",
    )
    parser.add_argument(
        "--no-header", action="store_true",
        help="with --subshard, leave out the table header line",
    )
    parser.add_argument(
        "-f", "--shard-file", metavar="PATH",
        help=f"path to the shard table JSON file "
             f"(default: {DEFAULT_SHARD_FILENAME} in the script's folder, "
             f"i.e. {DEFAULT_SHARD_PATH})",
    )
    parser.add_argument(
        "-q", "--quiet", action="store_true",
        help="suppress warnings and 'no shard found' messages on stderr",
    )
    parser.add_argument("--version", action="version",
                        version=f"%(prog)s {__version__}")
    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)

    # Load the shard table.
    try:
        table = load_shard_table(resolve_shard_file(args.shard_file))
    except ShardTableError as exc:
        print(f"error: {exc}", file=sys.stderr)
        return EXIT_DATA

    # Check every coordinate before printing anything.
    queries: List[Query] = []
    errors: List[str] = []
    for raw in args.coordinates:
        try:
            queries.append(parse_coordinate(raw, table, warn=not args.quiet))
        except CoordinateError as exc:
            errors.append(str(exc))
    if errors:
        parser.error("invalid coordinate(s):\n  " + "\n  ".join(errors))
        # parser.error() exits with status 2

    multi = len(queries) > 1
    exit_code = EXIT_OK

    if args.subshard and not args.no_header:
        print("query\tchrom\tsubshard_start\tsubshard_end\tshard\tsubshard")

    for q in queries:
        hits = table.find_overlaps(q.chrom, q.start, q.end)

        if not hits:
            exit_code = EXIT_NO_HITS
            if not args.quiet:
                where = f"{q.chrom}:{q.start}" + (f"-{q.end}" if q.is_range else "")
                print(f"warning: no shard found for {where} (it may fall in a "
                      f"gap, e.g. a centromere or between sub-shards)",
                      file=sys.stderr)
            if args.subshard:
                print(f"{q.raw}\t{q.chrom}\tNA\tNA\tNA\tNA")
            elif multi:
                print(f"{q.raw}\tNA")
            continue

        if args.subshard:
            for h in hits:
                print(f"{q.raw}\t{h.chrom}\t{h.start}\t{h.end}\t"
                      f"{h.shard}\t{h.subshard}")
        elif multi:
            print(f"{q.raw}\t{format_shard_list(hits)}")
        else:
            print(format_shard_list(hits))

    return exit_code


if __name__ == "__main__":
    sys.exit(main())