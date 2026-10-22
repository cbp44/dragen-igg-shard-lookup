# dragen-igg-shard-lookup

Command-line utility for mapping GRCh38 genomic coordinates to the corresponding Genomics England aggV3 / DRAGEN iGG shard. It uses the shard boundaries from the Genomics England aggV3 lookup tool and supports single positions and intervals.

## Requirements

- Python 3
- `aggv3_shards.json`, included in this repository

The script finds the JSON file next to itself by default, so no installation or extra Python packages are needed.

## Usage

Run the script with one or more coordinates:

```sh
python3 aggv3_shard_lookup.py chr14:48464880
python3 aggv3_shard_lookup.py 1:242434 chrX:5000000
python3 aggv3_shard_lookup.py chr1:242434-537588484
```

Coordinates use the form `CHROM:POS` or `CHROM:START-END`. The `chr` prefix is optional; chromosome names are case-insensitive. Commas are allowed in positions, and reversed interval endpoints are swapped automatically. Mitochondrial coordinates can use `MT` or `chrM`.

By default, the script prints the distinct overlapping shard numbers, comma-separated. With multiple coordinates, each line is prefixed by the coordinate and a tab. For example:

```text
$ python3 aggv3_shard_lookup.py chr1:242434 chrX:5000000
chr1:242434	1
chrX:5000000	94
```

Use `--subshard` to print one row for each overlapping sub-shard, including its boundaries:

```sh
python3 aggv3_shard_lookup.py chr1:242434 --subshard
```

The table columns are `query`, `chrom`, `subshard_start`, `subshard_end`, `shard`, and `subshard`. Use `--no-header` to omit the column names.

## Options

- `-h`, `--help`: show all options and examples.
- `-s`, `--subshard`: show sub-shard boundaries and numbers.
- `--no-header`: omit the header in `--subshard` output.
- `-f PATH`, `--shard-file PATH`: use a shard JSON file at another path.
- `-q`, `--quiet`: suppress warnings and no-shard messages on standard error.
- `--version`: print the script version.

The script exits with status `0` when every query overlaps a sub-shard, `1` when at least one query has no overlap, `2` for invalid arguments or coordinates, and `3` if the shard table is missing or invalid.
