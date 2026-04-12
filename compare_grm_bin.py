#!/usr/bin/env python3

from __future__ import annotations

import argparse
import math
import sys
from array import array
from dataclasses import dataclass
from pathlib import Path


@dataclass
class DifferenceRecord:
    index: int
    row: int
    col: int
    value1: float
    value2: float
    abs_diff: float
    rel_diff: float


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Compare two lower-triangular GRM binary files (*.grm.bin).")
    parser.add_argument("file1", type=Path, help="First .grm.bin file")
    parser.add_argument("file2", type=Path, help="Second .grm.bin file")
    parser.add_argument(
        "--atol",
        type=float,
        default=0.0,
        help="Absolute tolerance for math.isclose (default: 0.0)")
    parser.add_argument(
        "--rtol",
        type=float,
        default=0.0,
        help="Relative tolerance for math.isclose (default: 0.0)")
    parser.add_argument(
        "--limit",
        type=int,
        default=10,
        help="How many differing entries to print (default: 10)")
    return parser.parse_args()


def read_double_array(path: Path) -> array:
    if not path.is_file():
        raise FileNotFoundError(f"File does not exist: {path}")

    values = array("d")
    with path.open("rb") as fin:
        values.frombytes(fin.read())

    if sys.byteorder != "little":
        values.byteswap()
    return values


def infer_num_ids(num_entries: int) -> int:
    discriminant = 1 + 8 * num_entries
    root = math.isqrt(discriminant)
    if root * root != discriminant:
        raise ValueError(
            f"Entry count {num_entries} is not a valid lower-triangular matrix length.")

    num_ids = (root - 1) // 2
    if num_ids * (num_ids + 1) // 2 != num_entries:
        raise ValueError(
            f"Entry count {num_entries} does not map to n(n+1)/2 exactly.")
    return num_ids


def tri_index_to_row_col(index: int) -> tuple[int, int]:
    row = (math.isqrt(8 * index + 1) - 1) // 2
    while row * (row + 1) // 2 > index:
        row -= 1
    while (row + 1) * (row + 2) // 2 <= index:
        row += 1
    col = index - row * (row + 1) // 2
    return row, col


def relative_difference(value1: float, value2: float) -> float:
    denominator = max(abs(value1), abs(value2))
    if denominator == 0.0:
        return 0.0
    return abs(value1 - value2) / denominator


def compare_values(values1: array, values2: array, atol: float, rtol: float,
        limit: int) -> tuple[int, int, int | None, DifferenceRecord | None, list[DifferenceRecord]]:
    num_diff = 0
    num_outside_tol = 0
    first_diff_index: int | None = None
    max_abs_record: DifferenceRecord | None = None
    samples: list[DifferenceRecord] = []

    for index, (value1, value2) in enumerate(zip(values1, values2)):
        if value1 == value2:
            continue

        row, col = tri_index_to_row_col(index)
        abs_diff = abs(value1 - value2)
        rel_diff = relative_difference(value1, value2)
        record = DifferenceRecord(index, row, col, value1, value2, abs_diff, rel_diff)

        num_diff += 1
        if not math.isclose(value1, value2, rel_tol=rtol, abs_tol=atol):
            num_outside_tol += 1
        if first_diff_index is None:
            first_diff_index = index
        if max_abs_record is None or abs_diff > max_abs_record.abs_diff:
            max_abs_record = record
        if len(samples) < limit:
            samples.append(record)

    return num_diff, num_outside_tol, first_diff_index, max_abs_record, samples


def print_summary(path1: Path, path2: Path, values1: array, values2: array,
        atol: float, rtol: float, num_diff: int, num_outside_tol: int,
        first_diff_index: int | None, max_abs_record: DifferenceRecord | None,
        samples: list[DifferenceRecord]) -> None:
    num_entries = len(values1)
    num_ids = infer_num_ids(num_entries)

    print(f"File 1: {path1}")
    print(f"File 2: {path2}")
    print(f"Entries: {num_entries}")
    print(f"Inferred matrix size: {num_ids} x {num_ids}")
    print(f"Exact equal: {'yes' if num_diff == 0 else 'no'}")
    print(
        f"Within tolerance (atol={atol}, rtol={rtol}): "
        f"{'yes' if num_outside_tol == 0 else 'no'}")
    print(f"Differing entries: {num_diff}")
    print(f"Differing entries outside tolerance: {num_outside_tol}")

    if first_diff_index is not None:
        first_row, first_col = tri_index_to_row_col(first_diff_index)
        print(
            f"First difference: index={first_diff_index}, row={first_row}, col={first_col}, "
            f"byte_offset={first_diff_index * 8}")

    if max_abs_record is not None:
        print(
            f"Max abs diff: {max_abs_record.abs_diff:.18g} "
            f"(index={max_abs_record.index}, row={max_abs_record.row}, col={max_abs_record.col})")

    if samples:
        print("\nSample differences:")
        for record in samples:
            print(
                f"  index={record.index:>8} row={record.row:>5} col={record.col:>5} "
                f"value1={record.value1:.18g} value2={record.value2:.18g} "
                f"abs_diff={record.abs_diff:.18g} rel_diff={record.rel_diff:.18g}")


def main() -> int:
    args = parse_args()
    values1 = read_double_array(args.file1)
    values2 = read_double_array(args.file2)

    if len(values1) != len(values2):
        print(
            f"Length mismatch: {args.file1} has {len(values1)} doubles, "
            f"{args.file2} has {len(values2)} doubles.",
            file=sys.stderr)
        return 1

    num_diff, num_outside_tol, first_diff_index, max_abs_record, samples = compare_values(
        values1, values2, args.atol, args.rtol, args.limit)

    print_summary(
        args.file1,
        args.file2,
        values1,
        values2,
        args.atol,
        args.rtol,
        num_diff,
        num_outside_tol,
        first_diff_index,
        max_abs_record,
        samples)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
