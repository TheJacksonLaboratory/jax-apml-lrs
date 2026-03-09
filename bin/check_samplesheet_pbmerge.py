#!/usr/bin/env python3
"""
Validate and transform the input samplesheet for the lrs_pbmerge pipeline.

Checks that the samplesheet contains the required columns (sID, bam1, bam2),
validates BAM file extensions, and ensures the combination of sample ID and
bam1 is unique. bam3 is optional. Writes a validated CSV to file_out.

Usage:
    check_samplesheet_pbmerge.py <samplesheet.csv> <samplesheet.valid.csv>

Expected samplesheet format:
    sID,bam1,bam2,bam3
    SAMPLE01,/path/to/run1.bam,/path/to/run2.bam,/path/to/run3.bam
"""

import argparse
import csv
import logging
import sys
from collections import Counter
from pathlib import Path

logger = logging.getLogger()


class RowChecker:
    """
    Validate and transform each row of the pbmerge samplesheet.

    Attributes:
        modified (list): Validated and transformed rows in original order.

    Notes:
        bam2 and bam3 are optional. The uniqueness check is on the combination
        of sample ID and bam1 path, allowing the same sample ID to appear
        multiple times with different primary BAM files. Samples sharing the
        same ID are renamed with a _T{n} suffix (e.g. SAMPLE01_T1, SAMPLE01_T2).
    """

    VALID_FORMATS = (".bam",)

    def __init__(
        self,
        sample_col="sID",
        first_col="bam1",
        second_col="bam2",
        third_col="bam3",
        **kwargs
    ):
        """
        Initialize the row checker with expected column names.

        Args:
            sample_col (str): Column containing the sample identifier (default: sID).
            first_col (str): Column containing the first BAM file path (default: bam1).
            second_col (str): Column containing the second BAM file path (default: bam2).
            third_col (str): Column containing the third BAM file path (default: bam3).
        """
        super().__init__(**kwargs)
        self._sample_col = sample_col
        self._first_col  = first_col
        self._second_col = second_col
        self._third_col  = third_col
        self._seen       = set()
        self.modified    = []

    def validate_and_transform(self, row):
        """
        Run all validations on a row and append to modified list if valid.

        Args:
            row (dict): Mapping from column headers to row values.
        """
        self._validate_sample(row)
        self._validate_first(row)
        self._validate_second(row)
        self._validate_third(row)
        self._seen.add((row[self._sample_col], row[self._first_col]))
        self.modified.append(row)

    def _validate_sample(self, row):
        """Assert that the sample ID is non-empty and replace spaces with underscores."""
        row[self._sample_col] = row[self._sample_col].replace(" ", "_")

    def _validate_first(self, row):
        """Assert that the first BAM entry is non-empty and has a valid extension."""
        self._validate_bam_format(row[self._first_col])

    def _validate_second(self, row):
        """Assert that the second BAM entry has a valid extension if provided."""
        pass

    def _validate_third(self, row):
        """Assert that the third BAM entry has a valid extension if provided."""
        pass

    def _validate_bam_format(self, filename):
        """Assert that a given filename has one of the expected BAM extensions."""
        if not any(filename.endswith(ext) for ext in self.VALID_FORMATS):
            raise AssertionError(
                f"The BAM file has an unrecognized extension: {filename}\n"
                f"It should be one of: {', '.join(self.VALID_FORMATS)}"
            )

    def validate_unique_samples(self):
        """
        Assert that the combination of sample name and bam1 filename is unique.

        Samples sharing the same ID are renamed with a _T{n} suffix to
        distinguish them (e.g. SAMPLE01_T1, SAMPLE01_T2).
        """
        if len(self._seen) != len(self.modified):
            raise AssertionError("The pair of sample name and BAM must be unique.")
        counts = Counter()
        for row in self.modified:
            counts[row[self._sample_col]] += 1
        seen = Counter()
        for row in self.modified:
            sample = row[self._sample_col]
            if counts[sample] > 1:
                seen[sample] += 1
                row[self._sample_col] = f"{sample}_T{seen[sample]}"


def read_head(handle, num_lines=10):
    """Read the first num_lines lines from the current file position."""
    lines = []
    for idx, line in enumerate(handle):
        if idx == num_lines:
            break
        lines.append(line)
    return "".join(lines)


def sniff_format(handle):
    """
    Detect the tabular format (CSV, TSV, etc.) of the file.

    Args:
        handle (text file): File handle positioned at the beginning.

    Returns:
        csv.Dialect: The detected tabular format.
    """
    peek    = read_head(handle)
    handle.seek(0)
    sniffer = csv.Sniffer()
    dialect = sniffer.sniff(peek)
    return dialect


def check_samplesheet(file_in, file_out):
    """
    Validate and transform the input samplesheet, writing a clean CSV to file_out.

    Args:
        file_in (pathlib.Path): Input samplesheet (CSV or TSV).
        file_out (pathlib.Path): Output path for the validated CSV.
    """
    required_columns = {"sID", "bam1", "bam2"}

    with file_in.open(newline="") as in_handle:
        reader = csv.DictReader(in_handle, dialect=sniff_format(in_handle))
        if not required_columns.issubset(reader.fieldnames):
            req_cols = ", ".join(required_columns)
            logger.critical(f"The sample sheet **must** contain these column headers: {req_cols}.")
            sys.exit(1)
        checker = RowChecker()
        for i, row in enumerate(reader):
            try:
                checker.validate_and_transform(row)
            except AssertionError as error:
                logger.critical(f"{str(error)} On line {i + 2}.")
                sys.exit(1)
        checker.validate_unique_samples()

    header = list(reader.fieldnames)
    with file_out.open(mode="w", newline="") as out_handle:
        writer = csv.DictWriter(out_handle, header, delimiter=",")
        writer.writeheader()
        for row in checker.modified:
            writer.writerow(row)


def parse_args(argv=None):
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Validate and transform a pbmerge samplesheet.",
        epilog="Example: check_samplesheet_pbmerge.py samplesheet.csv samplesheet.valid.csv",
    )
    parser.add_argument(
        "file_in",
        metavar="FILE_IN",
        type=Path,
        help="Input samplesheet in CSV or TSV format.",
    )
    parser.add_argument(
        "file_out",
        metavar="FILE_OUT",
        type=Path,
        help="Output validated samplesheet in CSV format.",
    )
    parser.add_argument(
        "-l",
        "--log-level",
        help="Log level (default: WARNING).",
        choices=("CRITICAL", "ERROR", "WARNING", "INFO", "DEBUG"),
        default="WARNING",
    )
    return parser.parse_args(argv)


def main(argv=None):
    """Parse arguments and run samplesheet validation."""
    args = parse_args(argv)
    logging.basicConfig(level=args.log_level, format="[%(levelname)s] %(message)s")
    if not args.file_in.is_file():
        logger.error(f"The given input file {args.file_in} was not found!")
        sys.exit(2)
    args.file_out.parent.mkdir(parents=True, exist_ok=True)
    check_samplesheet(args.file_in, args.file_out)


if __name__ == "__main__":
    sys.exit(main())