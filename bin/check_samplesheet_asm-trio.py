#!/usr/bin/env python3
"""
Validate and transform the input samplesheet for the lrs_asm_trio pipeline.

Checks that the samplesheet contains the required columns (sID,
proband_hifi_fastq, mat_hifi_fastq, pat_hifi_fastq, HPO), validates file
extensions, and ensures sample IDs are unique. Writes a validated CSV to
file_out.

Usage:
    check_samplesheet_asm-trio.py <samplesheet.csv> <samplesheet.valid.csv>

Expected samplesheet format:
    sID,proband_hifi_fastq,mat_hifi_fastq,pat_hifi_fastq,HPO
    SAMPLE01,/path/to/proband.fastq.gz,/path/to/mat.fastq.gz,/path/to/pat.fastq.gz,/path/to/HPO.txt
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
    Validate and transform each row of the trio samplesheet.

    Attributes:
        modified (list): Validated and transformed rows in original order.

    Notes:
        If the same sample ID appears more than once, sample IDs are renamed
        with a _T{n} suffix (e.g. SAMPLE01_T1, SAMPLE01_T2). This affects
        output filenames downstream — ensure this behaviour is expected before
        running with multi-sample samplesheets containing duplicate IDs.
    """

    VALID_FASTQ_FORMATS = (".fastq", ".fastq.gz", ".fq", ".fq.gz")
    VALID_HPO_FORMATS   = (".txt",)

    def __init__(
        self,
        sample_col="sID",
        proband_col="proband_hifi_fastq",
        mat_col="mat_hifi_fastq",
        pat_col="pat_hifi_fastq",
        hpo_col="HPO",
        **kwargs
    ):
        """
        Initialize the row checker with expected column names.

        Args:
            sample_col (str): Column containing the sample identifier (default: sID).
            proband_col (str): Column containing the proband HiFi FASTQ path.
            mat_col (str): Column containing the maternal HiFi FASTQ path.
            pat_col (str): Column containing the paternal HiFi FASTQ path.
            hpo_col (str): Column containing the HPO terms file path.
        """
        super().__init__(**kwargs)
        self._sample_col  = sample_col
        self._proband_col = proband_col
        self._mat_col     = mat_col
        self._pat_col     = pat_col
        self._hpo_col     = hpo_col
        self._seen        = set()
        self.modified     = []

    def validate_and_transform(self, row):
        """
        Run all validations on a row and append to modified list if valid.

        Args:
            row (dict): Mapping from column headers to row values.
        """
        self._validate_sample(row)
        self._validate_fastq(row, self._proband_col, "proband")
        self._validate_fastq(row, self._mat_col, "maternal")
        self._validate_fastq(row, self._pat_col, "paternal")
        self._validate_hpo(row)
        self._seen.add(row[self._sample_col])
        self.modified.append(row)

    def _validate_sample(self, row):
        """Assert that the sample ID is non-empty and replace spaces with underscores."""
        if len(row[self._sample_col]) <= 0:
            raise AssertionError("Sample ID (sID) is required.")
        row[self._sample_col] = row[self._sample_col].replace(" ", "_")

    def _validate_fastq(self, row, col, label):
        """Assert that a FASTQ path is non-empty and has a valid extension."""
        if len(row[col]) <= 0:
            raise AssertionError(f"The {label} HiFi FASTQ file path is required.")
        if not any(row[col].endswith(ext) for ext in self.VALID_FASTQ_FORMATS):
            raise AssertionError(
                f"Unrecognized {label} FASTQ file extension: {row[col]}\n"
                f"Expected one of: {', '.join(self.VALID_FASTQ_FORMATS)}"
            )

    def _validate_hpo(self, row):
        """Assert that the HPO file path has a valid extension if provided."""
        if len(row[self._hpo_col]) > 0:
            if not any(row[self._hpo_col].endswith(ext) for ext in self.VALID_HPO_FORMATS):
                raise AssertionError(
                    f"Unrecognized HPO file extension: {row[self._hpo_col]}\n"
                    f"Expected one of: {', '.join(self.VALID_HPO_FORMATS)}"
                )

    def validate_unique_samples(self):
        """
        Assert that each sample ID is unique.

        If the same sample ID appears more than once, all instances are renamed
        with a _T{n} suffix to distinguish them.
        """
        if len(self._seen) != len(self.modified):
            raise AssertionError("Each sample ID must be unique.")
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
    required_columns = {"sID", "proband_hifi_fastq", "mat_hifi_fastq", "pat_hifi_fastq", "HPO"}

    with file_in.open(newline="") as in_handle:
        reader = csv.DictReader(in_handle, dialect=sniff_format(in_handle))
        if not required_columns.issubset(reader.fieldnames):
            req_cols = ", ".join(required_columns)
            logger.critical(f"Samplesheet must contain these column headers: {req_cols}.")
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
        description="Validate and transform a trio samplesheet.",
        epilog="Example: check_samplesheet_asm-trio.py samplesheet.csv samplesheet.valid.csv",
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
        logger.error(f"Input file not found: {args.file_in}")
        sys.exit(2)
    args.file_out.parent.mkdir(parents=True, exist_ok=True)
    check_samplesheet(args.file_in, args.file_out)


if __name__ == "__main__":
    sys.exit(main())