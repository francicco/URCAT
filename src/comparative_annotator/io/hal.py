from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from tempfile import NamedTemporaryFile
from typing import Callable, List, Sequence
import subprocess

from comparative_annotator.models.projection import ProjectionInterval


class HALError(RuntimeError):
    pass


@dataclass
class HALCommandResult:
    returncode: int
    stdout: str
    stderr: str


def run_command(cmd: Sequence[str]) -> HALCommandResult:
    """
    Default external command runner.
    Wrapped so it can be mocked in tests.
    """
    proc = subprocess.run(
        list(cmd),
        capture_output=True,
        text=True,
        check=False,
    )
    return HALCommandResult(
        returncode=proc.returncode,
        stdout=proc.stdout,
        stderr=proc.stderr,
    )


class HALAdapter:
    """
    Minimal wrapper around HAL tools.

    Version 1 uses halLiftover on genomic intervals.
    Later versions can add transcript-aware projection and
    exon-chain reconstruction.
    """

    def __init__(
        self,
        hal_path: str | Path,
        runner: Callable[[Sequence[str]], HALCommandResult] = run_command,
    ) -> None:
        self.hal_path = Path(hal_path)
        self.runner = runner

    def project_interval(
        self,
        source_species: str,
        target_species: str,
        seqid: str,
        start: int,
        end: int,
        strand: str,
        source_transcript: str,
    ) -> List[ProjectionInterval]:
        """
        Project one interval from source_species to target_species using halLiftover.

        Input coordinates are assumed to be 1-based inclusive.
        BED uses 0-based half-open, so conversion is needed.
        """
        if start > end:
            raise ValueError(f"Invalid interval: start {start} > end {end}")

        bed_line = self._to_bed_line(
            seqid=seqid,
            start=start,
            end=end,
            name=source_transcript,
            strand=strand,
        )

        with NamedTemporaryFile("w", delete=False) as in_bed:
            in_bed.write(bed_line)
            in_bed_path = in_bed.name

        with NamedTemporaryFile("w", delete=False) as out_bed:
            out_bed_path = out_bed.name

        cmd = [
            "halLiftover",
            str(self.hal_path),
            source_species,
            in_bed_path,
            target_species,
            out_bed_path,
        ]

        result = self.runner(cmd)

        if result.returncode != 0:
            raise HALError(
                f"halLiftover failed for {source_species}->{target_species}: {result.stderr}"
            )

        projections = self._parse_bed_output(
            Path(out_bed_path),
            target_species=target_species,
            source_species=source_species,
            source_transcript=source_transcript,
        )

        return projections

    @staticmethod
    def _to_bed_line(
        seqid: str,
        start: int,
        end: int,
        name: str,
        strand: str,
    ) -> str:
        """
        Convert 1-based inclusive coordinates to BED 0-based half-open.
        BED6: chrom, start, end, name, score, strand
        """
        bed_start = start - 1
        bed_end = end
        return f"{seqid}\t{bed_start}\t{bed_end}\t{name}\t0\t{strand}\n"

    @staticmethod
    def _parse_bed_output(
        bed_path: Path,
        target_species: str,
        source_species: str,
        source_transcript: str,
    ) -> List[ProjectionInterval]:
        """
        Parse BED6-like halLiftover output and convert back to 1-based inclusive.
        """
        intervals: List[ProjectionInterval] = []

        if not bed_path.exists():
            return intervals

        with bed_path.open("r", encoding="utf-8") as fh:
            for raw in fh:
                line = raw.strip()
                if not line:
                    continue

                parts = line.split("\t")
                if len(parts) < 6:
                    continue

                seqid, bed_start, bed_end, _name, _score, strand = parts[:6]

                start = int(bed_start) + 1
                end = int(bed_end)

                intervals.append(
                    ProjectionInterval(
                        species=target_species,
                        seqid=seqid,
                        start=start,
                        end=end,
                        strand=strand,
                        source_species=source_species,
                        source_transcript=source_transcript,
                    )
                )

        return intervals
