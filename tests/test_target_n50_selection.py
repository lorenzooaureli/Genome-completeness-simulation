import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from simulate_dna_completeness import select_breakpoints_for_target_n50


def test_returns_input_assembly_when_target_n50_exceeds_total_genome_length():
    state = select_breakpoints_for_target_n50(
        seq="ACGTACGT",
        prioritized_blocks=[],
        block_categories={},
        target_n50=12,
        seed=7,
        mandatory_breakpoints=[],
        min_contig_size=2,
        breakpoint_jitter=0,
    )

    assert state["breakpoints"] == []
    assert state["intervals"] == [(0, 8)]
    assert state["lengths"] == [8]
    assert state["n50"] == 8
    assert state["l50"] == 1
    assert state["contig_count"] == 1
    assert state["events"] == []
    assert state["random_breakpoints"] == 0


def test_returns_input_assembly_when_target_n50_exceeds_input_assembly_n50():
    state = select_breakpoints_for_target_n50(
        seq="ACGTACGTAC",
        prioritized_blocks=[],
        block_categories={},
        target_n50=8,
        seed=7,
        mandatory_breakpoints=[5],
        min_contig_size=2,
        breakpoint_jitter=0,
    )

    assert state["breakpoints"] == [5]
    assert state["intervals"] == [(0, 5), (5, 10)]
    assert state["lengths"] == [5, 5]
    assert state["n50"] == 5
    assert state["l50"] == 1
    assert state["contig_count"] == 2
    assert state["events"] == []
    assert state["random_breakpoints"] == 0


def test_rejects_non_positive_target_n50():
    with pytest.raises(ValueError, match="positive integer"):
        select_breakpoints_for_target_n50(
            seq="ACGTACGT",
            prioritized_blocks=[],
            block_categories={},
            target_n50=0,
            seed=7,
            mandatory_breakpoints=[],
            min_contig_size=2,
            breakpoint_jitter=0,
        )
