from typing import Pattern

def cleave(
    sequence: str,
    rule: str | Pattern[str],
    missed_cleavages: int,
    min_length: int | None,
    max_length: int | None,
    semi: bool,
) -> set[str]: ...
