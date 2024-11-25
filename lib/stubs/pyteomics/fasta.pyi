from collections import namedtuple
from collections.abc import Iterator
from typing import IO

Protein = namedtuple("Protein", ("description", "sequence"))

def read(source: str | IO | None, use_index: bool = False) -> Iterator[Protein]: ...
