from typing import IO
from collections.abc import Iterator

class FastaEntry:
    description: str
    sequence: str

def read(source: str | IO | None) -> Iterator[FastaEntry]: ...
