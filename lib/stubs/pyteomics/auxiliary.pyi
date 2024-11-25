class PyteomicsError(Exception):
    message: str

    def __init__(self, msg: str) -> None: ...
