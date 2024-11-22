watch:
    watchexec -e py,pyi,toml just test check fmt

test:
    cd lib && uv run pytest

check:
    cd lib && uv run mypy .
    cd lib && uv run ruff check

fmt:
    cd lib && uv run ruff format
