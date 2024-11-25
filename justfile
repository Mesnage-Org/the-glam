export PYTHONWARNINGS := "ignore::DeprecationWarning:glycowork"

watch-py:
    watchexec -e py,pyi,toml just test-py check-py fmt-py

test-py:
    cd lib && uv run pytest

check-py:
    cd lib && uv run mypy .
    cd lib && uv run ruff check

# Formats all code in both `lib/` and `web/`
# FIXME: Or will do... Eventually
fmt-py:
    cd lib && uv run ruff format

# Installs the local `glam` package as editable
install-py:
    cd lib && uv pip install -e .

watch-web:
    cd web && bun dev --open
