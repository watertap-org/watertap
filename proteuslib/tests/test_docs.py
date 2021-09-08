"""
Tests for documentation
"""
from pathlib import Path
import re

__authors__ = ["Dan Gunter"]


def get_autodoc(pth):
    """Get the name, location, and options for all autodoc directives
    in the file at 'pth'.

    Yields:
        tuple of length 4: type of thing, name of it, filename, line number (from 1)
    """
    mode = "text"
    auto_what, auto_name, auto_options, auto_where = None, None, None, None
    with open(pth, encoding="utf-8") as f:
        for i, line in enumerate(f):
            if mode == "text":
                m = re.match(r'..\s+auto(\w+)::\s*([a-zA-Z][a-zaA-Z_.0-9]*)', line)
                if m:
                    auto_what = m.group(1)
                    auto_name = m.group(2)
                    auto_where = f"{f.name}:{i + 1}"
                    auto_options = []
                    mode = "autodoc"
            elif mode == "autodoc":
                # blank or unindented line, back to text mode
                if not re.match(r'^\s\s', line):
                    yield auto_what, auto_name, auto_where, auto_options
                    auto_what = None
                    mode = "text"
                else:
                    m = re.match(r"\s*:(\w+):", line)
                    if m:
                        option = m.group(1)
                        auto_options.append(option)
            else:
                raise RuntimeError(f"unknown mode: {mode}")
        if auto_what:
            yield auto_what, auto_name, auto_where, auto_options


def test_autodoc_has_noindex():
    docs_dir = Path(__file__).parent.parent.parent / "docs"
    # this path may only work in developer install
    if docs_dir.exists() and docs_dir.is_dir():
        root = docs_dir / "technical_reference"
        bad = []
        for p in root.glob("**/*.rst"):
            for what, name, where, options in get_autodoc(p):
                if "noindex" not in options:
                    bad.append(where)
        nl = "\n"
        if len(bad) > 0:
            err = f"{len(bad)} 'autodoc' directive(s) missing 'noindex':\n" \
                  f"{nl.join(('%d) %s' % (i+1, n) for i, n in enumerate(bad)))}"
            assert False, err
