#!/usr/bin/env python3
"""Check that relative Markdown links in the repository's documentation resolve.

Usage (from anywhere):  python3 docs/check_links.py [--all]

By default checks the entry-point documents (CLAUDE.md, README.md, docs/*.md, notebooks/README.md,
python/**/NOTES.md, verification/**/README.md, data/**/README.md). ``--all`` checks every ``*.md`` file
outside ``.venv`` and ``legacy/``. Only relative links to files or directories are checked; ``http(s)``,
``mailto`` and pure ``#anchor`` links are skipped, as are links inside fenced code blocks and inline code.
Exit status 1 if any link is broken.
"""
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
LINK = re.compile(r"(?<!\!)\[[^\]\n]*\]\(([^)\s]+)(?:\s+\"[^\"]*\")?\)")
INLINE_CODE = re.compile(r"`[^`\n]*`")


def markdown_files(check_all):
    skip = {".venv", "legacy", ".git", "node_modules", "__pycache__"}
    for path in sorted(ROOT.rglob("*.md")):
        rel = path.relative_to(ROOT)
        if skip & set(rel.parts):
            continue
        if check_all or rel.parts[0] in {"docs", "notebooks", "python", "verification", "data"} or len(rel.parts) == 1:
            yield path


def broken_links(path):
    in_fence = False
    for number, line in enumerate(path.read_text(encoding="utf-8", errors="replace").splitlines(), 1):
        if line.lstrip().startswith("```"):
            in_fence = not in_fence
            continue
        if in_fence:
            continue
        for match in LINK.finditer(INLINE_CODE.sub("", line)):
            target = match.group(1)
            if re.match(r"^(https?:|mailto:|#)", target):
                continue
            target = target.split("#", 1)[0]
            if target and not (path.parent / target).exists():
                yield number, target


def main():
    check_all = "--all" in sys.argv[1:]
    total = bad = 0
    for path in markdown_files(check_all):
        total += 1
        for number, target in broken_links(path):
            bad += 1
            print(f"BROKEN  {path.relative_to(ROOT)}:{number}  ->  {target}")
    print(f"checked {total} Markdown files, {bad} broken relative link(s)")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
