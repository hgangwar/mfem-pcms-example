#!/bin/bash
# ============================================================
# format.sh - Run clang-format on all changed C++ files
# ============================================================

module use /opt/scorec/spack/rhel9/v0222_2/lmod/linux-rhel9-x86_64/Core/
module load llvm

CLANG_FORMAT=$(command -v clang-format 2>/dev/null)
CLANG_TIDY=$(command -v clang-tidy 2>/dev/null)

if [ -z "$CLANG_FORMAT" ]; then
    echo "clang-format not found in PATH."
    echo "Install via LLVM tarball or pip in your user space."
else
    echo "Using clang-format at: $CLANG_FORMAT"
fi

if [ -z "$CLANG_TIDY" ]; then
    echo "clang-tidy not found. Skipping static analysis."
else
    echo "Using clang-tidy at: $CLANG_TIDY"
fi

# ------------------------------------------
# Detect modified or staged tracked C++ files
# ------------------------------------------
FILES=$(
  {
    git diff --cached --name-only --diff-filter=ACM
    git diff --name-only --diff-filter=ACM
  } | sort -u | grep -E '\.(cpp|cc|cxx|hpp|hh|hxx|h)$' || true
)

if [ -z "$FILES" ]; then
    echo "No modified C++ files to format or check."
    exit 0
fi

echo "Found the following files:"
echo "$FILES"
echo "------------------------------------------------------------"

# ------------------------------------------
# Step 1: clang-format (style)
# ------------------------------------------
if [ -n "$CLANG_FORMAT" ]; then
    echo "Running clang-format..."
    for f in $FILES; do
        if [ -f "$f" ]; then
            "$CLANG_FORMAT" -i --style=file "$f"
            echo "Formatted: $f"
        fi
    done
    echo "✨ clang-format complete."
fi

# ------------------------------------------
# Step 2: clang-tidy (lint / static analysis)
# ------------------------------------------
if [ -n "$CLANG_TIDY" ]; then
    echo "Running clang-tidy..."
    for f in $FILES; do
        if [ -f "$f" ]; then
            # Run clang-tidy using project settings; silence if no compile_commands.json
            "$CLANG_TIDY" "$f" --quiet --use-color -- -std=c++17 >/dev/null 2>&1
            if [ $? -eq 0 ]; then
                echo "clang-tidy passed: $f"
            else
                echo "clang-tidy warnings in: $f"
            fi
        fi
    done
    echo "clang-tidy analysis complete."
fi

echo "------------------------------------------------------------"
echo "Done. Review with 'git diff' before committing."
exit 0