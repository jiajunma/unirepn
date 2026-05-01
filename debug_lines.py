#!/usr/bin/env python3

# Read the file
with open('BMSZ2-revised.tex', 'r') as f:
    lines = f.readlines()

# Print lines around 7705-7717 (0-indexed: 7704-7716)
print("Lines 7705-7718 (1-indexed):")
for i in range(7704, min(7718, len(lines))):
    print(f"{i+1}: {repr(lines[i])}")

# Now let's find "Finally," in the file
for i, line in enumerate(lines):
    if 'Finally,' in line and i > 7000:  # Find it around our area
        print(f"\nFound 'Finally,' at line {i+1}: {repr(line)}")
        # Check what's before and after
        if i+1 < len(lines):
            print(f"Line {i+2}: {repr(lines[i+1])}")
        break