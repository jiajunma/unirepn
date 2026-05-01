#!/usr/bin/env python3

# Read the file
with open('BMSZ2-revised.tex', 'r') as f:
    content = f.read()

lines = content.split('\n')

# Find the "Finally," that's at line 7705 (1-indexed)
# In 0-indexed, that's line 7704
target_line = 7704  # 0-indexed

# Build the old section from lines 7704 onwards
# We want lines 7705-7717 (1-indexed), which is 7704-7716 (0-indexed)
old_section_lines = lines[7704:7716+1]  # +1 because slice end is exclusive
print(f"Old section has {len(old_section_lines)} lines")
print(f"First line: {repr(old_section_lines[0])}")
print(f"Last line: {repr(old_section_lines[-1])}")

old_section = '\n'.join(old_section_lines) + '\n'
print(f"\nOld section length: {len(old_section)}")

# Check if this is in the content
if old_section in content:
    print("Found old section in content!")
else:
    print("Old section NOT found in content")
    # Try to find where "Finally," appears after position 7700
    pos = content.find("Finally,")
    while pos != -1:
        line_num = content[:pos].count('\n') + 1
        if line_num > 7700 and line_num < 7800:
            print(f"Found 'Finally,' at byte position {pos}, line {line_num}")
            print(f"Context: {repr(content[pos-50:pos+100])}")
            break
        pos = content.find("Finally,", pos + 1)