#!/usr/bin/env python3

with open('BMSZ2-revised.tex', 'r') as f:
    content = f.read()

# Find and replace the duplicate subsection
old = """\\end proof}
\\subsection{The generalized good descent case}\\label{sec:app-gd}

\\subsection{The generalized good descent case}\\label{sec:app-gd}"""

new = """\\end proof}
\\subsection{The generalized good descent case}\\label{sec:app-gd}"""

if old in content:
    content = content.replace(old, new)
    with open('BMSZ2-revised.tex', 'w') as f:
        f.write(content)
    print("Replacement successful!")
else:
    print("Old pattern not found, trying alternate pattern")
    # Try with subsection (subsec) instead
    old2 = """\\end proof}
\\subsection{The generalized good descent case}\\label{sec:app-gd}

\\subsection{The generalized good descent case}\\label{sec:app-gd}"""
    new2 = """\\end proof}
\\subsection{The generalized good descent case}\\label{sec:app-gd}"""
    if old2 in content:
        content = content.replace(old2, new2)
        with open('BMSZ2-revised.tex', 'w') as f:
            f.write(content)
        print("Replacement successful with alternate!")
    else:
        print("Alternate pattern not found either")
        # Let's find what's actually there
        idx = content.find("\\end proof}")
        if idx != -1:
            print("Context after \\end proof}:")
            print(repr(content[idx:idx+200]))