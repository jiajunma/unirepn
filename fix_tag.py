with open('/Users/hoxide/mydoc/unirepn/BMSZ2-revised.tex', 'r') as f:
    content = f.read()

old1 = '$$\nB_T(u,v):=\\langle Tu,v\\rangle,\n\\qquad\nX_{S,T}:=X_0+Y_{C_S,B(T)}.                            \\tag{3}\n$$'
new1 = '\\[\nB_T(u,v):=\\langle Tu,v\\rangle,\n\\qquad\nX_{S,T}:=X_0+Y_{C_S,B(T)}.                            \\tag{3}\n\\]'

if old1 in content:
    content = content.replace(old1, new1)
    print('First replacement done')
else:
    print('First pattern not found')
    # Debug: find the location
    idx = content.find('B_T(u,v)')
    if idx >= 0:
        print(content[idx-5:idx+150])

old2 = '$$\nr=\n\\begin{cases}\nl-1,&D=\\mathbb R,\\ \\epsilon=+1,\\ l\\text{ odd},\\\\\nl,&\\text{otherwise},\n\\end{cases}                                           \\tag{4}\n$$'
new2 = '\\[\nr=\n\\begin{cases}\nl-1,&D=\\mathbb R,\\ \\epsilon=+1,\\ l\\text{ odd},\\\\\nl,&\\text{otherwise},\n\\end{cases}                                           \\tag{4}\n\\]'

if old2 in content:
    content = content.replace(old2, new2)
    print('Second replacement done')
else:
    print('Second pattern not found')
    # Debug: find the location
    idx = content.find('r=\n\\begin{cases}')
    if idx >= 0:
        print(content[idx-10:idx+150])

with open('/Users/hoxide/mydoc/unirepn/BMSZ2-revised.tex', 'w') as f:
    f.write(content)
print('Done')
