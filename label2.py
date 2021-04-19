#!/usr/bin/python3
import re

re_comments = re.compile(r"(%.*$)")

token_pattern = r'(?P<label>\\(label)\{([^\}]*)\})'  \
    r'|(?P<ref>\\(ref|eqref|cref|Cref)\{([^\}]*)\})' \
    r'|(?P<cite>\\(cite)(\[[^\]]*\])?\{([^\}]*)\})' \
    r'|(?P<bibitem>\\(bibitem)(\[[^\]]*\])?\{([^\}]*)\})' \
    r'|(?P<bib>\\(bib)\{([^\}]*)\})' \
    r'|(?P<begin>\\(begin)\{(equation|align|eqarray)\})' \
    r'|(?P<end>\\(end)\{(equation|align|eqarray)\})'

tokidx = {'label':2, 'ref':2, 'cite':3, 'bibitem':3, 'bib':2, 'begin':2,
           'end':2}

token_re = re.compile(token_pattern)


def tokenize(f,line_n = 0):
    for line in f:
        line_n += 1
        # filter out comments
        line = re_comments.sub(r"",line)
        while True:
            m = token_re.search(line)
            if not m: break
            line = line[m.end():]
            tokname = m.lastgroup
            # line_number, token_name, label_name
            # print line_n, tokname, m.group(m.lastindex
            #                                +tokidx[tokname])
            kk =  m.group(m.lastindex+tokidx[tokname])
            for k in kk.split(','):
                k=k.strip()
                yield (line_n, tokname, k)


class  ct:
    def __init__(self, f,line_n=0):
        self.tt = tokenize(f)
        self.dd = dict()
        self.leq = dict()
        self.eq = list()
        self.cite = dict()
        self.treat()
    
    def treat(self, ineq=0):
        haslabel = False
        while True:
            dd = self.dd
            eq = self.eq
            leq = self.leq
            cite = self.cite
            try:
                # line_number, token_name, label_name
                ln, tok, lb = next(self.tt)
            except StopIteration:
                if ineq:
                    print ("Warrning unmatched equation at:", ln)
                break
            #print ln, tok, lb
            if tok == 'begin':
                self.treat(ineq=ln)
            elif tok == 'end':
                if ineq : break
                else:
                    print ("Warrning unmatched equation at:", ln)
            elif tok == 'label':
                haslabel = True
                if lb in dd:
                    if dd[lb][1]>=0:
                        dd[lb].append(ln)
                    else:
                        dd[lb][1]=ln
                else:
                    dd[lb] = [0,ln]
                if ineq:
                    if lb in leq:
                        if leq[lb][1]>=0:
                            leq[lb].append(ln)
                        else:
                            leq[lb][1]=ln
                    else:
                        leq[lb] = [0,ln]
                    # print "debug Eq", ineq, lb,ln
            elif tok == 'ref':
                if lb in dd:
                    dd[lb][0] = dd[lb][0]+1
                else:
                    dd[lb] = [1,-ln]
                if lb in leq:
                    leq[lb][0] = leq[lb][0]+1
            elif tok == 'cite':
                if lb in cite:
                    cite[lb][0] = cite[lb][0]+1
                else:
                    cite[lb] = [1,-ln]
            elif tok in {'bibitem','bib'}:
                if lb in cite:
                    cite[lb][1] = ln
                else:
                    cite[lb] = [0,ln]
        if ineq and (not haslabel):
            eq.append(ineq)
        return 


    def summary(self):
        dd, leq, eq, cite = self.dd, self.leq, self.eq, self.cite
        print ("Equation without label")
        for ee in eq:
            print ("Line: ", ee)
        print ("\nref without label")
        for k in sorted(
                filter(lambda k:dd[k][1] <= 0, dd),
                key=(lambda k:-dd[k][1])):
            print ("Line:\t", -dd[k][1], "\t",k)
        print ("\nLabel without reference")        
        for k in sorted(
                filter(lambda k:dd[k][1] > 0 and dd[k][0]==0, dd),
                key=(lambda k:dd[k][1])):
            print ("Line:\t", dd[k][1], "\t",k)

        print ("\nMulti labels")        
        for k in sorted(
                filter(lambda k:len(dd[k]) > 2, dd),
                key=(lambda k:dd[k][1])):
            print ("Line:\t", dd[k][1:], "\t",k)

        print ("\nEquation without reference")        
        for k in sorted(
                filter(lambda k:leq[k][1] > 0 and leq[k][0]==0, leq),
                key=(lambda k:leq[k][1])):
            print ("Line:\t", leq[k][1], "\t",k)
        
        print ("\nBibitem without cite")        
        for k in sorted(
                filter(lambda k:cite[k][1] > 0 and
                       cite[k][0]==0, cite),
                key=(lambda k:cite[k][1])):
            print ("Line:\t", cite[k][1], "\t",k)
        print ("\ncite without Bibitem")        
        for k in sorted(
                filter(lambda k:cite[k][1] <= 0, cite),
                key=(lambda k:-cite[k][1])):
            print ("Line:\t", -cite[k][1], "\t",k)
                
if __name__ == "__main__":
    from sys import argv
    if len(argv) != 2:
        print ("Usage: ", argv[0], "<filename>")
    else:
        fn = argv[1]
        a = ct(open(argv[1],'r'))
        a.summary()
        # for i in tokenize(open(fn,'r')):
        # pass
