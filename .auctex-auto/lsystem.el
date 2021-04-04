(TeX-add-style-hook
 "lsystem"
 (lambda ()
   (setq TeX-command-extra-options
         "-shell-escape")
   (TeX-add-to-alist 'LaTeX-provided-class-options
                     '(("subfiles" "ssunip")))
   (add-to-list 'LaTeX-verbatim-environments-local "lstlisting")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "path")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "url")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "nolinkurl")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "hyperbaseurl")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "hyperimage")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "hyperref")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "href")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "lstinline")
   (add-to-list 'LaTeX-verbatim-macros-with-delims-local "path")
   (add-to-list 'LaTeX-verbatim-macros-with-delims-local "lstinline")
   (TeX-run-style-hooks
    "latex2e"
    "subfiles"
    "subfiles10")
   (TeX-add-symbols
    "oAC"
    "pAC"
    "nAC")
   (LaTeX-add-labels
    "tail0"
    "tailtip"
    "prop:CC.bij"
    "eq:DD.CC"
    "eq:DD.CC1"
    "prop:D.sign"
    "prop:delta"
    "eq:delta"
    "eq:sign.D"
    "eq:delta.I"
    "eq:sign.GD"
    "eq:delta.1"
    "eq:ydelta"
    "eq:sign.GD2"
    "cor:D.inj1"
    "eq:D.BD"
    "cor:dpinj"
    "prop:LLS"
    "p:drcls.1"
    "eq:LS.dis"
    "it:LLS.4"
    "eq:up"
    "sec:pfDC.init"
    "c:init.CD"
    "sec:pf.ds.CD"
    "eq:LS.CD.inj"
    "eq:eps.CD1"
    "eq:LS.D.ds"
    "sec:pf.gd.CD"
    "eq:LS.taup"
    "eq:lsign.1"
    "eq:gd.ls"
    "sec:z.r"
    "eq:rr.c"
    "eq:gd.rr"
    "eq:rd.c"
    "sec:z.c"
    "eq:ss.c"
    "eq:gd.ss"
    "sec:z.d"
    "eq:ped.ssc"
    "eq:ped.ssd"
    "c:gd.C1"
    "eq:uptaupp.sign"
    "c:gd.C3"
    "c:gd.D1"
    "c:gd.C2"
    "c:gd.C2.1"
    "c:noticed.bij"
    "it:c:noticed.bij.4"
    "c:d+"
    "c:gd.pnoticed"
    "c:gd.pnoticed.p"
    "eq:pnoticed.1"
    "c:gd.noticed.inj"
    "c:gd.pnoticed.n"))
 :latex)

