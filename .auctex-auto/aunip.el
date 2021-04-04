(TeX-add-style-hook
 "aunip"
 (lambda ()
   (setq TeX-command-extra-options
         "-shell-escape")
   (TeX-add-to-alist 'LaTeX-provided-class-options
                     '(("amsart" "12pt" "a4paper")))
   (TeX-add-to-alist 'LaTeX-provided-package-options
                     '(("geometry" "margin=2.5cm" "marginpar=2cm") ("hyperref" "bookmarksopen" "bookmarksdepth=3" "hidelinks") ("cleveref" "nameinlink") ("showkeys" "color") ("xy" "all" "cmtip") ("xcolor" "rgb" "table" "dvipsnames") ("ulem" "normalem") ("ytableau" "centertableaux") ("standalone" "mode=buildnew")))
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
    "amsart"
    "amsart12"
    "geometry"
    "hyperref"
    "cleveref"
    "showkeys"
    "array"
    "amssymb"
    "amsmath"
    "mathrsfs"
    "mathbbol"
    "mathabx"
    "amsthm"
    "graphicx"
    "braket"
    "mathtools"
    "amsrefs"
    "xy"
    "rotating"
    "leftidx"
    "pifont"
    "xcolor"
    "imakeidx"
    "ulem"
    "ytableau"
    "enumitem"
    "xparse"
    "diagbox"
    "arydshln"
    "standalone"
    "tikz"
    "etoolbox"
    "upgreek"
    "listings"
    "subfiles")
   (TeX-add-symbols
    '("byhide" ["argument"] 1)
    '("trivial" ["argument"] 1)
    '("tytb" 1)
    '("cover" 1)
    '("sfrac" 2)
    '("intn" 1)
    '("wpair" 1)
    '("pair" 1)
    '("mjjc" 1)
    '("lokec" 1)
    '("circnuma" 1)
    "Rank"
    "cqq"
    "rsym"
    "rskew"
    "fraksp"
    "frakso"
    "frakm"
    "frakp"
    "pr"
    "rhopst"
    "Rad"
    "Res"
    "Hol"
    "WF"
    "AV"
    "AVC"
    "VC"
    "bfv"
    "depth"
    "wtM"
    "wtMone"
    "nullpp"
    "nullp"
    "bfonenp"
    "bfonepn"
    "bfone"
    "piSigma"
    "piSigmap"
    "sfVprime"
    "sfVdprime"
    "gminusone"
    "eva"
    "bcN"
    "iso"
    "riso"
    "BA"
    "BC"
    "BD"
    "BE"
    "BF"
    "BG"
    "BH"
    "BI"
    "BJ"
    "BK"
    "BL"
    "BM"
    "BN"
    "BO"
    "BP"
    "BQ"
    "BR"
    "BS"
    "BT"
    "BU"
    "BV"
    "BW"
    "BX"
    "BY"
    "BZ"
    "Bk"
    "CA"
    "CB"
    "CC"
    "CE"
    "CF"
    "CG"
    "CH"
    "CI"
    "CJ"
    "CK"
    "CL"
    "CM"
    "CN"
    "CO"
    "CP"
    "CQ"
    "CR"
    "CS"
    "CT"
    "CU"
    "CV"
    "CW"
    "CX"
    "CY"
    "CZ"
    "RA"
    "RB"
    "RC"
    "RD"
    "RE"
    "RF"
    "RG"
    "RH"
    "RI"
    "RJ"
    "RK"
    "RL"
    "RM"
    "RN"
    "RO"
    "RP"
    "RQ"
    "RS"
    "RT"
    "RU"
    "RV"
    "RW"
    "RX"
    "RY"
    "RZ"
    "cod"
    "cont"
    "cl"
    "cusp"
    "disc"
    "Gm"
    "I"
    "Jac"
    "PM"
    "new"
    "NS"
    "N"
    "ord"
    "rk"
    "rr"
    "rh"
    "Sel"
    "Sim"
    "wt"
    "wh"
    "ds"
    "ov"
    "incl"
    "lra"
    "imp"
    "bs"
    "Norma"
    "Ima"
    "con"
    "gr"
    "ad"
    "der"
    "dif"
    "pro"
    "Ev"
    "Invf"
    "Inv"
    "HC"
    "lef"
    "righ"
    "Diff"
    "diag"
    "sh"
    "sch"
    "open"
    "sgn"
    "triv"
    "Sh"
    "oN"
    "oc"
    "od"
    "os"
    "ol"
    "oL"
    "oJ"
    "oH"
    "oO"
    "oS"
    "oR"
    "oT"
    "oZ"
    "oD"
    "oW"
    "oE"
    "oP"
    "PD"
    "oU"
    "gC"
    "gl"
    "re"
    "g"
    "h"
    "p"
    "Z"
    "R"
    "Q"
    "M"
    "B"
    "V"
    "W"
    "F"
    "E"
    "Y"
    "FF"
    "HH"
    "ve"
    "aut"
    "ii"
    "jj"
    "kk"
    "la"
    "ra"
    "bp"
    "be"
    "ee"
    "LRleq"
    "noticed"
    "ess"
    "dotminus"
    "idxemph"
    "okay"
    "editc"
    "mjj"
    "mjjr"
    "mjjd"
    "mjjb"
    "mjje"
    "mjjcb"
    "mjjce"
    "sun"
    "sund"
    "mv"
    "delete"
    "Ueven"
    "Uodd"
    "ttau"
    "Wcp"
    "Kur"
    "Im"
    "gen"
    "inn"
    "ta"
    "tb"
    "binn"
    "innwi"
    "innw"
    "innv"
    "innbfv"
    "innvi"
    "innvp"
    "innp"
    "simrightarrow"
    "surj"
    "usecsname"
    "useLetter"
    "usedbletter"
    "mydefcirc"
    "mydefvec"
    "mydefdot"
    "mydefacute"
    "mydefbr"
    "mydefbar"
    "mydefhat"
    "mydefwh"
    "mydeft"
    "mydefu"
    "mydefr"
    "mydefb"
    "mydefwt"
    "mydefbf"
    "mydefc"
    "mydefsf"
    "mydefs"
    "mydefcks"
    "mydefckc"
    "mydefck"
    "abs"
    "norm"
    "fsl"
    "fsp"
    "YD"
    "SYD"
    "MK"
    "MYD"
    "AND"
    "deti"
    "AOD"
    "oAC"
    "owAC"
    "pac"
    "nac"
    "ttail"
    "Forall"
    "AC"
    "wAC"
    "ac"
    "lotimes"
    "KM"
    "ckbfG"
    "eDD"
    "eDDo"
    "DD"
    "DDc"
    "gDD"
    "gDDc"
    "flushl"
    "flushr"
    "flushmr"
    "cpc"
    "ccJ"
    "ccL"
    "wtbfK"
    "AbfV"
    "abfV"
    "afgg"
    "abfG"
    "half"
    "ihalf"
    "slt"
    "sltr"
    "slee"
    "slff"
    "slhh"
    "sleei"
    "slxx"
    "slyy"
    "slxxi"
    "slH"
    "Mop"
    "fggJ"
    "fggJp"
    "NilGC"
    "NilGCp"
    "Nilgp"
    "Nilg"
    "peNil"
    "dpeNil"
    "nNil"
    "eNil"
    "KS"
    "MM"
    "MMP"
    "IST"
    "tIST"
    "gpi"
    "bfWo"
    "bfWoo"
    "bfWg"
    "Xg"
    "Xo"
    "Xoo"
    "fppo"
    "fggo"
    "dliftv"
    "slift"
    "bbThetav"
    "tsign"
    "lsign"
    "bsign"
    "ssign"
    "dsign"
    "bcO"
    "scratch"
    "defpcmd"
    "KK"
    "A"
    "K"
    "G"
    "J"
    "L"
    "eps"
    "pp"
    "fggR"
    "rmtop"
    "dimo"
    "JW"
    "floor"
    "KSP"
    "UUC"
    "tUUC"
    "OmegabfW"
    "BB"
    "X"
    "Lslt"
    "Xslt"
    "eslt"
    "thetaO"
    "Thetav"
    "thetav"
    "Thetab"
    "cKaod"
    "mstar"
    "GVr"
    "tGVr"
    "GVpr"
    "tGVpr"
    "GVar"
    "tGVar"
    "GV"
    "GVp"
    "KVr"
    "tKVr"
    "KV"
    "KaV"
    "acO"
    "asO"
    "sp"
    "bfLz"
    "sOpe"
    "sOpeR"
    "sOR"
    "gdliftv"
    "gdlift"
    "bcOp"
    "bsO"
    "bsOp"
    "bfVpe"
    "bfEz"
    "bfVn"
    "bfEzp"
    "totimes"
    "dotbfV"
    "aod"
    "unip"
    "ssP"
    "ssD"
    "ssdd"
    "phik"
    "phikp"
    "bbfK"
    "brrho"
    "whAX"
    "mktvvp"
    "Piunip"
    "cf"
    "Groth"
    "Irr"
    "edrc"
    "drc"
    "drcs"
    "drcns"
    "LS"
    "LLS"
    "LSaod"
    "Unip"
    "lUnip"
    "tbfxx"
    "PBPe"
    "PBPes"
    "PBPesp"
    "pbp"
    "pbpst"
    "pbpssp"
    "pbpsns"
    "pbpsp"
    "pbpns"
    "DDn"
    "dsrcd"
    "taupna"
    "tauna"
    "enon"
    "ytb"
    "ckcOp"
    "ckcOpp"
    "cOp"
    "cOpp"
    "cLpp"
    "cLppp"
    "pUpsilon"
    "nUpsilon"
    "pcL"
    "ncL"
    "pcP"
    "ncP"
    "pcT"
    "ncT"
    "pcC"
    "ncC"
    "pcLp"
    "ncLp"
    "pcLpp"
    "ncLpp"
    "pcB"
    "ncB"
    "uptaup"
    "uptaupp"
    "uptauppp"
    "bdelta"
    "tcO"
    "tcOp"
    "tcOpp"
    "tuptau"
    "tuptaup"
    "tuptaupp"
    "tuptauppp"
    "taup"
    "taupp"
    "tauppp"
    "cpT"
    "cnT"
    "cpB"
    "cnB"
    "BOX"
    "ckDD"
    "deltas"
    "deltans"
    "PP"
    "uum"
    "uup"
    "LEG"
    "PBP"
    "BODY"
    "eee"
    "upp"
    "umm"
    "bipartl"
    "bipartr"
    "dsdiagl"
    "dsdiagr"
    "DDl"
    "DDr"
    "GLEz"
    "GLE"
    "wtGLE"
    "wtGLEz"
    "wtPE"
    "JU"
    "LU"
    "wtGU"
    "dsfss"
    "UU"
    "fggs"
    "fggsp"
    "fggspo"
    "fggspt"
    "fggspi"
    "fkks"
    "fkksp"
    "fkkspo"
    "fkkspt"
    "fkkspi"
    "fpps"
    "fppsp"
    "fppspo"
    "fppspt"
    "fppspi"
    "DDss"
    "DDsso"
    "Mss"
    "Ms"
    "Msp"
    "Gs"
    "Gsp"
    "Gspo"
    "Gspt"
    "CMs"
    "CMsp"
    "CMss"
    "Wss"
    "Woss"
    "CXss"
    "Xp"
    "Xpo"
    "ww"
    "wwo"
    "Vs"
    "Vsp"
    "Vspo"
    "Vspt"
    "Vspi"
    "Wssi"
    "Wsso"
    "Wsst"
    "mX"
    "kX"
    "ZX"
    "ZXO"
    "sqii"
    "St"
    "VV"
    "SLT"
    "SLTK"
    "GC"
    "CCSS"
    "DCO"
    "cEp"
    "acm"
    "acme"
    "opac"
    "onac"
    "tpac"
    "tnac"
    "dbM"
    "dbMM"
    "dbX"
    "dbfpp"
    "ZdbX"
    "aV"
    "ssfkkp"
    "sscX"
    "oocX"
    "zzcX"
    "ozcX"
    "zocX"
    "iicX"
    "Vker"
    "XET"
    "XST")
   (LaTeX-add-labels
    "sec:intro"
    "secsu"
    "typebcd"
    "ugz"
    "chico"
    "secbip"
    "eq:BOX"
    "eq:sp-nsp.C"
    "defpbp0"
    "eqbp"
    "thmcount"
    "subsec:comTOrep"
    "eq:DD.wp"
    "eq:def-pi"
    "thm1"
    "cor1"
    "defaod"
    "kgroupaod"
    "thmac0"
    "sec:comb"
    "lstarco"
    "eq:def.alphap"
    "lemDDn1"
    "lemDDn2"
    "sec:desc"
    "descb"
    "descb2"
    "descd1"
    "descd2"
    "prop:DD.BDinj"
    "sec:Nil"
    "lem:cartan"
    "comdet"
    "secmmap"
    "momentmap"
    "descko"
    "kkpo"
    "momentmap2"
    "kkpo2"
    "subsecass"
    "dualdesc"
    "thmac1"
    "thmac2"
    "thmac3"
    "thmac4"
    "thmac5"
    "thmpitau"
    "thmac7"
    "bijthm1"
    "bijthm2"
    "sec:Integrals"
    "secoscil"
    "vn"
    "deforos"
    "iota0"
    "isoggp"
    "identifyjg"
    "boundpsi"
    "int"
    "nuchi"
    "growthdp"
    "psi123"
    "matrico"
    "defn:CRcov"
    "convint00"
    "lemconv"
    "convint01"
    "thetab0"
    "boundm"
    "weaklycont"
    "defn:CR33"
    "positivity000"
    "positivity"
    "intxi2"
    "iota1"
    "isosfss1"
    "boundxx"
    "convint0004"
    "intpi0004"
    "convint0011"
    "doublelift"
    "thetabv00"
    "intt00"
    "sec:DP"
    "chi0"
    "chip"
    "midp"
    "convint00112"
    "opencell"
    "modulus"
    "convint001123"
    "isoipi"
    "imb"
    "intguv0"
    "imb2"
    "chid"
    "degens"
    "lem:coinv"
    "doublelift5"
    "bound"
    "sec:AC"
    "liftop"
    "lem:GDS.AC"
    "prop:GDS.AC"
    "lem:RDS.C"
    "lem:red"
    "eq:tt1"
    "eq:GG"
    "eq:GG2"
    "lem:GDS.sh"
    "it:GDS.sh.2"
    "it:GDS.sh.1"
    "it:GDS.sh.3"
    "eq:GDS.T2"
    "eq:GDS.T3"
    "eq:GDS.T5"
    "sec:proof"
    "eq:Vl.1"
    "eq:ssrel"
    "lem:KX1"
    "eg:MYD"
    "eq:LCD"
    "eq:LDC"
    "sec:tail"
    "tail0"
    "tailtip"
    "prop:CC.bij"
    "eq:DD.CC"
    "eq:DD.CC1"
    "lem:D.sign"
    "lem:delta"
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
    "sec:ACC"
    "lem:dlift"
    "eq:C"
    "eq:BD"
    "eq:TBD"
    "eq:BD2"
    "sec:ac"
    "lem:ac0"
    "lem:C"
    "lem:BD"
    "lem:C*"
    "lem:BD2"
    "lem:BD3"
    "eq:up"
    "sec:pfDC.init"
    "eq:dsign1"
    "convint004"
    "intpi004"
    "eq:UEE22"
    "embvsf"
    "eq:Vperp.dec"
    "est002"
    "bnuv"
    "defqi"
    "sec:DPandRQ"
    "dimu"
    "sec:def.AC"
    "sec:ATL"
    "cor:Cbound"
    "sec:EAC"
    "it:case.D"
    "it:case.GD"
    "defn:tlift.rho"
    "defn:DS.ch"
    "sec:dpmm"
    "clinv"
    "sec:lift.AC"
    "eq:DS.chc"
    "sec:rG"
    "def:J"
    "tab:realforms"
    "def:L"
    "defn:eespace"
    "def:Vsign"
    "tab:sign"
    "sec:MC"
    "eq:P_eps"
    "def:NilC"
    "eq:E.depth"
    "eq:dd"
    "defdo"
    "def:c"
    "subsec:SYD"
    "eq:k+p"
    "def:dec.rP"
    "dedd"
    "def:dec.sNG"
    "sec:KX"
    "lem:char.res"
    "isoo"
    "sec:descent"
    "def:DP"
    "sec:LD"
    "eq:def.LsO22"
    "def:LC"
    "def:GD"
    "eq:GD"
    "gendec"
    "eq:GD.min"
    "def:GD.good"
    "lem:GDS.set"
    "sec:alpha"
    "eq:alpha"
    "eq:alpha_l"
    "lem:alpha.e"
    "lem:char.surj"
    "sec:LVB"
    "idenkr"
    "eq:dec.KO"
    "eq:def.alpha1"
    "labst"
    "defn:glift.rho"
    "sec:aod"
    "def:admD"
    "idenip"
    "lem:Kaod"
    "eq:l.adm"
    "lem:admchar.surj"
    "eq:rho0.sp"
    "intgrability"
    "sec:MCI"
    "estosc"
    "estosc2"
    "defn:CR"
    "intpios"
    "intpi"
    "intpi2"
    "positivity0"
    "intpo"
    "prop:calas"
    "prop:Ch.eq"
    "eq:LCh"
    "boundch"
    "thetabv"
    "intt"
    "subsec:induced"
    "thm:Bar"
    "eq:dim-ine"
    "lem:indC"
    "it:indC.1"
    "it:indC.2"
    "lem:indR"
    "it:indR.1"
    "it:indR.2"
    "sec:PC.ro"
    "eq:decomind1"
    "chop1"
    "chop2"
    "sec:PC.rsp"
    "decomind2"
    "chops1"
    "eq:Sp.gdd"
    "chops2"
    "chopqs0"
    "chopqs1"
    "chopqs2"
    "decomindq"
    "chopq1"
    "chopq2"
    "sec:unipot"
    "sec:cons"
    "eq:eta"
    "intunip"
    "thmunip"
    "eqchj"
    "eq:pij"
    "anatheta"
    "algtheta"
    "sec:GM"
    "sec:F.M"
    "lem:F.cl"
    "sec:Sdes"
    "lem:DS.sh"
    "it:DS.G1"
    "it:DS.G3"
    "it:DS.G2"
    "eq:DS.T1"
    "lem:DS.U"
    "it:DS.U1"
    "it:DS.U2"
    "eq:ES.UD"
    "lem:GDS.U"
    "it:GDS.U1"
    "it:GDS.U2"
    "sec:pf.indR"
    "tstar"
    "eq:stab.X1"
    "eq:indorb")
   (LaTeX-add-index-entries
    "#1")
   (LaTeX-add-amsthm-newtheorems
    "thm"
    "thml"
    "lem"
    "obs"
    "lemt"
    "whyp"
    "prop"
    "prpt"
    "prpl"
    "cor"
    "claim"
    "defn"
    "dfnl"
    "IndH"
    "eg"
    "remark"
    "remarks"
    "Example")
   (LaTeX-add-xcolor-definecolors
    "srcol")
   (LaTeX-add-enumitem-newlists
    '("enumC" "enumerate")
    '("enumT" "enumerate")
    '("enumPF" "enumerate")
    '("enumS" "enumerate")
    '("enumI" "enumerate")
    '("enumIL" "enumerate*")
    '("enumR" "enumerate")
    '("des" "enumerate")
    '("enuma" "enumerate"))
   (LaTeX-add-xparse-macros
    '("cent" "o m ")
    '("C" "")
    '("NilP" "t'")
    '("KTW" "o g")
    '("CHI" "o g")
    '("PR" "g")
    '("XX" "g")
    '("PP" "g")
    '("LL" "g")
    '("ZZ" "g")
    '("WW" "g")
    '("KK" "g")
    '("XXo" "d()")
    '("ZZo" "g")
    '("bcO" "t'")
    '("oliftc" "g")
    '("oliftr" "g")
    '("olift" "g")
    '("tlift" "g")
    '("NN" "g")
    '("RR" "m m")
    '("sign" "m")
    '("lnn" "t+ t- g")
    '("KC" "s o e{_}")
    '("LW" "g")
    '("dlift" "O{\\sfss'} O{\\sfss}")))
 :latex)

