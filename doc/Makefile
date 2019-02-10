
doc = slspack.pdf

SRC_TEX   = $(wildcard *.tex)
TEXFILES = $(SRC_TEX)

LATEX = latex
DVIPDF = dvipdf

default: $(doc)

slspack.pdf: slspack.dvi
	$(DVIPDF) slspack.dvi

slspack.dvi: $(TEXFILES)
	$(LATEX) slspack.tex
	$(LATEX) slspack.tex

.PHONY: clean

clean:
	@rm -fv *.aux slspack.dvi *.log *.out *.toc slspack.idx slspack.ind slspack.ilg
