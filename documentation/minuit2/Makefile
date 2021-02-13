MDFILES=titlepage.md Minuit2.md

OPTIONS=-s -f markdown+smart --toc --toc-depth=2 --top-level-division=chapter --number-sections

OPTIONSHTML=-H css/github.css --mathjax
OPTIONSEPUB=--mathml --epub-cover-image=figures/cover.png
OPTIONSPDFA4=--variable=geometry:a4paper -V geometry:margin=1.5cm
OPTIONSPDFLETTER=--variable=geometry:letterpaper -V geometry:margin=1.5cm

.PHONY: pdfa4 pdfletter html tex epub

all: pdfa4 pdfletter html tex epub

pdfa4: $(MDFILES)
	pandoc $(OPTIONS) $(OPTIONSPDFA4) $(MDFILES) -o Minuit2.pdf

pdfletter: $(MDFILES)
	pandoc $(OPTIONS) $(OPTIONSPDFLETTER) $(MDFILES) -o Minuit2Letter.pdf

html: $(MDFILES) title.txt
	pandoc $(OPTIONS) $(OPTIONSHTML) -t html title.txt $(MDFILES) -o Minuit2.html

tex: $(MDFILES)
	pandoc $(OPTIONS) -t latex $(MDFILES) -o Minuit2.tex

epub: $(MDFILES) metadata.xml title.txt
	pandoc $(OPTIONS) $(OPTIONSEPUB) --epub-metadata=metadata.xml \
	   title.txt $(MDFILES) -o Minuit2.epub

clean:
	rm -f Minuit2.pdf Minuit2Letter.pdf Minuit2.html \
	   Minuit2.tex Minuit2.epub
