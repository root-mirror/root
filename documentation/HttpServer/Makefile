MDFILES=titlepage.md HttpServer.md

OPTIONS=-s -f markdown+smart --toc --toc-depth=2 --top-level-division=chapter --number-sections

OPTIONSHTML=-H css/github.css --mathjax
OPTIONSEPUB=--mathml --epub-cover-image=figures/cover.png
OPTIONSPDFA4=--variable=geometry:a4paper -V geometry:margin=1.5cm
OPTIONSPDFLETTER=--variable=geometry:letterpaper -V geometry:margin=1.5cm

.PHONY: pdfa4 pdfletter html tex epub

all: pdfa4 pdfletter html tex epub

pdfa4: $(MDFILES)
	pandoc $(OPTIONS) $(OPTIONSPDFA4) $(MDFILES) -o HttpServer.pdf

pdfletter: $(MDFILES)
	pandoc $(OPTIONS) $(OPTIONSPDFLETTER) $(MDFILES) -o HttpServerLetter.pdf

html: $(MDFILES) title.txt
	pandoc $(OPTIONS) $(OPTIONSHTML) -t html title.txt $(MDFILES) -o HttpServer.html

tex: $(MDFILES)
	pandoc $(OPTIONS) -t latex $(MDFILES) -o HttpServer.tex

epub: $(MDFILES) metadata.xml title.txt
	pandoc $(OPTIONS) $(OPTIONSEPUB) --epub-metadata=metadata.xml \
	   title.txt $(MDFILES) -o HttpServer.epub

clean:
	rm -f HttpServer.pdf HttpServerLetter.pdf HttpServer.html \
	   HttpServer.tex HttpServer.epub
