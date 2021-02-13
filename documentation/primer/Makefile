MDFILES=abstract.md Introduction.md ROOT_as_calculator.md \
	your_first_ROOT_macro.md graphs.md histograms.md \
	functions_and_parameter_estimation.md filio.md root_in_python.md \
	concludingRemarks.md references.md

MDFILESMOD=abstract_preprocessed.md Introduction_preprocessed.md ROOT_as_calculator_preprocessed.md \
	   your_first_ROOT_macro_preprocessed.md graphs_preprocessed.md histograms_preprocessed.md \
	   functions_and_parameter_estimation_preprocessed.md filio_preprocessed.md root_in_python_preprocessed.md \
	   concludingRemarks_preprocessed.md references_preprocessed.md

OPTIONS=-s -f markdown+smart --toc --toc-depth=2 --top-level-division=chapter --number-sections --bibliography=bibliography.bib

OPTIONSHTML=-H css/github.css --mathjax
OPTIONSEPUB=--mathml --epub-cover-image=figures/Cover.png
OPTIONSPDFA4=--variable=geometry:a4paper -V geometry:margin=1.5cm --template template.tex
OPTIONSPDFLETTER=--variable=geometry:letterpaper -V geometry:margin=1.5cm --template template.tex

.PHONY: pdfa4 pdfletter html tex epub

all: pdfa4 pdfletter html tex epub

%_preprocessed.md: %.md
	./rootDocPreprocessor.py $<

pdfa4: $(MDFILES) bibliography.bib $(MDFILESMOD)
	pandoc $(OPTIONS) $(OPTIONSPDFA4) $(MDFILESMOD) -o ROOTPrimer.pdf

pdfletter: $(MDFILES) bibliography.bib $(MDFILESMOD)
	pandoc $(OPTIONS) $(OPTIONSPDFLETTER) $(MDFILESMOD) -o ROOTPrimerLetter.pdf

html: $(MDFILES) bibliography.bib $(MDFILESMOD) title.txt
	pandoc $(OPTIONS) $(OPTIONSHTML) -t html title.txt $(MDFILESMOD) -o ROOTPrimer.html

tex: $(MDFILES) bibliography.bib $(MDFILESMOD)
	pandoc $(OPTIONS) -t latex $(MDFILESMOD) -o ROOTPrimer.tex

epub: $(MDFILES) metadata.xml title.txt bibliography.bib $(MDFILESMOD)
	pandoc $(OPTIONS) $(OPTIONSEPUB) --epub-metadata=metadata.xml \
	   title.txt $(MDFILESMOD) -o ROOTPrimer.epub

clean:
	rm -f ROOTPrimer.pdf ROOTPrimerLetter.pdf ROOTPrimer.html \
	   ROOTPrimer.tex ROOTPrimer.epub $(MDFILESMOD)
