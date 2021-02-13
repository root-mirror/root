OUTDIR=output

MDFILES=ROOTUsersGuide.md Preface.md Introduction.md GettingStarted.md \
	Histograms.md Graphs.md FittingHistograms.md ALittleC++.md Cling.md \
	ObjectOwnership.md Graphics.md FoldersTasks.md InputOutput.md \
	Trees.md MathLibraries.md LinearAlgebra.md AddingaClass.md \
	CollectionClasses.md PhysicsVectors.md Geometry.md PythonRuby.md \
	TutorialsandTests.md ExampleAnalysis.md Networking.md Threads.md \
	PROOF.md WritingGUI.md SignalSlot.md HTMLDoc.md InstallandBuild.md

HTMLFILES=$(addprefix $(OUTDIR)/,$(MDFILES:.md=.html))
PDFFILES=$(addprefix $(OUTDIR)/,$(MDFILES:.md=.pdf))

LATEXA4PREFIX=ROOTUsersGuideA4
LATEXLETTERPREFIX=ROOTUsersGuideLetter

OPTIONS=-s -f markdown+smart --toc --toc-depth=2 --top-level-division=chapter --number-sections

OPTIONSHTML=-H css/github.css --mathjax
OPTIONSEPUB=--mathml
OPTIONSLATEXA4=--template=./template.tex --variable=geometry:a4paper -V geometry:margin=1.5cm
OPTIONSLATEXLETTER=--template=./template.tex --variable=geometry:letterpaper -V geometry:margin=1.5cm
OPTIONSPDF=--variable=geometry:a4paper -V geometry:margin=1.5cm

.PHONY: pdfa4 pdfletter html epub pdf_chapters html_chapters

all: pdfa4 pdfletter html epub pdf_chapters html_chapters

$(LATEXA4PREFIX).tex: $(MDFILES)
	pandoc $(OPTIONS) $(OPTIONSLATEXA4) -t latex $(MDFILES) \
	-o $(LATEXA4PREFIX).tex

# on OSX pdflatex returns 1 even in case of success !?
$(LATEXA4PREFIX).toc: $(LATEXA4PREFIX).tex
	pdflatex -interaction=batchmode $(LATEXA4PREFIX).tex || true

$(LATEXA4PREFIX).ind: $(LATEXA4PREFIX).toc ROOTUsersGuide.ist
	makeindex -s ROOTUsersGuide.ist $(LATEXA4PREFIX).idx

pdfa4: $(LATEXA4PREFIX).tex $(LATEXA4PREFIX).ind
	@mkdir -p $(OUTDIR)
	pdflatex -interaction=batchmode -output-directory=$(OUTDIR) $(LATEXA4PREFIX).tex || true

$(LATEXLETTERPREFIX).tex: $(MDFILES)
	pandoc $(OPTIONS) $(OPTIONSLATEXLETTER) -t latex $(MDFILES) \
	-o $(LATEXLETTERPREFIX).tex

$(LATEXLETTERPREFIX).toc: $(LATEXLETTERPREFIX).tex
	pdflatex -interaction=batchmode $(LATEXLETTERPREFIX).tex || true

$(LATEXLETTERPREFIX).ind: $(LATEXLETTERPREFIX).toc ROOTUsersGuide.ist
	makeindex -s ROOTUsersGuide.ist $(LATEXLETTERPREFIX).idx

pdfletter: $(MDFILES) $(LATEXLETTERPREFIX).toc $(LATEXLETTERPREFIX).ind
	@mkdir -p $(OUTDIR)
	pdflatex -interaction=batchmode -output-directory=$(OUTDIR) $(LATEXLETTERPREFIX).tex || true

html: $(MDFILES)
	@mkdir -p $(OUTDIR)
	@mkdir -p $(OUTDIR)/pictures
	@cp -r pictures/* $(OUTDIR)/pictures
	pandoc --quiet $(OPTIONS) $(OPTIONSHTML) -t html $(MDFILES) \
	-o $(OUTDIR)/ROOTUsersGuide.html

epub: $(MDFILES) metadata.xml title.txt
	@mkdir -p $(OUTDIR)
	pandoc $(OPTIONS) $(OPTIONSEPUB) --epub-metadata=metadata.xml \
	title.txt $(MDFILES) -o $(OUTDIR)/ROOTUsersGuide.epub

# --number-offset unfortunately not supported in PDF output :(
$(OUTDIR)/%.pdf: %.md
	@mkdir -p $(OUTDIR)
	@echo Generating PDF chapter: $@
	@pandoc $(OPTIONS) $(OPTIONSPDF) $^ -o $@

pdf_chapters: $(PDFFILES)

$(OUTDIR)/%.html: %.md
	@mkdir -p $(OUTDIR)
	@echo Generating HTML chapter: $@
	@echo % Chapter: $^ | sed -e "s/.md//" | sed -e "s/and/ and /" > $@.txt
	@( OFFSET=$$( for F in $(MDFILES) ; do echo $$F ; done | \
	   grep -n "^$^$$" | head -n1 | cut -d: -f1 ) && \
	   OFFSET=`expr $$OFFSET - 1` && \
	   pandoc $(OPTIONS) $(OPTIONSHTML) -t html $@.txt --number-offset $$OFFSET $^ -o $@ )

html_chapters: $(HTMLFILES)

clean:
	rm -rf $(LATEXA4PREFIX).pdf $(LATEXLETTERPREFIX).pdf \
	ROOTUsersGuide.html ROOTUsersGuide.epub \
	$(LATEXA4PREFIX).{tex,dvi,ind,toc,out,log,idx,aux,ilg} \
	$(LATEXLETTERPREFIX).{tex,dvi,ind,toc,out,log,idx,aux,ilg} \
	$(PDFFILES) $(HTMLFILES) $(OUTDIR)
