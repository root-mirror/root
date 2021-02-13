import argparse

def get_argparse():
	DESCRIPTION = """hist2workspace is a utility to create RooFit/RooStats workspace from histograms
"""
	parser = argparse.ArgumentParser(prog='hist2workspace',
	description = DESCRIPTION)
	parser.add_argument("-standard_form", help="""default  model,  which  creates  an  extended PDF that interpolates between RooHistFuncs
This is much faster for models with many bins and uses significantly less memory""")
	parser.add_argument("-number_counting_form", help="""This was the original model in 5.28 (without patches). It uses a Poisson for each
bin of the histogram.  This can become slow and memory intensive when there are many bins.
""")
	parser.add_argument("-v", help="""Switch HistFactory message stream to INFO level.""", action='store_true')
	parser.add_argument("-vv", help="""Switch HistFactory message stream to DEBUG level.""", action='store_true')
	return parser
