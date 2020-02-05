project: JOREK
project_website: http://jorek.eu
summary: JOREK is a non-linear extended MHD code for toroidal X-point geometries
author: the JOREK team
project_dir: ./
src_dir: ./
output_dir: ./doc
exclude_dir: ./doc
	     ./.git
	     ./models/model001
	     ./models/model199
	     ./models/model302
	     ./models/model305
	     ./models/model306
	     ./models/model333
	     ./models/model400
	     ./models/model500
             ./models/model501
	     ./models/model555
	     ./models/model701
	     ./models/model710
docmark: <
docmark_alt: *
predocmark: >
predocmark_alt: #
extra_filetypes: c //
graph: true
coloured_edges: true

This file serves as the input file for FORD to generate documentation on the code
(use `make docs` to get this)
