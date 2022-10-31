default:

install:
	pip install -r requirements.txt
	pip install -e .

test:
	bash tests/smalltest.sh
	pytest

format:
	black gctree tests
	docformatter --in-place gctree/*.py

lint:
	# stop the build if there are Python syntax errors or undefined names
	flake8 . --count --select=E9,F63,F7,F82 --show-source --statistics
	# exit-zero treats all errors as warnings. The GitHub editor is 127 chars wide
	flake8 . --count --max-complexity=30 --max-line-length=127 --statistics

docs:
	wget -O HS5F_Mutability.csv https://bitbucket.org/kleinstein/shazam/raw/ba4b30fc6791e2cfd5712e9024803c53b136e664/data-raw/HS5F_Mutability.csv
	wget -O HS5F_Substitution.csv https://bitbucket.org/kleinstein/shazam/raw/ba4b30fc6791e2cfd5712e9024803c53b136e664/data-raw/HS5F_Substitution.csv
	rm -f docs/outfile docs/outtree # remove phylip dnapars output
	make -C docs html
	make -C docs html # need it again to get images from quickstart commands

.PHONY: install test format lint deploy docs
