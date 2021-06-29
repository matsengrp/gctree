default:

install:
	pip install -r requirements.txt
	pip install -e .

test:
	# pytest
	gctree test

format:
	black gctree
	docformatter --in-place gctree/*.py

lint:
	pylint **/[^_]*.py && echo "LINTING PASS"

docs:
	make -C docs html

deploy:
	make docs
	git checkout gh-pages
	cp -a docs/_build/html/* .
	git add .
	git commit --amend -av -m "update docs"
	git push -f

.PHONY: install test format lint deploy docs
