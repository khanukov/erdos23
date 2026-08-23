.PHONY: check artifacts replay paper

check: artifacts replay

artifacts:
	./scripts/check_artifacts.sh

replay:
	./scripts/run_full_replay.sh

paper:
	mkdir -p output/pdf
	pdflatex -interaction=nonstopmode -halt-on-error -output-directory=output/pdf paper/erdos23_full_solution_draft.tex
	pdflatex -interaction=nonstopmode -halt-on-error -output-directory=output/pdf paper/erdos23_full_solution_draft.tex
