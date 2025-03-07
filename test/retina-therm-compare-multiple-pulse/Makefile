SHELL := /bin/zsh

all: plot

output/: output/reedbed/

output/reedbed/:
	mkdir -p output/reedbed/
	powerconf render CONFIG.yml CONFIG-reedbed.yml /dev/stdout | time cargo run --release -- start -f json	-t json > output/reedbed/Tvst.txt 2> reedbed.stderr

.PHONY: check
check:
	@which powerconf >/dev/null || echo "powerconf cli not found. Please install (e.g. pipx install powerconf) and/or make sure the directory with the executable is in your PATH"
	@which awk >/dev/null || echo "awk not found. Please install."
	@which sed >/dev/null || echo "sed not found. Please install."

.PHONY: run
run: output/

.PHONY: clean
clean:
	rm -rf output *.stderr

.PHONY: plot
plot: run
	cd output && gnuplot --persist -e 'plot "reedbed/Tvst.txt" using 1:2 with line title "rise"'
