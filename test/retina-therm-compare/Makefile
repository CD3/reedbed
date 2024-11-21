SHELL := /bin/zsh

all: plot

output/: output/reedbed/ output/retina-therm/

output/reedbed/:
	mkdir -p output/reedbed/
	powerconf render CONFIG.yml CONFIG-reedbed.yml /dev/stdout | time cargo run --release -- temperature-rise -f json	-t json > output/reedbed/Tvst.txt 2> reedbed.stderr

output/retina-therm/:
	mkdir -p output/retina-therm/
	powerconf render CONFIG.yml CONFIG-retina-therm.yml /dev/stdout | time retina-therm temperature-rise /dev/stdin 2> retina-therm.stderr

.PHONY: check
check:
	@which powerconf >/dev/null || echo "powerconf cli not found. Please install (e.g. pipx install powerconf) and/or make sure the directory with the executable is in your PATH"
	@which retina-therm >/dev/null || echo "retina-therm cli not found. Please install (e.g. pipx install retina-therm) and/or make sure the directory with the executable is in your PATH"
	@which awk >/dev/null || echo "awk not found. Please install."
	@which sed >/dev/null || echo "sed not found. Please install."

.PHONY: run
run: output/

.PHONY: clean
clean:
	rm -rf output \*.stderr

.PHONY: plot
plot: run
	cd output && gnuplot --persist -e 'plot "reedbed/Tvst.txt", "retina-therm/Tvst.txt"'
